#!/usr/bin/env python3
"""Inputs of the figure scripts <Test>/plot_<test>.py that exist only on ps1, and the command line the scripts share.

    python3 tests/MHD_COMP_ISOTH/<Test>/plot_<test>.py [plot] [name ...]   figures -> <Test>/fig/<name>.pdf
    python3 tests/MHD_COMP_ISOTH/<Test>/plot_<test>.py fetch               missing inputs from ps1 into runs/
    python3 tests/MHD_COMP_ISOTH/<Test>/plot_<test>.py list                figures and their inputs
    options: --runs DIR (default: tests/MHD_COMP_ISOTH/runs/, from any directory), --host HOST (ps1),
             --remote DIR (~/work/calliope_isoth)

fetch writes nothing on ps1 and never overwrites a local file (an existing target is skipped). It
    copies small files (summaries, out1d lines) with tar over ssh to the same path under the local runs/;
    extracts arrays from large outputs by running this file on ps1,
        ssh ps1 'cd ~/work/calliope_isoth && python3 -B - <md5> <op> <arg> ...' < figdata.py > runs/figdata/...npz,
    which imports the analysis modules of the worktree there (tests/MHD_COMP_ISOTH/stages), stops unless their md5
    equal those of the local stages/ (MODS), and writes an npz to stdout. The ops run in REMOTE_RUNS there, so a
    <case> is relative to the runs/ of the tests like the copied files. The ops:
        plane <case> <y0|z0> <t> <nlx,nly,nlz> <q,...>   quantities q of gen_ic.PLANE_Q (float32) of the plane
                                                          record nearest to t (analyze_stage4.plane_records)
        table <case> <file>                               a table file of the run (cfl.dat, ...) per output
                                                          directory (analyze_stage1.load_runs) and its header
        planar <case>                                     the last line of analyze_stage4.planar_snaps
Python 3.9 (ps1) and later.
"""
import os
import sys
import io
import json
import hashlib
import shlex
import subprocess
import tarfile
import numpy as np

MODS = ('analyze_stage1.py', 'analyze_stage4.py', 'gen_ic.py', 'gen_stage4.py', 'riemann_exact.py',
        'cp_stability.py')
REMOTE_MOD_DIR = 'tests/MHD_COMP_ISOTH/stages'
REMOTE_RUNS = 'tests/MHD_COMP_ISOTH/runs'


def md5(path):
    with open(path, 'rb') as f:
        return hashlib.md5(f.read()).hexdigest()


# --------------------------------------------------------------------------------------- remote side (ps1)
def op_plane(case, plane, t, dims, qs):
    import analyze_stage4 as a4
    import gen_ic
    recs = a4.plane_records(case, plane)
    if not recs:
        raise SystemExit('no plane records in %s' % case)
    tr, step, path, rec = min(recs, key=lambda r: abs(r[0] - float(t)))
    f = a4.plane_at(gen_ic.Grid(*[int(n) for n in dims.split(',')]), path, rec, plane)
    out = dict(t=np.float64(tr), step=np.int64(step), times=np.array([r[0] for r in recs]))
    out.update({q: f[q].astype(np.float32) for q in qs.split(',')})
    return out


def op_table(case, name):
    import analyze_stage1 as a1
    parts = a1.load_runs(case, name)
    if not parts:
        raise SystemExit('no %s in %s' % (name, case))
    hdr = ''
    for d in a1.run_dirs(case):
        fn = os.path.join(d, name)
        if os.path.exists(fn) and os.path.getsize(fn) > 0:
            with open(fn) as f:
                hdr = f.readline().strip()
            break
    out = {'part%d' % i: a for i, a in enumerate(parts)}
    out.update(n=np.int64(len(parts)), header=np.array(hdr))
    return out


def op_planar(case):
    import analyze_stage4 as a4
    m = a4.case_info(case)
    snaps, steps, dbn = a4.planar_snaps(case, m)
    if not snaps:
        raise SystemExit('no lines in %s' % case)
    t, W, zeta, eta, D = snaps[-1]
    info = {k: m[k] for k in ('prob', 'dim', 'kvec', 'kn', 'nrun', 'dt')}
    info['N'] = m['case']['N']
    return dict(t=np.float64(t), W=W, zeta=zeta, eta=eta, D=D, times=np.array([s[0] for s in snaps]),
                dBn=np.float64(dbn), info=np.array(json.dumps(info)))


OPS = dict(plane=op_plane, table=op_table, planar=op_planar)


def remote(argv):
    """python3 - <name:md5,...> <op> <arg> ... in the worktree on the remote host: the npz to stdout."""
    out = os.fdopen(os.dup(1), 'wb')
    os.dup2(2, 1)                                # anything printed goes to stderr, not into the npz
    mod = REMOTE_MOD_DIR
    if not os.path.exists(os.path.join(mod, 'analyze_stage4.py')):
        raise SystemExit('analysis modules not found in %s' % os.path.join(os.getcwd(), mod))
    want = dict(s.split(':') for s in argv[0].split(','))
    bad = [n for n, h in want.items() if md5(os.path.join(mod, n)) != h]
    if bad:
        raise SystemExit('md5 of %s in %s differ from the local stages/' % (', '.join(bad), mod))
    sys.path.insert(0, os.path.abspath(mod))
    os.chdir(REMOTE_RUNS)
    buf = io.BytesIO()
    np.savez(buf, **OPS[argv[1]](*argv[2:]))
    out.write(buf.getvalue())
    out.flush()


# --------------------------------------------------------------------------------------------- local side
FILE = globals().get('__file__', '<stdin>')       # '<stdin>' under python3 - on the remote host
HERE = None if FILE == '<stdin>' else os.path.dirname(os.path.abspath(FILE))


class Copy:
    """Files copied as they are: paths relative to runs/."""

    def __init__(self, *paths):
        self.paths = paths

    def targets(self):
        return list(self.paths)


class Extract:
    """An npz from an op run on the remote host: out relative to the local runs/, case to REMOTE_RUNS."""

    def __init__(self, out, op, *args):
        self.out, self.op, self.args = out, op, args

    def targets(self):
        return [self.out]


def ssh_cmd(opt, cmd):
    return ['ssh', '-o', 'RemoteCommand=none', '-o', 'RequestTTY=no', '-o', 'BatchMode=yes', '-o',
            'ConnectTimeout=20', opt.host, cmd]


def remote_dir(opt):
    r = opt.remote
    return r if not r.startswith('~/') else '"$HOME"/' + shlex.quote(r[2:])


def fetch_copy(opt, item):
    todo = [p for p in item.paths if not os.path.exists(os.path.join(opt.runs, p))]
    for p in item.paths:
        if p not in todo:
            print('exists, skipped:', os.path.join(opt.runs, p))
    if not todo:
        return
    cmd = 'cd %s/%s && tar cf - -- %s' % (remote_dir(opt), REMOTE_RUNS, ' '.join(shlex.quote(p) for p in todo))
    proc = subprocess.Popen(ssh_cmd(opt, cmd), stdout=subprocess.PIPE)
    want, got = set(todo), set()
    with tarfile.open(fileobj=proc.stdout, mode='r|') as tf:
        for m in tf:
            if m.name not in want or not m.isfile():
                continue                         # only the requested regular files, at their own paths
            dst = os.path.join(opt.runs, m.name)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            tmp = dst + '.part'
            with open(tmp, 'wb') as f:
                f.write(tf.extractfile(m).read())
            os.replace(tmp, dst)
            got.add(m.name)
            print('copied %s (%d bytes)' % (dst, m.size))
    if proc.wait() != 0 or got != want:
        raise SystemExit('copy from %s failed; missing: %s' % (opt.host, ', '.join(sorted(want - got))))


def fetch_extract(opt, item):
    dst = os.path.join(opt.runs, item.out)
    if os.path.exists(dst):
        print('exists, skipped:', dst)
        return
    sums = ','.join('%s:%s' % (n, md5(os.path.join(HERE, n))) for n in MODS)
    cmd = 'cd %s && python3 -B - %s' % (remote_dir(opt), ' '.join(shlex.quote(a) for a in (sums, item.op) + item.args))
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    tmp = dst + '.part'
    with open(FILE, 'rb') as src, open(tmp, 'wb') as f:
        r = subprocess.run(ssh_cmd(opt, cmd), stdin=src, stdout=f)
    if r.returncode != 0:
        os.remove(tmp)
        raise SystemExit('%s %s failed on %s' % (item.op, ' '.join(item.args), opt.host))
    with np.load(tmp) as d:                      # a complete npz
        d.files
    os.replace(tmp, dst)
    print('extracted %s (%d bytes)' % (dst, os.path.getsize(dst)))


def main(test_dir, figs, fetches, argv=None):
    """Command line of a plot_<test>.py: figs {name: f(runs)} drawn in order; fetches [Copy | Extract]. A figure
    whose input is missing (FileNotFoundError) is skipped and listed at the end."""
    import argparse
    import paths
    cmds = ('plot', 'fetch', 'list')
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] not in cmds and not argv[0].startswith('-'):
        argv.insert(0, 'plot')                   # names alone draw those figures
    ap = argparse.ArgumentParser(description='figures of %s' % os.path.basename(test_dir))
    ap.add_argument('cmd', nargs='?', default='plot', choices=cmds)
    ap.add_argument('names', nargs='*', help='figures to draw (default: all)')
    ap.add_argument('--runs', default=paths.RUNS)
    ap.add_argument('--host', default='ps1')
    ap.add_argument('--remote', default='~/work/calliope_isoth')
    opt = ap.parse_args(argv)
    if opt.cmd == 'list':
        for name, f in figs.items():
            print('%-18s %s' % (name, ' '.join((f.__doc__ or '').split()).split('. ')[0]))
        for item in fetches:
            print('fetch: %s' % ', '.join(item.targets()))
        return
    if opt.cmd == 'fetch':
        for item in fetches:
            (fetch_copy if isinstance(item, Copy) else fetch_extract)(opt, item)
        return
    unknown = [n for n in opt.names if n not in figs]
    if unknown:
        raise SystemExit('unknown figures: %s (list shows them)' % ', '.join(unknown))
    skipped = []
    for name in opt.names or list(figs):
        try:
            figs[name](opt.runs)
        except FileNotFoundError as e:
            skipped.append((name, e.filename or str(e)))
            print('skipped %s: missing %s' % (name, e.filename or e))
    if skipped:
        print('skipped%s: %s' % (' (run fetch first)' if fetches else '', ', '.join(n for n, _ in skipped)))


if __name__ == '__main__' and HERE is None:
    remote(sys.argv[1:])
