#!/usr/bin/env python3
"""Deterministic Gallope <-> Calliope comparison for MHD_INCOMP.

The two codes store the same six spectral fields, but on different half
spaces, because the real transform is taken along a different axis:

  Gallope   real transform along z.  Global Fortran array (nkz, nky, nkx)
            with nkz = nlz/2+1, nky = nly, nkx = nlx.
  Calliope  real transform along y.  Global Fortran array (nkx, nkz, nky)
            with nkx = nlx, nkz = nlz, nky = nly/2+1.

Neither is a transpose of the other, so this script reconstructs the full
nlx x nly x nlz complex cube from each code by Hermitian expansion and
compares mode by mode.  No FFT and no normalization convention is involved,
so the comparison is exact.

Self-check: the reconstruction is only correct if the axis mapping is right.
Both initial conditions are exactly solenoidal, and the equations preserve
that, so <|k.u|>/(<|k|><|u|>) computed from the reconstructed cube must come
out at roundoff.  A mis-mapped axis shows up immediately as an O(1) value.
This is reported for every field set and should be inspected first.

In a shearing box the constraint is imposed on the sheared wavenumber
kx(t) = kx + q*tsc*ky, not on kx (shearing_box.F90:63 in Gallope, :62 in
Calliope, identical in both), tsc being the time since the last remap as
recorded in restart/time.dat.  Pass --q to have the check use it; with the
static kx the same field sets come out at 1e-8 rather than 1e-24, which is
small enough to look healthy and is nevertheless the wrong quantity.

Usage:
    compare.py --a DIR --a-code gallope --b DIR --b-code calliope \\
               --nlx 128 --nly 128 --nlz 128 [--lx 1.0 --ly 1.0 --lz 1.0]
               [--q 1.5]
"""

import argparse
import os
import sys

import numpy as np

FIELDS = ("ux", "uy", "uz", "bx", "by", "bz")
ITEMSIZE = 16  # complex(8)


def code_wavenumbers(n):
    """Signed wavenumber index, matching the Fortran convention in both codes:
    ikx(i) = i-1 for i <= n/2+1, else i-n-1  (i is 1-based)."""
    idx = np.arange(n)
    return np.where(idx <= n // 2, idx, idx - n)


def read_raw(path, shape_fortran):
    """Read an MPI_ORDER_FORTRAN complex(8) dump and return it with the
    Fortran index order preserved as a C-indexed array."""
    expected = ITEMSIZE * int(np.prod(shape_fortran))
    actual = os.path.getsize(path)
    if actual != expected:
        raise SystemExit(
            f"{path}: size {actual} B but the declared grid needs {expected} B "
            f"(shape {shape_fortran}, complex(8)). Check nlx/nly/nlz, and for "
            f"Calliope check that the run used `pruned = .false.`."
        )
    flat = np.fromfile(path, dtype=np.complex128)
    return flat.reshape(shape_fortran, order="F")


def load_gallope(restart_dir, nlx, nly, nlz):
    """-> dict of arrays indexed [ikx, iky, ikz] with ikz = 0 .. nlz/2."""
    nkz = nlz // 2 + 1
    out = {}
    for f in FIELDS:
        a = read_raw(os.path.join(restart_dir, f + ".dat"), (nkz, nly, nlx))
        out[f] = np.ascontiguousarray(a.transpose(2, 1, 0))  # -> [ikx, iky, ikz]
    return out


def load_calliope(restart_dir, nlx, nly, nlz):
    """-> dict of arrays indexed [ikx, iky, ikz] with iky = 0 .. nly/2."""
    nky = nly // 2 + 1
    out = {}
    for f in FIELDS:
        a = read_raw(os.path.join(restart_dir, f + ".dat"), (nlx, nlz, nky))
        out[f] = np.ascontiguousarray(a.transpose(0, 2, 1))  # -> [ikx, iky, ikz]
    return out


def expand_hermitian(half, axis, nlx, nly, nlz):
    """Fill the full (nlx, nly, nlz) cube from a half space that is complete
    on every axis but `axis`, where it holds indices 0 .. n/2.

    Uses F[k] = conj(F[-k]), which holds because the physical fields are real.
    """
    full = np.zeros((nlx, nly, nlz), dtype=np.complex128)
    n = (nlx, nly, nlz)[axis]
    nhalf = n // 2 + 1

    sl = [slice(None)] * 3
    sl[axis] = slice(0, nhalf)
    full[tuple(sl)] = half

    # Mirror indices: the partner of full index m is (-m) % n on every axis.
    mx = (-np.arange(nlx)) % nlx
    my = (-np.arange(nly)) % nly
    mz = (-np.arange(nlz)) % nlz
    src = [mx, my, mz]
    # On `axis`, only indices nhalf .. n-1 are missing; their partners
    # (-m) % n land in 1 .. n/2, which the half space does hold.
    src[axis] = (-np.arange(nhalf, n)) % n

    sl[axis] = slice(nhalf, n)
    full[tuple(sl)] = np.conj(half[np.ix_(*src)])
    return full


def full_cube(fields, code, nlx, nly, nlz):
    axis = 2 if code == "gallope" else 1  # halved axis: z for Gallope, y for Calliope
    return {f: expand_hermitian(fields[f], axis, nlx, nly, nlz) for f in FIELDS}


def kgrids(nlx, nly, nlz, lx, ly, lz):
    """Physical wavenumbers. lx/ly/lz are the input values, i.e. box = 2*pi*lx."""
    kx = code_wavenumbers(nlx) / lx
    ky = code_wavenumbers(nly) / ly
    kz = code_wavenumbers(nlz) / lz
    return (kx[:, None, None], ky[None, :, None], kz[None, None, :])


def div_measure(cube, prefix, kx, ky, kz):
    """<|k.f|>/(<|k|><|f|>), the quantity is_div_free prints, but evaluated on
    the full cube so that it is identical for the two codes."""
    fx, fy, fz = (cube[prefix + c] for c in "xyz")
    div = np.abs(kx * fx + ky * fy + kz * fz).sum()
    amp = np.sqrt(np.abs(fx) ** 2 + np.abs(fy) ** 2 + np.abs(fz) ** 2).sum()
    kabs = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2).sum()
    if amp == 0.0 or kabs == 0.0:
        return float("nan")
    return div / (amp * kabs)


def energy(cube):
    """W = 0.5 <|u|^2 + |b|^2>. With field_k = FFT(field_r)/N, Parseval gives
    <field_r^2> = sum_k |field_k|^2 over the full cube."""
    tot = 0.0
    for f in FIELDS:
        tot += float((np.abs(cube[f]) ** 2).sum())
    return 0.5 * tot


def read_restart_time(restart_dir):
    path = os.path.join(restart_dir, "time.dat")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        lines = [ln for ln in fh.read().split("\n") if ln.strip()]
    if len(lines) < 2:
        return None
    return [float(v) for v in lines[-1].split()]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--a", required=True, help="run directory A (contains restart/)")
    p.add_argument("--b", required=True, help="run directory B (contains restart/)")
    p.add_argument("--a-code", required=True, choices=("gallope", "calliope"))
    p.add_argument("--b-code", required=True, choices=("gallope", "calliope"))
    p.add_argument("--nlx", type=int, default=128)
    p.add_argument("--nly", type=int, default=128)
    p.add_argument("--nlz", type=int, default=128)
    p.add_argument("--lx", type=float, default=1.0)
    p.add_argument("--ly", type=float, default=1.0)
    p.add_argument("--lz", type=float, default=1.0)
    p.add_argument("--q", type=float, default=0.0,
                   help="shear rate; when nonzero the solenoidal check uses "
                        "kx + q*tsc*ky with tsc taken from restart/time.dat")
    args = p.parse_args()

    n = (args.nlx, args.nly, args.nlz)
    loaders = {"gallope": load_gallope, "calliope": load_calliope}

    cubes, times = {}, {}
    for tag, d, code in (("A", args.a, args.a_code), ("B", args.b, args.b_code)):
        rd = os.path.join(d, "restart")
        if not os.path.isdir(rd):
            raise SystemExit(f"{rd} does not exist")
        cubes[tag] = full_cube(loaders[code](rd, *n), code, *n)
        times[tag] = read_restart_time(rd)

    kx, ky, kz = kgrids(*n, args.lx, args.ly, args.lz)

    print(f"grid {args.nlx} x {args.nly} x {args.nlz}, "
          f"box (2pi*{args.lx}, 2pi*{args.ly}, 2pi*{args.lz})")
    print(f"A = {args.a}  [{args.a_code}]")
    print(f"B = {args.b}  [{args.b_code}]")
    print()

    # --- gate 1: the two runs must be at the same time -----------------------
    print("-- restart time (tt, tsc) --")
    for tag in ("A", "B"):
        print(f"  {tag}: {times[tag]}")
    if times["A"] and times["B"]:
        dtt = abs(times["A"][0] - times["B"][0])
        rel = dtt / abs(times["A"][0]) if times["A"][0] else dtt
        verdict = "OK" if rel < 1e-12 else "MISMATCH -- dt diverged, comparison invalid"
        print(f"  |tt_A - tt_B| = {dtt:.3e}  ({verdict})")
    print()

    # --- gate 2: axis mapping self-check ------------------------------------
    print("-- solenoidal error on the reconstructed cube --")
    print("   (must be at roundoff; an O(1) value means the axis mapping is wrong)")
    for tag in ("A", "B"):
        # The constraint is on the sheared wavenumber, so each run is checked
        # against its own tsc rather than against a common grid.
        tsc = times[tag][1] if (args.q and times[tag] and len(times[tag]) > 1) else 0.0
        kxt = kx + args.q * tsc * ky
        du = div_measure(cubes[tag], "u", kxt, ky, kz)
        db = div_measure(cubes[tag], "b", kxt, ky, kz)
        note = f"  [kx(t) = kx + {args.q}*{tsc:.6g}*ky]" if args.q else ""
        print(f"  {tag}: <|k.u|>/(<|k|><|u|>) = {du:.3e}   "
              f"<|k.b|>/(<|k|><|b|>) = {db:.3e}{note}")
    print()

    # --- the comparison -----------------------------------------------------
    print("-- relative field difference, ||f_A - f_B||_2 / ||f_B||_2 --")
    num2 = den2 = 0.0
    for f in FIELDS:
        a, b = cubes["A"][f], cubes["B"][f]
        dn = float(np.linalg.norm((a - b).ravel()))
        bn = float(np.linalg.norm(b.ravel()))
        num2 += dn ** 2
        den2 += bn ** 2
        rel = dn / bn if bn else float("nan")
        amax = float(np.abs(a - b).max())
        bmax = float(np.abs(b).max())
        print(f"  {f}: L2 {rel:.6e}   max|df|/max|f| "
              f"{(amax / bmax if bmax else float('nan')):.6e}")
    print(f"  all six fields combined: {np.sqrt(num2 / den2):.6e}")
    print()

    print("-- energy W = 0.5<|u|^2+|b|^2> from the spectral fields --")
    wa, wb = energy(cubes["A"]), energy(cubes["B"])
    print(f"  A: {wa:.15e}")
    print(f"  B: {wb:.15e}")
    print(f"  relative difference: {abs(wa - wb) / wb:.6e}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
