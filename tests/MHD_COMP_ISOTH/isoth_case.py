#!/usr/bin/env python3
"""Generate calliope.in for MHD_COMP_ISOTH test cases.

    from isoth_case import write_case
    write_case('runs/stage1/sound', nproc=1,
               box=dict(nlx=64), time=dict(dt=1e-3, nstep=1000),
               ic=('sound_wave', dict(dir="'x'", amp=1e-3)))

Each keyword (box, time, phys, num, shock, diag, isoth_diag, ...) overrides
entries of the corresponding namelist below. Values are written verbatim,
so strings must carry their own quotes. Python bools become T/F.
"""
import os
import copy

DEFAULTS = {
    'mpi_settings':           dict(dim_decomp=2),
    'box_parameters':         dict(lx=1.0, ly=1.0, lz=1.0, nlx=64, nly=2, nlz=2, pruned=False),
    'time_parameters':        dict(dt=1e-3, nstep=1000, cfl=0.5, reset_method="'stop'", increase_dt=1e10),
    'scheme_parameters':      dict(time_step_scheme="'eSSPIFRK3'"),
    'operation_parameters':   dict(nonlinear=True),
    'physical_parameters':    dict(nu=0.0, nu_h=0.0, nu_h_exp=4, eta=0.0, eta_h=0.0, eta_h_exp=4,
                                   lmd=0.0, lmd_h=0.0, lmd_h_exp=4, beta0=2.0, rho_min=1e-2,
                                   shear=False, q=0.0),
    'numerics_parameters':    dict(formulation="'C'", visc_form="'flux'", filter_lnrho=False,
                                   filter_mom=True, filter_mag=True, cfl_visc=1.0,
                                   mass_correction=False, s_abort=100.0, debug_nan_step=-1),
    'shock_parameters':       dict(C_shock=2.0, shock_nmax=1, shock_nsmooth=1, shock_dl=0.0,
                                   shock_selftest=False),
    'svv_parameters':         dict(svv_eps0=0.5, svv_mN=0.5, svv_above_kc=True, svv_form="'entropy'",
                                   cfl_svv=2.0, svv_pi_init=50, svv_pi_check=100),
    'isoth_diagnostics':      dict(write_intvl_1D=0.0, write_intvl_cons=0.0, budget=True),
    'force':                  dict(driven=False),
    # below the 30 min of A_dev, so that a long run ends by itself (and is flagged by analyze_stage1)
    'terminate':              dict(max_wall_time=0.45, terminate_file_name="'terminate'"),
    'diagnostics_parameters': dict(write_intvl=1e10, write_intvl_2D=0.0, write_intvl_3D=0.0),
    'restart_parameters':     dict(save_restart_intvl=1e10, restart_dir="'./restart/'"),
}

# short keyword -> namelist
ALIAS = {
    'mpi': 'mpi_settings', 'box': 'box_parameters', 'time': 'time_parameters',
    'scheme': 'scheme_parameters', 'op': 'operation_parameters', 'phys': 'physical_parameters',
    'num': 'numerics_parameters', 'shock': 'shock_parameters', 'svv': 'svv_parameters',
    'isoth_diag': 'isoth_diagnostics',
    'force': 'force', 'term': 'terminate', 'diag': 'diagnostics_parameters',
    'restart': 'restart_parameters',
}

# init_type -> namelist of its parameters
IC_NML = {'sound_wave': 'wave1d_params', 'simple_wave': 'wave1d_params',
          'riemann': 'riemann_params', 'single_mode': 'initial_condition_params',
          'random': 'initial_condition_params', 'profile1d': 'profile1d_params', 'zero': None}


def fmt(v):
    if isinstance(v, bool):
        return 'T' if v else 'F'
    if isinstance(v, float):
        return repr(v).replace('e', 'd') if 'e' in repr(v) else repr(v) + 'd0'
    if isinstance(v, (list, tuple)):
        return ' '.join(fmt(x) for x in v)
    return str(v)


def render(nml):
    lines = []
    for name, entries in nml.items():
        lines.append('&' + name)
        for k, v in entries.items():
            lines.append('  %s = %s' % (k, fmt(v)))
        lines.append('/')
    return '\n'.join(lines) + '\n'


def make_nml(ic=('zero', {}), **kw):
    nml = copy.deepcopy(DEFAULTS)
    for key, over in kw.items():
        nml[ALIAS.get(key, key)].update(over)
    init_type, params = ic
    nml['initial_condition'] = dict(init_type="'%s'" % init_type)
    if IC_NML[init_type] is not None:
        nml[IC_NML[init_type]] = dict(params)
    return nml


def write_case(path, nproc=1, header='', **kw):
    os.makedirs(path, exist_ok=True)
    with open(os.path.join(path, 'calliope.in'), 'w') as f:
        if header:
            f.write(''.join('! %s\n' % l for l in header.splitlines()))
        f.write(render(make_nml(**kw)))
    with open(os.path.join(path, 'nproc'), 'w') as f:
        f.write('%d\n' % nproc)
    return path
