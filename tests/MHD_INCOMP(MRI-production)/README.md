# MHD_INCOMP(MRI-production) — nonlinear MRI production run, Calliope half

The CPU-side half of the shearing-box MRI comparison behind Fig. 4 and
Section 5 of COMPHY-D-26-00110, rerun after the shearing-box pressure fix
(`92f4047` here, `49e52a2` in Gallope).  The GPU half is in
`gallope_dev/tests/MHD_INCOMP(MRI-production)/`, whose README carries the
physical rationale; this file records only what is specific to Calliope.

## Run definition

`calliope.in` is parameter-for-parameter identical to the Gallope input —
512 x 512 x 256 on a `2*pi * (2, 4, 1)` shearing box, `q = 1.5`,
`b0 = (0, 0, 5e-2)`, random solenoidal seed with `|k_j| <= 3` and amplitude
5e-2, `eSSPIFRK3` at `cfl = 2.0`, and `nu_h = eta_h = 2.00e2` with
`nu_h_exp = eta_h_exp = 2`.  The only intended differences are the
`&mpi_settings` group (`dim_decomp = 2`) and an empty `&forced_fields`, both of
which Gallope does not have.

## Placement

`job.pbs-prod-calliope`: 16 Miyabi-G nodes, `mpiprocs=72:ompthreads=1`,
1152 MPI ranks, queue `regular-g`.  Measured 0.149 s/step at this grid — 1.7x
the cost of Gallope on the same 16 nodes, which is the performance statement
the comparison supports.

Like the Gallope job it is wall-time bounded and archives each finished leg
into `legs/leg<N>/` before resubmitting, because the NetCDF diagnostics are
opened with clobber semantics and no restart branch.  There is no `legs/` here:
the run finished inside a single leg.

## Outcome

Clean exit: 1 127 444 steps to `t = 1181.37`, graceful `max_wall_time`
shutdown, `rc = 0`.  (The Gallope half hung at `t = 920` and was killed by PBS;
its README records the forensics.  Both runs cover the averaging window with
room to spare, so the comparison is unaffected.)

Section 5 quantities, time-averaged over `t` in `[200, 900]` — the same window
as Gallope, which is a change from the first revision
(`diagnostics/MHD_INCOMP/time_average.txt`, indices 395–1789):

| | Calliope | Gallope |
|---|---|---|
| `<W>` | 3.792 +- 0.808 | 4.043 +- 1.150 |
| `<W_dot>` | 1.236e-3 +- 1.334e-1 | 5.055e-3 +- 1.724e-1 |
| `<P>` | 1.945 +- 0.404 | 1.987 +- 0.562 |
| `<D>` | 1.943 +- 0.409 | 1.980 +- 0.574 |

`<W_dot>` vanishes to within its scatter and `<P>` matches `<D>` to 0.1 %, so
the window is in the statistically steady state and the budget closes.  The two
codes agree to 6.6 % in `<W>`, 2.2 % in `<P>` and 1.9 % in `<D>`, all inside
1 sigma; before the pressure fix they disagreed by 33 % in both `<W>` and
`<P>`.

## What is in this directory

| | |
|---|---|
| kept, tracked | `calliope.in`, `job.pbs-prod-calliope`, `out.std`, this file |
| kept, untracked | `rms.dat`, `cfl.dat`, `calliope.out*.nc`, `out2d/`, the figures under `diagnostics/` |
| not copied from Miyabi | `restart/` (930 MB), `calliope.modes.out` (1.2 GB), and all but the last snapshot of `out2d/` |

`out2d/` is **only the final snapshot**, 48 MB out of the 4.2 GB on Miyabi;
the Gallope README explains how it was truncated and why that is all Fig. 4(b)
needs.  `plot_fields.py` reproduces the panel at `t = 1181.37` from it.

`diagnostics/` is a copy of the analysis tree as it stood when the figures were
made, with `time_average.py`'s `avg_start_time`/`avg_end_time` pinned to the
window above.  `load.py` opens `<runname>.out.nc` unconditionally, so the
NetCDF file — untracked, but present here — is what makes the archive
self-contained: `time_average.py` and `plot_fields.py` both run to completion
from this directory alone, with no access to Miyabi.
