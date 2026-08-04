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
opened with clobber semantics and no restart branch.

## What is in this directory

| | |
|---|---|
| kept, tracked | `calliope.in`, `job.pbs-prod-calliope`, `out.std`, this file |
| kept, untracked | `rms.dat`, `cfl.dat`, `calliope.out*.nc`, the figures under `diagnostics/` |
| not copied from Miyabi | `restart/` (930 MB), `out2d/` (4.2 GB), `calliope.modes.out` (1.2 GB) |
