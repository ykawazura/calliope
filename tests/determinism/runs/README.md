# runs/ — the Calliope half of the determinism comparison, as executed

Copied verbatim from Miyabi-G (`~/work/determinism/`) before the machine was
decommissioned on 2026-08-13.  The protocol document, the result tables and the
Gallope half are in `gallope_dev/tests/determinism/`; that README is the one to
read first.

## Layout

| directory | what ran | count |
|---|---|---|
| `rc/` | MHD_INCOMP, `run_{shear,noshear}_n{1,10,100,1000}_c16` plus `run_shear_n2000_c16`, 16 MPI ranks | 9 |
| `probe/` | MHD_INCOMP convergence probes: `{probe,order,lin}_dt{1e3,5e4,25e5}_calliope` plus `probe_lin_calliope` | 10 |
| `rmhd/` | RMHD, `run_n{1,10,100,1000}_c16` plus `probe_dt{1e3,5e4,25e5}_calliope` | 7 |

Each run directory holds the input it was started from, its `out.std`, and its
`rms.dat` / `cfl.dat` time series.  The PBS job logs sit one level up.

The comparison outputs (`compare_all_fixed.out`, `probe_all_fixed.out`,
`compare_rmhd.out`) are copied here as well, but they are joint results — they
describe a Gallope/Calliope pair, not Calliope alone.

## What was not copied

`restart/*.dat` — the spectral field dumps the comparison reads — is about
100 MB per run and was left on Miyabi, along with `out2d/`, `out3d/`,
`*.modes.out`, the NetCDF diagnostics and the run binaries.
