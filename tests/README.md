# tests/ — archived runs from the Gallope/Calliope comparison

Calliope is the CPU code that Gallope (`github.com/ykawazura/gallope`) was
ported from, and the two are validated against each other.  This directory is
the Calliope half of the runs that comparison consists of; each subdirectory
is mirrored by one of the same name under `gallope_dev/tests/`.

| directory | comparison | model |
|---|---|---|
| `determinism/runs/` | deterministic agreement from a shared 3D Orszag-Tang state, 128^3 | MHD_INCOMP and RMHD |
| `MHD_INCOMP(MRI-production)/` | nonlinear MRI production run, 512 x 512 x 256 shearing box | MHD_INCOMP |

Both were run on Miyabi-G (JCAHPC, project gr96) for the revision of
COMPHY-D-26-00110 and were copied here before that machine was decommissioned
on 2026-08-13.

## What is and is not kept

Inputs, job scripts, stdout and PBS logs, comparison results and the analysis
scripts are committed.  The bulk output is not: `rms.dat` and `cfl.dat` are
present in the working tree but ignored, and `restart/`, `out2d/`,
`*.modes.out` and the NetCDF diagnostics were never copied off Miyabi at all.
Each subdirectory carries a `.gitignore` that states its own rule and a README
that says what was dropped and why.

The protocol behind the determinism comparison — how the two codes' spectral
layouts are reconciled, what is being compared and what the numbers mean — is
documented once, in `gallope_dev/tests/determinism/README.md`, rather than
duplicated here.
