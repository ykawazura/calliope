"""Where the runs of the tests live: tests/MHD_COMP_ISOTH/runs/. RUNS is relative to the current directory, so it is
tests/MHD_COMP_ISOTH/runs from the worktree root, where the gen scripts print CASES and the jobs are submitted.

Until 2026-11-13 the runs were in runs/ of the worktree root, which is now left to production runs; runs/... in the
REPORTs, the docstrings, the case lists and the json records means this directory.
"""
import os

RUNS = os.path.relpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'runs'))
