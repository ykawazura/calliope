#!/usr/bin/env python3
"""Normalized energy-balance residual of a determinism run, read through
ncdump because netCDF4 is not installed on this machine.

Residual = dW/dt + D - P
         = (u2dot + b2dot) + (u2dissip + b2dissip)
           - (p_ext_ene + p_re + p_ma)
exactly as diagnostics/MHD_INCOMP/plot_energy.py defines it.
"""
import subprocess, sys, re
import numpy as np

VARS = ("tt", "u2_sum", "b2_sum", "u2dot_sum", "b2dot_sum",
        "u2dissip_sum", "b2dissip_sum", "p_ext_ene_sum", "p_re_sum", "p_ma_sum")


def read(path):
    out = subprocess.run(["ncdump", "-v", ",".join(VARS), path],
                         capture_output=True, text=True).stdout
    body = out.split("data:", 1)[1]
    d = {}
    for v in VARS:
        m = re.search(r"^ %s = (.*?);" % re.escape(v), body, re.S | re.M)
        d[v] = np.array([float(x) for x in m.group(1).replace("\n", "").split(",")])
    return d


for path in sys.argv[1:]:
    d = read(path)
    W = 0.5 * (d["u2_sum"] + d["b2_sum"])
    dWdt = d["u2dot_sum"] + d["b2dot_sum"]
    D = d["u2dissip_sum"] + d["b2dissip_sum"]
    P = d["p_ext_ene_sum"] + d["p_re_sum"] + d["p_ma_sum"]
    r = dWdt + D - P
    # normalize by the dissipation, as in the response letter's part (d)
    scale = np.abs(D).mean()
    print(f"{path}")
    print(f"   samples {len(r)},  <D> = {scale:.6e},  <W> = {W.mean():.6e}")
    print(f"   residual rms = {np.sqrt((r ** 2).mean()):.3e}  "
          f"({np.sqrt((r ** 2).mean()) / scale * 100:.3f} % of <D>),  "
          f"max = {np.abs(r).max():.3e} ({np.abs(r).max() / scale * 100:.3f} %)")
