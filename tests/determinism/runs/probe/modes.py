#!/usr/bin/env python3
"""Which modes carry the Gallope/Calliope difference in the linear shear probe?"""
import sys, numpy as np
sys.path.insert(0, ".")
import compare as C

N = (128, 128, 128)
L = (1.0, 1.0, 1.0)
A = C.full_cube(C.load_gallope("probe_lin_gallope/restart", *N), "gallope", *N)
B = C.full_cube(C.load_calliope("probe_lin_calliope/restart", *N), "calliope", *N)

kxv, kyv, kzv = (C.code_wavenumbers(n) for n in N)
for f in ("ux", "uy", "uz", "bx", "by", "bz"):
    d = np.abs(A[f] - B[f])
    tot = np.linalg.norm(A[f] - B[f])
    print(f"--- {f}:  ||df|| = {tot:.6e},  ||f|| = {np.linalg.norm(B[f]):.6e}")
    if tot == 0.0:
        continue
    idx = np.dstack(np.unravel_index(np.argsort(d.ravel())[::-1][:6], d.shape))[0]
    for i, j, k in idx:
        print(f"      k=({kxv[i]:+4d},{kyv[j]:+4d},{kzv[k]:+4d})  "
              f"|df|={d[i,j,k]:.4e}  |f_A|={abs(A[f][i,j,k]):.4e}  "
              f"|f_B|={abs(B[f][i,j,k]):.4e}  "
              f"arg(A/B)={np.angle(A[f][i,j,k]/B[f][i,j,k]):+.4e}"
              if abs(B[f][i, j, k]) > 0 else
              f"      k=({kxv[i]:+4d},{kyv[j]:+4d},{kzv[k]:+4d})  |df|={d[i,j,k]:.4e}  B=0")

# how many modes are non-negligible at all?
nz = np.abs(B["ux"]) > 1e-14 * np.abs(B["ux"]).max()
print(f"\nmodes with |ux| above 1e-14 of the peak: {nz.sum()}")
