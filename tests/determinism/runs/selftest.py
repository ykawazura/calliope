#!/usr/bin/env python3
"""Self-test for compare.py.  Requires no Gallope/Calliope run and no GPU.

compare.py reconstructs the full spectral cube from two different half spaces
(Gallope halves z, Calliope halves y).  If that reconstruction were wrong the
comparison would be meaningless, so it is validated here against synthetic
data whose answer is known exactly.

The test builds one exactly-solenoidal, exactly-Hermitian field set, writes it
in both on-disk layouts, and checks three things:

  1. round trip -- compare.py must recover the source cube from either layout
     at roundoff, and must report the two layouts as identical;
  2. sensitivity -- a single mode perturbed by 1e-12 must rise above the
     roundoff floor;
  3. negative control -- a shifted axis and a swapped vector component must
     both be caught, the latter by the solenoidal self-check alone.

The same three checks are then run for RMHD, whose three fields (phi, omg,
psi) exercise the identical reconstruction machinery.  There the self-check
is omg = -kprp^2*phi rather than k.u = 0; the negative control is again a
shifted x axis, which that identity must catch because kprp^2 involves kx.

The grid is deliberately unequal (16 x 12 x 8) with an unequal box aspect, so
a transposed axis cannot hide behind a cubic grid.

Nyquist planes are excluded from the synthetic field.  The signed-wavenumber
convention shared by both codes assigns k = +n/2 to index n/2, whose mirror
index is itself, so k(mirror) = -k fails there and a spectrally projected
field would not be Hermitian on those planes.  That is a property of this
synthetic construction, not of the codes: in a real run the 2/3 rule zeroes
every mode with |k| >= (2/3) k_max and the Nyquist plane sits at |k| = k_max.

Usage:  python3 selftest.py            # exits non-zero on failure
"""

import os
import shutil
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare as C  # noqa: E402

NLX, NLY, NLZ = 16, 12, 8
LX, LY, LZ = 1.0, 2.0, 0.5
FIELDS = C.MODEL_FIELDS["MHD_INCOMP"]
RMHD_FIELDS = C.MODEL_FIELDS["RMHD"]
TOL = 1e-12          # relative, against the field scale


def keep_mask():
    """Modes retained in the synthetic field: no Nyquist plane, no mean."""
    kx, ky, kz = C.kgrids(NLX, NLY, NLZ, LX, LY, LZ)
    k2 = kx ** 2 + ky ** 2 + kz ** 2
    kxv, kyv, kzv = (C.code_wavenumbers(n) for n in (NLX, NLY, NLZ))
    return ((np.abs(kxv) != NLX // 2)[:, None, None]
            & (np.abs(kyv) != NLY // 2)[None, :, None]
            & (np.abs(kzv) != NLZ // 2)[None, None, :]) & (k2 > 0.0)


def hermitian_scalar(rng, keep):
    """A Hermitian spectral scalar: the transform of a real field, truncated."""
    c = np.fft.fftn(rng.standard_normal((NLX, NLY, NLZ)))
    c[~keep] = 0.0
    return c


def solenoidal_triple(rng):
    kx, ky, kz = C.kgrids(NLX, NLY, NLZ, LX, LY, LZ)
    k2 = kx ** 2 + ky ** 2 + kz ** 2
    keep = keep_mask()

    f = [np.fft.fftn(rng.standard_normal((NLX, NLY, NLZ))) for _ in range(3)]
    for c in f:
        c[~keep] = 0.0
    k2safe = np.where(k2 == 0.0, 1.0, k2)
    dot = kx * f[0] + ky * f[1] + kz * f[2]
    out = [f[0] - kx * dot / k2safe,
           f[1] - ky * dot / k2safe,
           f[2] - kz * dot / k2safe]
    for c in out:
        c[~keep] = 0.0
    return out


def write_fortran(arr, path):
    """MPI_ORDER_FORTRAN complex(8) dump: first index fastest.  read_raw()
    reverses this exactly."""
    np.asarray(arr, dtype=np.complex128).ravel(order="F").tofile(path)


def write_layout(cubes, code, restart_dir):
    os.makedirs(restart_dir, exist_ok=True)
    for f, c in cubes.items():
        if code == "gallope":                   # (nkz, nly, nlx), real FFT on z
            arr = c[:, :, : NLZ // 2 + 1].transpose(2, 1, 0)
        else:                                   # (nlx, nlz, nky), real FFT on y
            arr = c[:, : NLY // 2 + 1, :].transpose(0, 2, 1)
        write_fortran(arr, os.path.join(restart_dir, f + ".dat"))
    with open(os.path.join(restart_dir, "time.dat"), "w") as fh:
        fh.write("   tt                            tsc\n")
        fh.write("  1.0000000000000000E+00  0.0000000000000000E+00\n")


def rel_l2(a, b):
    nb = float(np.linalg.norm(b.ravel()))
    return float(np.linalg.norm((a - b).ravel())) / nb if nb else float("nan")


def main():
    rng = np.random.default_rng(20260804)
    n = (NLX, NLY, NLZ)
    kx, ky, kz = C.kgrids(*n, LX, LY, LZ)

    cubes = {}
    for prefix in ("u", "b"):
        fx, fy, fz = solenoidal_triple(rng)
        cubes[prefix + "x"], cubes[prefix + "y"], cubes[prefix + "z"] = fx, fy, fz
    scale = max(float(np.abs(cubes[f]).max()) for f in FIELDS)

    failures = []

    def check(label, value, limit, expect_above=False):
        ok = (value > limit) if expect_above else (value <= limit)
        print(f"  [{'ok ' if ok else 'FAIL'}] {label:44s} {value:.3e}")
        if not ok:
            failures.append(label)

    print(f"grid {NLX}x{NLY}x{NLZ}, box 2pi*({LX},{LY},{LZ}), field scale {scale:.3e}")

    print("\n-- source field --")
    for prefix in ("u", "b"):
        div = np.abs(kx * cubes[prefix + "x"] + ky * cubes[prefix + "y"]
                     + kz * cubes[prefix + "z"]).max()
        check(f"max|k.{prefix}| of the synthetic field", div, TOL * scale)

    root = tempfile.mkdtemp(prefix="gallope-determinism-selftest-")
    try:
        for tag, code in (("gallope_run", "gallope"), ("calliope_run", "calliope")):
            write_layout(cubes, code, os.path.join(root, tag, "restart"))

        loaded = {
            "gallope": C.load_gallope(
                os.path.join(root, "gallope_run", "restart"), *n, FIELDS),
            "calliope": C.load_calliope(
                os.path.join(root, "calliope_run", "restart"), *n, FIELDS),
        }

        # 1. round trip, split so a failure localises to one stage or the other
        print("\n-- round trip --")
        for code, axis in (("gallope", 2), ("calliope", 1)):
            worst_io = worst_exp = 0.0
            for f in FIELDS:
                sl = [slice(None)] * 3
                sl[axis] = slice(0, n[axis] // 2 + 1)
                half = cubes[f][tuple(sl)]
                worst_io = max(worst_io, float(np.abs(loaded[code][f] - half).max()))
                worst_exp = max(worst_exp, float(
                    np.abs(C.expand_hermitian(half, axis, *n) - cubes[f]).max()))
            check(f"{code}: max|loaded - source half|", worst_io, TOL * scale)
            check(f"{code}: max|expanded - source|", worst_exp, TOL * scale)

        A = C.full_cube(loaded["gallope"], "gallope", *n)
        B = C.full_cube(loaded["calliope"], "calliope", *n)

        print("\n-- the two layouts must agree --")
        check("solenoidal error, gallope cube",
              C.div_measure(A, "u", kx, ky, kz), TOL)
        check("solenoidal error, calliope cube",
              C.div_measure(B, "u", kx, ky, kz), TOL)
        combined = np.sqrt(sum(np.linalg.norm((A[f] - B[f]).ravel()) ** 2
                               for f in FIELDS)
                           / sum(np.linalg.norm(B[f].ravel()) ** 2 for f in FIELDS))
        check("relative field difference, all six fields", float(combined), TOL)
        wa = C.energy(A, "MHD_INCOMP", kx, ky)
        wb = C.energy(B, "MHD_INCOMP", kx, ky)
        check("relative energy difference", abs(wa - wb) / wb, TOL)

        # 2. sensitivity: the floor must sit far below a physically small error
        print("\n-- sensitivity (must rise above the roundoff floor) --")
        floor = rel_l2(A["ux"], B["ux"])
        pert = {f: v.copy() for f, v in loaded["calliope"].items()}
        pert["ux"][3, 2, 1] *= 1.0 + 1e-12
        signal = rel_l2(A["ux"], C.full_cube(pert, "calliope", *n)["ux"])
        check("one mode perturbed by 1e-12, signal / floor",
              signal / floor if floor else float("inf"), 10.0, expect_above=True)

        # 3. negative control: corruptions must not pass
        print("\n-- negative control (must be caught) --")
        shifted = C.full_cube({f: np.roll(v, 1, axis=0)
                               for f, v in loaded["calliope"].items()}, "calliope", *n)
        check("x shifted by one index, L2(ux)",
              rel_l2(A["ux"], shifted["ux"]), 1e-2, expect_above=True)
        check("x shifted by one index, solenoidal error",
              C.div_measure(shifted, "u", kx, ky, kz), 1e-10, expect_above=True)
        swapped = C.full_cube({**loaded["calliope"],
                               "uy": loaded["calliope"]["uz"],
                               "uz": loaded["calliope"]["uy"]}, "calliope", *n)
        check("uy/uz swapped, solenoidal error",
              C.div_measure(swapped, "u", kx, ky, kz), 1e-10, expect_above=True)

        # --- RMHD: same machinery, three fields, a different self-check -----
        print("\n-- RMHD (phi, omg, psi) --")
        keep = keep_mask()
        phi = hermitian_scalar(rng, keep)
        psi = hermitian_scalar(rng, keep)
        rcubes = {"phi": phi, "omg": -(kx ** 2 + ky ** 2) * phi, "psi": psi}
        rscale = max(float(np.abs(v).max()) for v in rcubes.values())

        for tag, code in (("rmhd_gallope", "gallope"), ("rmhd_calliope", "calliope")):
            write_layout(rcubes, code, os.path.join(root, tag, "restart"))
        rloaded = {
            "gallope": C.load_gallope(
                os.path.join(root, "rmhd_gallope", "restart"), *n, RMHD_FIELDS),
            "calliope": C.load_calliope(
                os.path.join(root, "rmhd_calliope", "restart"), *n, RMHD_FIELDS),
        }
        RA = C.full_cube(rloaded["gallope"], "gallope", *n)
        RB = C.full_cube(rloaded["calliope"], "calliope", *n)

        worst = max(float(np.abs(RA[f] - rcubes[f]).max()) for f in RMHD_FIELDS)
        check("round trip, max|gallope - source|", worst, TOL * rscale)
        worst = max(float(np.abs(RB[f] - rcubes[f]).max()) for f in RMHD_FIELDS)
        check("round trip, max|calliope - source|", worst, TOL * rscale)
        check("omg + kprp^2 phi, gallope cube",
              C.laplacian_identity(RA, kx, ky), TOL)
        check("omg + kprp^2 phi, calliope cube",
              C.laplacian_identity(RB, kx, ky), TOL)
        rcombined = np.sqrt(sum(np.linalg.norm((RA[f] - RB[f]).ravel()) ** 2
                                for f in RMHD_FIELDS)
                            / sum(np.linalg.norm(RB[f].ravel()) ** 2
                                  for f in RMHD_FIELDS))
        check("relative field difference, all three fields", float(rcombined), TOL)
        rwa = C.energy(RA, "RMHD", kx, ky)
        rwb = C.energy(RB, "RMHD", kx, ky)
        check("relative energy difference", abs(rwa - rwb) / rwb, TOL)

        rfloor = rel_l2(RA["phi"], RB["phi"])
        rpert = {f: v.copy() for f, v in rloaded["calliope"].items()}
        rpert["phi"][3, 2, 1] *= 1.0 + 1e-12
        rsignal = rel_l2(RA["phi"], C.full_cube(rpert, "calliope", *n)["phi"])
        check("one mode perturbed by 1e-12, signal / floor",
              rsignal / rfloor if rfloor else float("inf"), 10.0, expect_above=True)

        rshifted = C.full_cube({f: np.roll(v, 1, axis=0)
                                for f, v in rloaded["calliope"].items()},
                               "calliope", *n)
        check("x shifted by one index, L2(phi)",
              rel_l2(RA["phi"], rshifted["phi"]), 1e-2, expect_above=True)
        check("x shifted by one index, omg + kprp^2 phi",
              C.laplacian_identity(rshifted, kx, ky), 1e-10, expect_above=True)
    finally:
        shutil.rmtree(root, ignore_errors=True)

    print()
    if failures:
        print(f"FAILED ({len(failures)}): " + "; ".join(failures))
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
