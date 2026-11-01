/* LD_PRELOAD shim for the FFTW legacy Fortran planner calls of p3dfft (dfftw_plan_many_*_).
   Logs the planner flags of each plan to stderr; with FFTW_SHIM_ESTIMATE=1 it replaces
   FFTW_MEASURE (0), FFTW_PATIENT and FFTW_EXHAUSTIVE by FFTW_ESTIMATE, whose plan does not
   depend on timings.  Test use only: with FFTW_MEASURE two runs of the same input differ at
   ~1e-13 (the planner times candidate algorithms), so a bit-identity test of the restart needs
   the same plans in both runs.  Built and preloaded by job.pbs-ps-fes-dev for cases with a file
   "fftw_estimate":  gcc -O2 -shared -fPIC -o fftw_shim.so fftw_shim.c -ldl */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>

#define FFTW_EXHAUSTIVE (1 << 3)
#define FFTW_PATIENT (1 << 5)
#define FFTW_ESTIMATE (1 << 6)

static int newflags(const char *name, const int *n, const int *howmany, const int *flags)
{
  int f = *flags;
  const char *e = getenv("FFTW_SHIM_ESTIMATE");
  if (e && e[0] == '1') f = (f & ~(FFTW_PATIENT | FFTW_EXHAUSTIVE)) | FFTW_ESTIMATE;
  fprintf(stderr, "fftw_shim: %-4s n = %d howmany = %d flags %d -> %d\n", name, n[0], *howmany, *flags, f);
  return f;
}

static void *next(const char *sym)
{
  void *p = dlsym(RTLD_NEXT, sym);
  if (!p) { fprintf(stderr, "fftw_shim: %s not found\n", sym); abort(); }
  return p;
}

typedef void (*dft_t)(void *, int *, const int *, int *, void *, const int *, int *, int *,
                      void *, const int *, int *, int *, int *, int *);
typedef void (*rc_t)(void *, int *, const int *, int *, void *, const int *, int *, int *,
                     void *, const int *, int *, int *, int *);

void dfftw_plan_many_dft_(void *p, int *rank, const int *n, int *howmany, void *in, const int *inembed,
                          int *istride, int *idist, void *out, const int *onembed, int *ostride, int *odist,
                          int *sign, int *flags)
{
  static dft_t f0 = 0;
  if (!f0) f0 = (dft_t)next("dfftw_plan_many_dft_");
  int f = newflags("dft", n, howmany, flags);
  f0(p, rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, sign, &f);
}

void dfftw_plan_many_r2r_(void *p, int *rank, const int *n, int *howmany, void *in, const int *inembed,
                          int *istride, int *idist, void *out, const int *onembed, int *ostride, int *odist,
                          int *kind, int *flags)
{
  static dft_t f0 = 0;
  if (!f0) f0 = (dft_t)next("dfftw_plan_many_r2r_");
  int f = newflags("r2r", n, howmany, flags);
  f0(p, rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, kind, &f);
}

void dfftw_plan_many_dft_r2c_(void *p, int *rank, const int *n, int *howmany, void *in, const int *inembed,
                              int *istride, int *idist, void *out, const int *onembed, int *ostride, int *odist,
                              int *flags)
{
  static rc_t f0 = 0;
  if (!f0) f0 = (rc_t)next("dfftw_plan_many_dft_r2c_");
  int f = newflags("r2c", n, howmany, flags);
  f0(p, rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, &f);
}

void dfftw_plan_many_dft_c2r_(void *p, int *rank, const int *n, int *howmany, void *in, const int *inembed,
                              int *istride, int *idist, void *out, const int *onembed, int *ostride, int *odist,
                              int *flags)
{
  static rc_t f0 = 0;
  if (!f0) f0 = (rc_t)next("dfftw_plan_many_dft_c2r_");
  int f = newflags("c2r", n, howmany, flags);
  f0(p, rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, &f);
}
