#!/usr/bin/env python3
"""Initial conditions of MHD_COMP_ISOTH in the restart format, and readers of its output.

The fields are given in real space on the grid of the code and written as the
Fourier coefficients that fields.F90 (restart) reads:

    restart_dir/{lnrho,mx,my,mz,bx,by,bz}.dat  complex*16, global array (nkx, nkz, nky)
                                               (Fortran order), f_k = FFT(f)/(nlx nly nlz)
    restart_dir/time.dat                       1 header line, then tt tsc in es30.21

Index convention of this module: real arrays are numpy arrays of shape
(nlx, nlz, nly) indexed [i, k, j] (the C view of the Fortran array (nly, nlz, nlx)),
and Fourier arrays are of shape (nkx, nkz, nky) indexed [i, k, j], as in the code.
y is the r2c direction (nky = nly/2 + 1).

    g = Grid(nlx=64, nly=2, nlz=64)            # 2D in the x-z plane
    lnrho_k = g.fft(0.1*np.sin(g.X + 2*g.Z))
    ...
    write_restart('runs/x/ic', g, lnrho_k, (mx_k, my_k, mz_k), (bx_k, by_k, bz_k))
"""
import os
import numpy as np


class Grid:
    """Grid and wavenumbers of grid.F90 (pruned = F), the 2/3 mask of dealias.F90 and the
    derivative wavenumbers k_d of shock_viscosity.F90 (Nyquist set to 0).
    lx, ly, lz are the input values (the box length is 2 pi lx)."""

    def __init__(self, nlx, nly, nlz, lx=1.0, ly=1.0, lz=1.0):
        self.nl = (nlx, nly, nlz)
        self.nlx, self.nly, self.nlz = nlx, nly, nlz
        self.Lx, self.Ly, self.Lz = 2*np.pi*lx, 2*np.pi*ly, 2*np.pi*lz
        self.nkx, self.nky, self.nkz = nlx, nly//2 + 1, nlz
        self.ntot = nlx*nly*nlz

        # grid.F90: xx(i) = lx (i-1)/nlx (after lx = 2 pi lx)
        self.x = self.Lx*np.arange(nlx)/nlx
        self.y = self.Ly*np.arange(nly)/nly
        self.z = self.Lz*np.arange(nlz)/nlz
        self.X, self.Z, self.Y = np.meshgrid(self.x, self.z, self.y, indexing='ij')

        # grid.F90: ikx(i) = i-1 for i <= nkx/2+1, else i-nkx-1 (the Nyquist mode is +nlx/2)
        def full(n):
            i = np.arange(1, n + 1)
            return np.where(i <= n//2 + 1, i - 1, i - n - 1)
        self.ikx, self.ikz = full(self.nkx), full(self.nkz)
        self.iky = np.arange(self.nky)
        kx = 2*np.pi*self.ikx/self.Lx
        ky = 2*np.pi*self.iky/self.Ly
        kz = 2*np.pi*self.ikz/self.Lz
        self.kx, self.ky, self.kz = (a.reshape(s) for a, s in
                                     ((kx, (-1, 1, 1)), (ky, (1, 1, -1)), (kz, (1, -1, 1))))

        # dealias.F90: zero if |k_d| >= (2/3) 2 pi (nl/2)/L in any direction
        kmx = 2*np.pi*(nlx//2)/self.Lx
        kmy = 2*np.pi*(nly//2)/self.Ly
        kmz = 2*np.pi*(nlz//2)/self.Lz
        self.mask = ((np.abs(self.kx) < kmx*2/3) & (np.abs(self.ky) < kmy*2/3)
                     & (np.abs(self.kz) < kmz*2/3)).astype(float)

        # shock_viscosity.F90: k_d = k with the Nyquist mode set to 0
        self.kxd = np.where(np.abs(self.ikx) == nlx//2, 0.0, kx).reshape(-1, 1, 1)
        self.kyd = np.where(np.abs(self.iky) == nly//2, 0.0, ky).reshape(1, 1, -1)
        self.kzd = np.where(np.abs(self.ikz) == nlz//2, 0.0, kz).reshape(1, -1, 1)

    def zeros(self):
        return np.zeros((self.nlx, self.nlz, self.nly))

    def fft(self, f):
        """Real (nlx, nlz, nly) -> Fourier (nkx, nkz, nky), normalized as the code."""
        return np.fft.rfftn(np.broadcast_to(f, (self.nlx, self.nlz, self.nly)), axes=(0, 1, 2))/self.ntot

    def ifft(self, fk):
        return np.fft.irfftn(fk*self.ntot, s=(self.nlx, self.nlz, self.nly), axes=(0, 1, 2))

    def project_divfree(self, bk):
        """Remove the k-parallel part with the k of the code (div_free in advance.F90)."""
        bx, by, bz = bk
        k2 = self.kx**2 + self.ky**2 + self.kz**2
        k2inv = np.where(k2 > 0, 1/np.where(k2 > 0, k2, 1), 0.0)
        phi = (self.kx*bx + self.ky*by + self.kz*bz)*k2inv
        return bx - self.kx*phi, by - self.ky*phi, bz - self.kz*phi

    def div_d(self, vk):
        """div_N v = F^-1 i k_d . v_k in real space (no mask applied)."""
        vx, vy, vz = vk
        return self.ifft(1j*(self.kxd*vx + self.kyd*vy + self.kzd*vz))


def state_from_real(g, lnrho, u, b, divfree=True):
    """Fourier state (lnrho_k, m_k, b_k) of real lnrho, u = (ux, uy, uz) and b = (bx, by, bz),
    with the 2/3 mask applied to m and B (as apply_filter with filter_mom = filter_mag = T)
    and B projected on div B = 0. lnrho is not filtered (filter_lnrho = F)."""
    rho = np.exp(lnrho)
    lk = g.fft(lnrho)
    mk = tuple(g.fft(rho*ui)*g.mask for ui in u)
    bk = tuple(g.fft(bi)*g.mask for bi in b)
    if divfree:
        bk = tuple(bi*g.mask for bi in g.project_divfree(bk))
    return lk, mk, bk


def real_from_state(g, lk, mk, bk):
    """rho, u, B in real space as the code computes them (rho = exp(F^-1 lnrho_k), u = m/rho)."""
    rho = np.exp(g.ifft(lk))
    u = tuple(g.ifft(mi)/rho for mi in mk)
    b = tuple(g.ifft(bi) for bi in bk)
    return rho, u, b


def _write_k(path, fk):
    # (nkx, nkz, nky) [i, k, j] -> Fortran order: i fastest
    np.ascontiguousarray(fk.transpose(2, 1, 0)).astype('<c16').tofile(path)


def _read_k(path, g):
    a = np.fromfile(path, dtype='<c16')
    a = a.reshape(-1, g.nky, g.nkz, g.nkx).transpose(0, 3, 2, 1)
    return a


NAMES = ('lnrho', 'mx', 'my', 'mz', 'bx', 'by', 'bz')


def write_restart(path, g, lk, mk, bk, tt=0.0, tsc=0.0):
    os.makedirs(path, exist_ok=True)
    for name, fk in zip(NAMES, (lk,) + tuple(mk) + tuple(bk)):
        assert fk.shape == (g.nkx, g.nkz, g.nky), (name, fk.shape)
        _write_k(os.path.join(path, name + '.dat'), fk)
    with open(os.path.join(path, 'time.dat'), 'w') as f:
        f.write('   tt                            tst\n')
        f.write('%30.21E%30.21E\n' % (tt, tsc))


def read_restart(path, g):
    """(lnrho_k, m_k, b_k, tt) of a restart directory."""
    fk = {n: _read_k(os.path.join(path, n + '.dat'), g)[0] for n in NAMES}
    with open(os.path.join(path, 'time.dat')) as f:
        f.readline()
        tt = float(f.readline().split()[0])
    return fk['lnrho'], (fk['mx'], fk['my'], fk['mz']), (fk['bx'], fk['by'], fk['bz']), tt


def _read_out3d_k(path, g, nrec):
    # mpiio_write_var (src/mpiio_write_var.F90) takes bytes = MPI_TYPE_SIZE(MPI_DOUBLE_COMPLEX) = 16 and then
    # adds the complex part once more, so record r starts at 2 r records and the gaps between are never
    # written (zero); a file without the gaps is accepted as well
    a = np.fromfile(path, dtype='<c16').reshape(-1, g.nky*g.nkz*g.nkx)
    if len(a) == 2*nrec - 1 and nrec > 1:
        assert not np.any(a[1::2]), path + ': nonzero data between the records'
        a = a[::2]
    assert len(a) == nrec, '%s: %d records for %d times' % (path, len(a), nrec)
    return a.reshape(-1, g.nky, g.nkz, g.nkx).transpose(0, 3, 2, 1)


def read_out3d(path, g):
    """Snapshots of out3d/ (Fourier, as restart; rho.dat holds rho, not ln rho):
    dict name -> array (nrec, nkx, nkz, nky) and 't' -> times."""
    t = np.atleast_2d(np.loadtxt(os.path.join(path, 'time.dat')))[:, 0]
    out = {'t': t}
    for n in ('rho', 'mx', 'my', 'mz', 'bx', 'by', 'bz'):
        out[n] = _read_out3d_k(os.path.join(path, n + '.dat'), g, len(t))
    return out


PLANE_Q = ('rho', 'ux', 'uy', 'uz', 'bx', 'by', 'bz', 'zeta', 'D', 'eta')


def read_planes(path, g, plane='y0'):
    """Records of out2d/plane_<plane>.dat written by write_planes (advance.F90):
    dict q -> array (nrec, nlx, nlz) for y0 or (nrec, nlx, nly) for z0, and 't', 'step'."""
    tab = np.atleast_2d(np.loadtxt(os.path.join(path, 'plane_time.dat')))
    n1, n2 = (g.nlx, g.nlz) if plane == 'y0' else (g.nlx, g.nly)
    a = np.fromfile(os.path.join(path, 'plane_%s.dat' % plane), dtype='<f8')
    a = a.reshape(-1, len(PLANE_Q), n2, n1).transpose(0, 1, 3, 2)  # x fastest in the file
    out = {q: a[:, iq] for iq, q in enumerate(PLANE_Q)}
    out['t'], out['step'] = tab[:len(a), 1], tab[:len(a), 2].astype(int)
    return out
