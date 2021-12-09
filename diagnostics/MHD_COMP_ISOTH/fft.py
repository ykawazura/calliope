# -*- coding: utf-8 -*-
import numpy as np
from numba import jit
from numpy import fft
from load  import *

##### Function to define inverse fft properly, with hermitian symmetries
##### the input will be f(kx,ky,kz), kx, ky and kz; the output will be the FT and
##### the axes x and y

@jit
def fft_fwd3d(fld): #we assume ky and kz have positive and negative wavenumbers, 
                    #kx only positive => add conjugates

  #create the fftfreq-like array
  kx_full = np.hstack([kx, -kx[:0:-1]])
  nkx_full = kx_full.size
    
  ###enforce hermitian symmetry
  fld_full = np.zeros(nkx_full*nky*nkz, dtype=complex).reshape(nkx_full, nky, nkz)
  # fill kx >= 0 entries
  fld_full[:nkx, :]=fld[:nkx, :]

  # generate kx < 0 entries using Hermiticity
  # i.e. fld(-k) = conjg(fld(k))
  for i in range(-nkx + 1, 0):
    # treat special case of ky=0 & ky=0
    fld_full[i + nkx_full, 0, 0] = np.conj(fld[-i, 0, 0])
    # treat special case of ky=0
    for k in range(1, nkz):
      fld_full[i + nkx_full, 0, k] = np.conj(fld[-i, 0, nkz - k])
    # treat special case of kz = 0
    for j in range(1, nky):
      fld_full[i + nkx_full, j, 0] = np.conj(fld[-i, nky - j, 0])
    # next treat general case of kx /= 0 and ky /= 0
    for j in range(1, nky):
      for k in range(1, nkz):
        fld_full[i + nkx_full, j, k] = np.conj(fld[-i, nky - j, nkz - k])
    
  #calculate ifftn
  fft3d = np.real(np.fft.ifftn(fld_full))

  #normalize
  fft3d = fft3d*(nkx_full*nky*nkz)

  return fft3d

@jit
def fft_fwd3d_t(fld): #we assume ky and kz have positive and negative wavenumbers, 
                      #kx only positive => add conjugates

  ntfield = fld.shape[-1]
  #create the fftfreq-like array
  kx_full = np.hstack([kx, -kx[:0:-1]])
  nkx_full = kx_full.size
    
  ###enforce hermitian symmetry
  fld_full = np.zeros(nkx_full*nky*nkz*ntfield, dtype=complex).reshape(nkx_full, nky, nkz, ntfield)
  # fill kx >= 0 entries
  fld_full[:nkx, :]=fld[:nkx, :]

  # generate kx < 0 entries using Hermiticity
  # i.e. fld(-k) = conjg(fld(k))
  for i in range(-nkx + 1, 0):
    # treat special case of ky=0 & ky=0
    fld_full[i + nkx_full, 0, 0, :] = np.conj(fld[-i, 0, 0, :])
    # treat special case of ky=0
    for k in range(1, nkz):
      fld_full[i + nkx_full, 0, k, :] = np.conj(fld[-i, 0, nkz - k, :])
    # treat special case of kz = 0
    for j in range(1, nky):
      fld_full[i + nkx_full, j, 0, :] = np.conj(fld[-i, nky - j, 0, :])
    # next treat general case of kx /= 0 and ky /= 0
    for j in range(1, nky):
      for k in range(1, nkz):
        fld_full[i + nkx_full, j, k, :] = np.conj(fld[-i, nky - j, nkz - k, :])
    
  #calculate ifftn
  fft3d = np.real(np.fft.ifftn(fld_full, axes=(0, 1, 2)))

  #normalize
  fft3d = fft3d*(nkx_full*nky*nkz)

  return fft3d

# @jit
def sum_negative_kz2d(fld):
  fld_new = np.zeros([nt, int(nkz/2), nkpolar])

  fld_new[:, 0,  :] = fld[:, 0, :]
  for i in np.arange(1, int(nkz/2)):
    fld_new[:, i, :] = 0.5*(fld[:, i, :] + fld[:, nkz - i, :])
  return fld_new
