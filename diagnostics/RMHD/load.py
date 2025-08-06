# -*- coding: utf-8 -*-
import warnings
warnings.filterwarnings('ignore')

import numpy as np
from scipy.io import netcdf

####### time index for final cut #########
input_dir = '../../'
runname = 'calliope'
restart_num = ''
final_idx     = -1
final_fld_idx = -1

if final_idx != -1:
  print ('\n!!! CAUTION: final_idx = %d !!!\n' % final_idx)

####### movie or final cut #########
ismovie = False

####### ignore these points #########
ignored_points     = [0]
ignored_points_fld = [0]


####### load netcdf file #########
ncfile = netcdf.netcdf_file(input_dir+runname+'.out.nc'+restart_num, 'r')   

# Load parameters

# Load coordinate
tt  = np.copy(ncfile.variables['tt' ][:]); tt  = np.delete(tt , ignored_points, axis = 0)
xx  = np.copy(ncfile.variables['xx' ][:])
yy  = np.copy(ncfile.variables['yy' ][:])
zz  = np.copy(ncfile.variables['zz' ][:])
kx  = np.copy(ncfile.variables['kx' ][:])
ky  = np.copy(ncfile.variables['ky' ][:])
kz  = np.copy(ncfile.variables['kz' ][:])
kpbin  = np.copy(ncfile.variables['kpbin'][:])

nt  = tt.size
nlx = xx.size
nly = yy.size
nlz = zz.size
nkx = kx.size
nky = ky.size
nkz = kz.size
nkpolar = kpbin.size

if nkz <= 2:
  is2D = True
else:
  is2D = False

# Load total energies
upe2_sum        = np.copy(ncfile.variables['upe2_sum'       ][:]); upe2_sum        = np.delete(upe2_sum       , ignored_points, axis = 0)
bpe2_sum        = np.copy(ncfile.variables['bpe2_sum'       ][:]); bpe2_sum        = np.delete(bpe2_sum       , ignored_points, axis = 0)
zppe2_sum       = np.copy(ncfile.variables['zppe2_sum'      ][:]); zppe2_sum       = np.delete(zppe2_sum      , ignored_points, axis = 0)
zmpe2_sum       = np.copy(ncfile.variables['zmpe2_sum'      ][:]); zmpe2_sum       = np.delete(zmpe2_sum      , ignored_points, axis = 0)
upe2dot_sum     = np.copy(ncfile.variables['upe2dot_sum'    ][:]); upe2dot_sum     = np.delete(upe2dot_sum    , ignored_points, axis = 0)
bpe2dot_sum     = np.copy(ncfile.variables['bpe2dot_sum'    ][:]); bpe2dot_sum     = np.delete(bpe2dot_sum    , ignored_points, axis = 0)
upe2dissip_sum  = np.copy(ncfile.variables['upe2dissip_sum' ][:]); upe2dissip_sum  = np.delete(upe2dissip_sum , ignored_points, axis = 0)
bpe2dissip_sum  = np.copy(ncfile.variables['bpe2dissip_sum' ][:]); bpe2dissip_sum  = np.delete(bpe2dissip_sum , ignored_points, axis = 0)
p_phi_sum       = np.copy(ncfile.variables['p_phi_sum'      ][:]); p_phi_sum       = np.delete(p_phi_sum      , ignored_points, axis = 0)
p_psi_sum       = np.copy(ncfile.variables['p_psi_sum'      ][:]); p_psi_sum       = np.delete(p_psi_sum      , ignored_points, axis = 0)
p_xhl_sum       = np.copy(ncfile.variables['p_xhl_sum'      ][:]); p_xhl_sum       = np.delete(p_xhl_sum      , ignored_points, axis = 0)

# Load binned spectra
upe2_bin        = np.copy(ncfile.variables['upe2_bin'       ][:]); upe2_bin        = np.delete(upe2_bin       , ignored_points, axis = 0)
bpe2_bin        = np.copy(ncfile.variables['bpe2_bin'       ][:]); bpe2_bin        = np.delete(bpe2_bin       , ignored_points, axis = 0)
zppe2_bin       = np.copy(ncfile.variables['zppe2_bin'      ][:]); zppe2_bin       = np.delete(zppe2_bin      , ignored_points, axis = 0)
zmpe2_bin       = np.copy(ncfile.variables['zmpe2_bin'      ][:]); zmpe2_bin       = np.delete(zmpe2_bin      , ignored_points, axis = 0)
hp_bin          = np.copy(ncfile.variables['hp_bin'         ][:]); hp_bin          = np.delete(hp_bin         , ignored_points, axis = 0)
hm_bin          = np.copy(ncfile.variables['hm_bin'         ][:]); hm_bin          = np.delete(hm_bin         , ignored_points, axis = 0)

ncfile.close()

tlab  = r'$(v_\mathrm{A}/L_\|)t$'
zlab  = r'$z/L_\|$'
kzlab = r'$L_\| k_z$'

# Load rms
import os
filename = input_dir+'rms.dat'+restart_num
if os.path.isfile(filename):
  data = np.loadtxt(filename).T

  # Reduce multiple points at the same time because of multistep update (RK or Gear)
  _, indices = np.unique(data[0], return_index=True)
  data = np.asarray([data.T[i] for i in indices]).T
  tt_rms, ux_rms, uy_rms, bx_rms, by_rms = data

  # Pick up where tt_rms = tt
  indices = [np.abs(tt_rms - t).argmin() for t in tt]
  tt_rms, ux_rms, uy_rms, bx_rms, by_rms = np.asarray([data.T[i] for i in indices]).T
else:
  tt_rms, ux_rms, uy_rms, bx_rms, by_rms = [np.array([0])] * 5

