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
aaa2_sum        = np.copy(ncfile.variables['aaa2_sum'       ][:]); aaa2_sum        = np.delete(aaa2_sum       , ignored_points, axis = 0)
upe2dot_sum     = np.copy(ncfile.variables['upe2dot_sum'    ][:]); upe2dot_sum     = np.delete(upe2dot_sum    , ignored_points, axis = 0)
bpe2dot_sum     = np.copy(ncfile.variables['bpe2dot_sum'    ][:]); bpe2dot_sum     = np.delete(bpe2dot_sum    , ignored_points, axis = 0)
aaa2dot_sum     = np.copy(ncfile.variables['aaa2dot_sum'    ][:]); aaa2dot_sum     = np.delete(aaa2dot_sum    , ignored_points, axis = 0)
upe2dissip_sum  = np.copy(ncfile.variables['upe2dissip_sum' ][:]); upe2dissip_sum  = np.delete(upe2dissip_sum , ignored_points, axis = 0)
bpe2dissip_sum  = np.copy(ncfile.variables['bpe2dissip_sum' ][:]); bpe2dissip_sum  = np.delete(bpe2dissip_sum , ignored_points, axis = 0)
aaa2dissip_sum  = np.copy(ncfile.variables['aaa2dissip_sum' ][:]); aaa2dissip_sum  = np.delete(aaa2dissip_sum , ignored_points, axis = 0)
p_sum           = np.copy(ncfile.variables['p_sum'          ][:]); p_sum           = np.delete(p_sum          , ignored_points, axis = 0)
zpep2_sum       = np.copy(ncfile.variables['zpep2_sum'      ][:]); zpep2_sum       = np.delete(zpep2_sum      , ignored_points, axis = 0)
zpem2_sum       = np.copy(ncfile.variables['zpem2_sum'      ][:]); zpem2_sum       = np.delete(zpem2_sum      , ignored_points, axis = 0)

# Load binned spectra
upe2_bin        = np.copy(ncfile.variables['upe2_bin'       ][:]); upe2_bin        = np.delete(upe2_bin       , ignored_points, axis = 0)
bpe2_bin        = np.copy(ncfile.variables['bpe2_bin'       ][:]); bpe2_bin        = np.delete(bpe2_bin       , ignored_points, axis = 0)
aaa2_bin        = np.copy(ncfile.variables['aaa2_bin'       ][:]); aaa2_bin        = np.delete(aaa2_bin       , ignored_points, axis = 0)
ux2_bin         = np.copy(ncfile.variables['ux2_bin'        ][:]); ux2_bin         = np.delete(ux2_bin        , ignored_points, axis = 0)
uy2_bin         = np.copy(ncfile.variables['uy2_bin'        ][:]); uy2_bin         = np.delete(uy2_bin        , ignored_points, axis = 0)
bx2_bin         = np.copy(ncfile.variables['bx2_bin'        ][:]); bx2_bin         = np.delete(bx2_bin        , ignored_points, axis = 0)
by2_bin         = np.copy(ncfile.variables['by2_bin'        ][:]); by2_bin         = np.delete(by2_bin        , ignored_points, axis = 0)
p_bin           = np.copy(ncfile.variables['p_bin'          ][:]); p_bin           = np.delete(p_bin          , ignored_points, axis = 0)

zpep2_bin       = np.copy(ncfile.variables['zpep2_bin'      ][:]); zpep2_bin       = np.delete(zpep2_bin      , ignored_points, axis = 0)
zpem2_bin       = np.copy(ncfile.variables['zpem2_bin'      ][:]); zpem2_bin       = np.delete(zpem2_bin      , ignored_points, axis = 0)

ncfile.close()

tlab  = r'$Nt$'
zlab  = r'$zN/v_\mathrm{A}$'
kzlab = r'$k_z v_\mathrm{A}/N$'
