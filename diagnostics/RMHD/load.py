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
ismovie = True

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
upe2dot_sum     = np.copy(ncfile.variables['upe2dot_sum'    ][:]); upe2dot_sum     = np.delete(upe2dot_sum    , ignored_points, axis = 0)
bpe2dot_sum     = np.copy(ncfile.variables['bpe2dot_sum'    ][:]); bpe2dot_sum     = np.delete(bpe2dot_sum    , ignored_points, axis = 0)
upe2dissip_sum  = np.copy(ncfile.variables['upe2dissip_sum' ][:]); upe2dissip_sum  = np.delete(upe2dissip_sum , ignored_points, axis = 0)
bpe2dissip_sum  = np.copy(ncfile.variables['bpe2dissip_sum' ][:]); bpe2dissip_sum  = np.delete(bpe2dissip_sum , ignored_points, axis = 0)
p_omg_sum       = np.copy(ncfile.variables['p_omg_sum'      ][:]); p_omg_sum       = np.delete(p_omg_sum      , ignored_points, axis = 0)
p_psi_sum       = np.copy(ncfile.variables['p_psi_sum'      ][:]); p_psi_sum       = np.delete(p_psi_sum      , ignored_points, axis = 0)

# Load binned spectra
upe2_bin        = np.copy(ncfile.variables['upe2_bin'       ][:]); upe2_bin        = np.delete(upe2_bin       , ignored_points, axis = 0)
bpe2_bin        = np.copy(ncfile.variables['bpe2_bin'       ][:]); bpe2_bin        = np.delete(bpe2_bin       , ignored_points, axis = 0)

ncfile.close()

# Load 2D cut
ncfile = netcdf.netcdf_file(input_dir+runname+'.out.fields_section.nc'+restart_num, 'r')   

tt_fld  = np.copy(ncfile.variables['tt' ][:]); tt_fld  = np.delete(tt_fld , ignored_points_fld, axis = 0)
xx_fld  = np.copy(ncfile.variables['xx' ][:])
yy_fld  = np.copy(ncfile.variables['yy' ][:])
zz_fld  = np.copy(ncfile.variables['zz' ][:])

phi_r_z0        = np.copy(ncfile.variables['phi_r_z0'       ][:]); phi_r_z0        = np.delete(phi_r_z0       , ignored_points_fld, axis = 0)
phi_r_x0        = np.copy(ncfile.variables['phi_r_x0'       ][:]); phi_r_x0        = np.delete(phi_r_x0       , ignored_points_fld, axis = 0)
phi_r_y0        = np.copy(ncfile.variables['phi_r_y0'       ][:]); phi_r_y0        = np.delete(phi_r_y0       , ignored_points_fld, axis = 0)

psi_r_z0        = np.copy(ncfile.variables['psi_r_z0'       ][:]); psi_r_z0        = np.delete(psi_r_z0       , ignored_points_fld, axis = 0)
psi_r_x0        = np.copy(ncfile.variables['psi_r_x0'       ][:]); psi_r_x0        = np.delete(psi_r_x0       , ignored_points_fld, axis = 0)
psi_r_y0        = np.copy(ncfile.variables['psi_r_y0'       ][:]); psi_r_y0        = np.delete(psi_r_y0       , ignored_points_fld, axis = 0)

omg_r_z0        = np.copy(ncfile.variables['omg_r_z0'       ][:]); omg_r_z0        = np.delete(omg_r_z0       , ignored_points_fld, axis = 0)
omg_r_x0        = np.copy(ncfile.variables['omg_r_x0'       ][:]); omg_r_x0        = np.delete(omg_r_x0       , ignored_points_fld, axis = 0)
omg_r_y0        = np.copy(ncfile.variables['omg_r_y0'       ][:]); omg_r_y0        = np.delete(omg_r_y0       , ignored_points_fld, axis = 0)

jpa_r_z0        = np.copy(ncfile.variables['jpa_r_z0'       ][:]); jpa_r_z0        = np.delete(jpa_r_z0       , ignored_points_fld, axis = 0)
jpa_r_x0        = np.copy(ncfile.variables['jpa_r_x0'       ][:]); jpa_r_x0        = np.delete(jpa_r_x0       , ignored_points_fld, axis = 0)
jpa_r_y0        = np.copy(ncfile.variables['jpa_r_y0'       ][:]); jpa_r_y0        = np.delete(jpa_r_y0       , ignored_points_fld, axis = 0)

ux_r_z0         = np.copy(ncfile.variables[ 'ux_r_z0'       ][:]);  ux_r_z0        = np.delete( ux_r_z0       , ignored_points_fld, axis = 0)
ux_r_x0         = np.copy(ncfile.variables[ 'ux_r_x0'       ][:]);  ux_r_x0        = np.delete( ux_r_x0       , ignored_points_fld, axis = 0)
ux_r_y0         = np.copy(ncfile.variables[ 'ux_r_y0'       ][:]);  ux_r_y0        = np.delete( ux_r_y0       , ignored_points_fld, axis = 0)

uy_r_z0         = np.copy(ncfile.variables[ 'uy_r_z0'       ][:]);  uy_r_z0        = np.delete( uy_r_z0       , ignored_points_fld, axis = 0)
uy_r_x0         = np.copy(ncfile.variables[ 'uy_r_x0'       ][:]);  uy_r_x0        = np.delete( uy_r_x0       , ignored_points_fld, axis = 0)
uy_r_y0         = np.copy(ncfile.variables[ 'uy_r_y0'       ][:]);  uy_r_y0        = np.delete( uy_r_y0       , ignored_points_fld, axis = 0)

bx_r_z0         = np.copy(ncfile.variables[ 'bx_r_z0'       ][:]);  bx_r_z0        = np.delete( bx_r_z0       , ignored_points_fld, axis = 0)
bx_r_x0         = np.copy(ncfile.variables[ 'bx_r_x0'       ][:]);  bx_r_x0        = np.delete( bx_r_x0       , ignored_points_fld, axis = 0)
bx_r_y0         = np.copy(ncfile.variables[ 'bx_r_y0'       ][:]);  bx_r_y0        = np.delete( bx_r_y0       , ignored_points_fld, axis = 0)

by_r_z0         = np.copy(ncfile.variables[ 'by_r_z0'       ][:]);  by_r_z0        = np.delete( by_r_z0       , ignored_points_fld, axis = 0)
by_r_x0         = np.copy(ncfile.variables[ 'by_r_x0'       ][:]);  by_r_x0        = np.delete( by_r_x0       , ignored_points_fld, axis = 0)
by_r_y0         = np.copy(ncfile.variables[ 'by_r_y0'       ][:]);  by_r_y0        = np.delete( by_r_y0       , ignored_points_fld, axis = 0)

ncfile.close()

tlab  = r'(v_\rmA/L_\|)t'
zlab  = r'z/L_\|'
kzlab = r'L_\| k_z'
