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
beta   = np.copy(ncfile.variables['beta'].data)
gamma  = np.copy(ncfile.variables['gamma'].data)
q      = np.copy(ncfile.variables['q'].data)
va2cs2_plus_1 = 2.0/(beta*gamma) + 1.0

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
upa2_sum        = np.copy(ncfile.variables['upa2_sum'       ][:]); upa2_sum        = np.delete(upa2_sum       , ignored_points, axis = 0)
bpa2_sum        = np.copy(ncfile.variables['bpa2_sum'       ][:]); bpa2_sum        = np.delete(bpa2_sum       , ignored_points, axis = 0)
upe2dot_sum     = np.copy(ncfile.variables['upe2dot_sum'    ][:]); upe2dot_sum     = np.delete(upe2dot_sum    , ignored_points, axis = 0)
bpe2dot_sum     = np.copy(ncfile.variables['bpe2dot_sum'    ][:]); bpe2dot_sum     = np.delete(bpe2dot_sum    , ignored_points, axis = 0)
upa2dot_sum     = np.copy(ncfile.variables['upa2dot_sum'    ][:]); upa2dot_sum     = np.delete(upa2dot_sum    , ignored_points, axis = 0)
bpa2dot_sum     = np.copy(ncfile.variables['bpa2dot_sum'    ][:]); bpa2dot_sum     = np.delete(bpa2dot_sum    , ignored_points, axis = 0)
upe2dissip_sum  = np.copy(ncfile.variables['upe2dissip_sum' ][:]); upe2dissip_sum  = np.delete(upe2dissip_sum , ignored_points, axis = 0)
bpe2dissip_sum  = np.copy(ncfile.variables['bpe2dissip_sum' ][:]); bpe2dissip_sum  = np.delete(bpe2dissip_sum , ignored_points, axis = 0)
upa2dissip_sum  = np.copy(ncfile.variables['upa2dissip_sum' ][:]); upa2dissip_sum  = np.delete(upa2dissip_sum , ignored_points, axis = 0)
bpa2dissip_sum  = np.copy(ncfile.variables['bpa2dissip_sum' ][:]); bpa2dissip_sum  = np.delete(bpa2dissip_sum , ignored_points, axis = 0)
p_aw_sum        = np.copy(ncfile.variables['p_aw_sum'       ][:]); p_aw_sum        = np.delete(p_aw_sum       , ignored_points, axis = 0)
p_compr_sum     = np.copy(ncfile.variables['p_compr_sum'    ][:]); p_compr_sum     = np.delete(p_compr_sum    , ignored_points, axis = 0)
p_compr_sum     = np.copy(ncfile.variables['p_compr_sum'    ][:]); p_compr_sum     = np.delete(p_compr_sum    , ignored_points, axis = 0)
zpep2_sum       = np.copy(ncfile.variables['zpep2_sum'      ][:]); zpep2_sum       = np.delete(zpep2_sum      , ignored_points, axis = 0)
zpem2_sum       = np.copy(ncfile.variables['zpem2_sum'      ][:]); zpem2_sum       = np.delete(zpem2_sum      , ignored_points, axis = 0)
zpap2_sum       = np.copy(ncfile.variables['zpap2_sum'      ][:]); zpap2_sum       = np.delete(zpap2_sum      , ignored_points, axis = 0)
zpam2_sum       = np.copy(ncfile.variables['zpam2_sum'      ][:]); zpam2_sum       = np.delete(zpam2_sum      , ignored_points, axis = 0)

# Load binned spectra
upe2_bin        = np.copy(ncfile.variables['upe2_bin'       ][:]); upe2_bin        = np.delete(upe2_bin       , ignored_points, axis = 0)
bpe2_bin        = np.copy(ncfile.variables['bpe2_bin'       ][:]); bpe2_bin        = np.delete(bpe2_bin       , ignored_points, axis = 0)
upa2_bin        = np.copy(ncfile.variables['upa2_bin'       ][:]); upa2_bin        = np.delete(upa2_bin       , ignored_points, axis = 0)
bpa2_bin        = np.copy(ncfile.variables['bpa2_bin'       ][:]); bpa2_bin        = np.delete(bpa2_bin       , ignored_points, axis = 0)
ux2_bin         = np.copy(ncfile.variables['ux2_bin'        ][:]); ux2_bin         = np.delete(ux2_bin        , ignored_points, axis = 0)
uy2_bin         = np.copy(ncfile.variables['uy2_bin'        ][:]); uy2_bin         = np.delete(uy2_bin        , ignored_points, axis = 0)
bx2_bin         = np.copy(ncfile.variables['bx2_bin'        ][:]); bx2_bin         = np.delete(bx2_bin        , ignored_points, axis = 0)
by2_bin         = np.copy(ncfile.variables['by2_bin'        ][:]); by2_bin         = np.delete(by2_bin        , ignored_points, axis = 0)
p_aw_bin        = np.copy(ncfile.variables['p_aw_bin'       ][:]); p_aw_bin        = np.delete(p_aw_bin       , ignored_points, axis = 0)
p_compr_bin     = np.copy(ncfile.variables['p_compr_bin'    ][:]); p_compr_bin     = np.delete(p_compr_bin    , ignored_points, axis = 0)
dissip_aw_bin   = np.copy(ncfile.variables['dissip_aw_bin'   ][:]); dissip_aw_bin    = np.delete(dissip_aw_bin   , ignored_points, axis = 0)
dissip_compr_bin= np.copy(ncfile.variables['dissip_compr_bin'][:]); dissip_compr_bin = np.delete(dissip_compr_bin, ignored_points, axis = 0)
ntrans_upe_upe_l_bin = np.copy(ncfile.variables['ntrans_upe_upe_l_bin'][:]); ntrans_upe_upe_l_bin = np.delete(ntrans_upe_upe_l_bin, ignored_points, axis = 0)
ntrans_bpe_upe_l_bin = np.copy(ncfile.variables['ntrans_bpe_upe_l_bin'][:]); ntrans_bpe_upe_l_bin = np.delete(ntrans_bpe_upe_l_bin, ignored_points, axis = 0)
ntrans_bpe_bpe_l_bin = np.copy(ncfile.variables['ntrans_bpe_bpe_l_bin'][:]); ntrans_bpe_bpe_l_bin = np.delete(ntrans_bpe_bpe_l_bin, ignored_points, axis = 0)
ntrans_upe_bpe_l_bin = np.copy(ncfile.variables['ntrans_upe_bpe_l_bin'][:]); ntrans_upe_bpe_l_bin = np.delete(ntrans_upe_bpe_l_bin, ignored_points, axis = 0)
ntrans_upa_upa_l_bin = np.copy(ncfile.variables['ntrans_upa_upa_l_bin'][:]); ntrans_upa_upa_l_bin = np.delete(ntrans_upa_upa_l_bin, ignored_points, axis = 0)
ntrans_bpa_upa_l_bin = np.copy(ncfile.variables['ntrans_bpa_upa_l_bin'][:]); ntrans_bpa_upa_l_bin = np.delete(ntrans_bpa_upa_l_bin, ignored_points, axis = 0)
ntrans_bpa_bpa_l_bin = np.copy(ncfile.variables['ntrans_bpa_bpa_l_bin'][:]); ntrans_bpa_bpa_l_bin = np.delete(ntrans_bpa_bpa_l_bin, ignored_points, axis = 0)
ntrans_upa_bpa_l_bin = np.copy(ncfile.variables['ntrans_upa_bpa_l_bin'][:]); ntrans_upa_bpa_l_bin = np.delete(ntrans_upa_bpa_l_bin, ignored_points, axis = 0)
ntrans_upe_upe_g_bin = np.copy(ncfile.variables['ntrans_upe_upe_g_bin'][:]); ntrans_upe_upe_g_bin = np.delete(ntrans_upe_upe_g_bin, ignored_points, axis = 0)
ntrans_bpe_upe_g_bin = np.copy(ncfile.variables['ntrans_bpe_upe_g_bin'][:]); ntrans_bpe_upe_g_bin = np.delete(ntrans_bpe_upe_g_bin, ignored_points, axis = 0)
ntrans_bpe_bpe_g_bin = np.copy(ncfile.variables['ntrans_bpe_bpe_g_bin'][:]); ntrans_bpe_bpe_g_bin = np.delete(ntrans_bpe_bpe_g_bin, ignored_points, axis = 0)
ntrans_upe_bpe_g_bin = np.copy(ncfile.variables['ntrans_upe_bpe_g_bin'][:]); ntrans_upe_bpe_g_bin = np.delete(ntrans_upe_bpe_g_bin, ignored_points, axis = 0)
ntrans_upa_upa_g_bin = np.copy(ncfile.variables['ntrans_upa_upa_g_bin'][:]); ntrans_upa_upa_g_bin = np.delete(ntrans_upa_upa_g_bin, ignored_points, axis = 0)
ntrans_bpa_upa_g_bin = np.copy(ncfile.variables['ntrans_bpa_upa_g_bin'][:]); ntrans_bpa_upa_g_bin = np.delete(ntrans_bpa_upa_g_bin, ignored_points, axis = 0)
ntrans_bpa_bpa_g_bin = np.copy(ncfile.variables['ntrans_bpa_bpa_g_bin'][:]); ntrans_bpa_bpa_g_bin = np.delete(ntrans_bpa_bpa_g_bin, ignored_points, axis = 0)
ntrans_upa_bpa_g_bin = np.copy(ncfile.variables['ntrans_upa_bpa_g_bin'][:]); ntrans_upa_bpa_g_bin = np.delete(ntrans_upa_bpa_g_bin, ignored_points, axis = 0)

zpep2_bin       = np.copy(ncfile.variables['zpep2_bin'      ][:]); zpep2_bin       = np.delete(zpep2_bin      , ignored_points, axis = 0)
zpem2_bin       = np.copy(ncfile.variables['zpem2_bin'      ][:]); zpem2_bin       = np.delete(zpem2_bin      , ignored_points, axis = 0)
zpap2_bin       = np.copy(ncfile.variables['zpap2_bin'      ][:]); zpap2_bin       = np.delete(zpap2_bin      , ignored_points, axis = 0)
zpam2_bin       = np.copy(ncfile.variables['zpam2_bin'      ][:]); zpam2_bin       = np.delete(zpam2_bin      , ignored_points, axis = 0)

ncfile.close()

# Load 2D cut
ncfile = netcdf.netcdf_file(input_dir+runname+'.out.2D.nc'+restart_num, 'r')   

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

upa_r_z0        = np.copy(ncfile.variables['upa_r_z0'       ][:]); upa_r_z0        = np.delete(upa_r_z0       , ignored_points_fld, axis = 0)
upa_r_x0        = np.copy(ncfile.variables['upa_r_x0'       ][:]); upa_r_x0        = np.delete(upa_r_x0       , ignored_points_fld, axis = 0)
upa_r_y0        = np.copy(ncfile.variables['upa_r_y0'       ][:]); upa_r_y0        = np.delete(upa_r_y0       , ignored_points_fld, axis = 0)

bpa_r_z0        = np.copy(ncfile.variables['bpa_r_z0'       ][:]); bpa_r_z0        = np.delete(bpa_r_z0       , ignored_points_fld, axis = 0)
bpa_r_x0        = np.copy(ncfile.variables['bpa_r_x0'       ][:]); bpa_r_x0        = np.delete(bpa_r_x0       , ignored_points_fld, axis = 0)
bpa_r_y0        = np.copy(ncfile.variables['bpa_r_y0'       ][:]); bpa_r_y0        = np.delete(bpa_r_y0       , ignored_points_fld, axis = 0)

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

tlab  = r'\Omega t'
zlab  = r'z\Omega/v_\rmA'
kzlab = r'k_z v_\rmA /\Omega'
