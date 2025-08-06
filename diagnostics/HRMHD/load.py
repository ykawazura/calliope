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
sgm = np.copy(ncfile.variables['sgm'].data)

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
p_xhl_sum       = np.copy(ncfile.variables['p_xhl_sum'      ][:]); p_xhl_sum       = np.delete(p_xhl_sum      , ignored_points, axis = 0)
zppe2_sum       = np.copy(ncfile.variables['zppe2_sum'      ][:]); zppe2_sum       = np.delete(zppe2_sum      , ignored_points, axis = 0)
zmpe2_sum       = np.copy(ncfile.variables['zmpe2_sum'      ][:]); zmpe2_sum       = np.delete(zmpe2_sum      , ignored_points, axis = 0)
zppa2_sum       = np.copy(ncfile.variables['zppa2_sum'      ][:]); zppa2_sum       = np.delete(zppa2_sum      , ignored_points, axis = 0)
zmpa2_sum       = np.copy(ncfile.variables['zmpa2_sum'      ][:]); zmpa2_sum       = np.delete(zmpa2_sum      , ignored_points, axis = 0)
upe2_KAW_dissip_sum  = np.copy(ncfile.variables['upe2_KAW_dissip_sum' ][:]); upe2_KAW_dissip_sum  = np.delete(upe2_KAW_dissip_sum , ignored_points, axis = 0)
bpe2_KAW_dissip_sum  = np.copy(ncfile.variables['bpe2_KAW_dissip_sum' ][:]); bpe2_KAW_dissip_sum  = np.delete(bpe2_KAW_dissip_sum , ignored_points, axis = 0)
upa2_KAW_dissip_sum  = np.copy(ncfile.variables['upa2_KAW_dissip_sum' ][:]); upa2_KAW_dissip_sum  = np.delete(upa2_KAW_dissip_sum , ignored_points, axis = 0)
bpa2_KAW_dissip_sum  = np.copy(ncfile.variables['bpa2_KAW_dissip_sum' ][:]); bpa2_KAW_dissip_sum  = np.delete(bpa2_KAW_dissip_sum , ignored_points, axis = 0)
upe2_ICW_dissip_sum  = np.copy(ncfile.variables['upe2_ICW_dissip_sum' ][:]); upe2_ICW_dissip_sum  = np.delete(upe2_ICW_dissip_sum , ignored_points, axis = 0)
bpe2_ICW_dissip_sum  = np.copy(ncfile.variables['bpe2_ICW_dissip_sum' ][:]); bpe2_ICW_dissip_sum  = np.delete(bpe2_ICW_dissip_sum , ignored_points, axis = 0)
upa2_ICW_dissip_sum  = np.copy(ncfile.variables['upa2_ICW_dissip_sum' ][:]); upa2_ICW_dissip_sum  = np.delete(upa2_ICW_dissip_sum , ignored_points, axis = 0)
bpa2_ICW_dissip_sum  = np.copy(ncfile.variables['bpa2_ICW_dissip_sum' ][:]); bpa2_ICW_dissip_sum  = np.delete(bpa2_ICW_dissip_sum , ignored_points, axis = 0)

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
dissip_KAW_bin  = np.copy(ncfile.variables['dissip_KAW_bin'  ][:]); dissip_KAW_bin   = np.delete(dissip_KAW_bin  , ignored_points, axis = 0)
dissip_ICW_bin  = np.copy(ncfile.variables['dissip_ICW_bin'  ][:]); dissip_ICW_bin   = np.delete(dissip_ICW_bin  , ignored_points, axis = 0)
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

zppe2_bin       = np.copy(ncfile.variables['zppe2_bin'      ][:]); zppe2_bin       = np.delete(zppe2_bin      , ignored_points, axis = 0)
zmpe2_bin       = np.copy(ncfile.variables['zmpe2_bin'      ][:]); zmpe2_bin       = np.delete(zmpe2_bin      , ignored_points, axis = 0)
zppa2_bin       = np.copy(ncfile.variables['zppa2_bin'      ][:]); zppa2_bin       = np.delete(zppa2_bin      , ignored_points, axis = 0)
zmpa2_bin       = np.copy(ncfile.variables['zmpa2_bin'      ][:]); zmpa2_bin       = np.delete(zmpa2_bin      , ignored_points, axis = 0)

upe2_KAW_bin    = np.copy(ncfile.variables['upe2_KAW_bin'   ][:]); upe2_KAW_bin    = np.delete(upe2_KAW_bin   , ignored_points, axis = 0)
bpe2_KAW_bin    = np.copy(ncfile.variables['bpe2_KAW_bin'   ][:]); bpe2_KAW_bin    = np.delete(bpe2_KAW_bin   , ignored_points, axis = 0)
upa2_KAW_bin    = np.copy(ncfile.variables['upa2_KAW_bin'   ][:]); upa2_KAW_bin    = np.delete(upa2_KAW_bin   , ignored_points, axis = 0)
bpa2_KAW_bin    = np.copy(ncfile.variables['bpa2_KAW_bin'   ][:]); bpa2_KAW_bin    = np.delete(bpa2_KAW_bin   , ignored_points, axis = 0)
upe2_ICW_bin    = np.copy(ncfile.variables['upe2_ICW_bin'   ][:]); upe2_ICW_bin    = np.delete(upe2_ICW_bin   , ignored_points, axis = 0)
bpe2_ICW_bin    = np.copy(ncfile.variables['bpe2_ICW_bin'   ][:]); bpe2_ICW_bin    = np.delete(bpe2_ICW_bin   , ignored_points, axis = 0)
upa2_ICW_bin    = np.copy(ncfile.variables['upa2_ICW_bin'   ][:]); upa2_ICW_bin    = np.delete(upa2_ICW_bin   , ignored_points, axis = 0)
bpa2_ICW_bin    = np.copy(ncfile.variables['bpa2_ICW_bin'   ][:]); bpa2_ICW_bin    = np.delete(bpa2_ICW_bin   , ignored_points, axis = 0)

ncfile.close()

tlab  = r'(v_\mathrm{A}/L_\|) t'
zlab  = r'z/L_\|'
kzlab = r'k_z L_\|'
