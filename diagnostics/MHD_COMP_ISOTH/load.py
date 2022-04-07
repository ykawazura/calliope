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
ignored_points     = [ ]
ignored_points_fld = [ ]


####### load netcdf file #########
ncfile = netcdf.netcdf_file(input_dir+runname+'.out.nc'+restart_num, 'r')   

# Load parameters
shear_flg = np.copy(ncfile.variables['shear_flg'].data)
q         = np.copy(ncfile.variables['q'].data)
beta0     = np.copy(ncfile.variables['beta0'].data)

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
wkin_sum         = np.copy(ncfile.variables['wkin_sum'       ][:]); wkin_sum        = np.delete(wkin_sum       , ignored_points, axis = 0)
wmag_sum         = np.copy(ncfile.variables['wmag_sum'       ][:]); wmag_sum        = np.delete(wmag_sum       , ignored_points, axis = 0)
wrho_sum         = np.copy(ncfile.variables['wrho_sum'       ][:]); wrho_sum        = np.delete(wrho_sum       , ignored_points, axis = 0)
wkin_dot_sum     = np.copy(ncfile.variables['wkin_dot_sum'   ][:]); wkin_dot_sum    = np.delete(wkin_dot_sum   , ignored_points, axis = 0)
wmag_dot_sum     = np.copy(ncfile.variables['wmag_dot_sum'   ][:]); wmag_dot_sum    = np.delete(wmag_dot_sum   , ignored_points, axis = 0)
wrho_dot_sum     = np.copy(ncfile.variables['wrho_dot_sum'   ][:]); wrho_dot_sum    = np.delete(wrho_dot_sum   , ignored_points, axis = 0)
wkin_dissip_sum  = np.copy(ncfile.variables['wkin_dissip_sum'][:]); wkin_dissip_sum = np.delete(wkin_dissip_sum, ignored_points, axis = 0)
wmag_dissip_sum  = np.copy(ncfile.variables['wmag_dissip_sum'][:]); wmag_dissip_sum = np.delete(wmag_dissip_sum, ignored_points, axis = 0)
wrho_dissip_sum  = np.copy(ncfile.variables['wrho_dissip_sum'][:]); wrho_dissip_sum = np.delete(wrho_dissip_sum, ignored_points, axis = 0)
p_ext_sum        = np.copy(ncfile.variables['p_ext_sum'      ][:]); p_ext_sum       = np.delete(p_ext_sum      , ignored_points, axis = 0)
p_re_sum         = np.copy(ncfile.variables['p_re_sum'       ][:]); p_re_sum        = np.delete(p_re_sum       , ignored_points, axis = 0)
p_ma_sum         = np.copy(ncfile.variables['p_ma_sum'       ][:]); p_ma_sum        = np.delete(p_ma_sum       , ignored_points, axis = 0)
zp2_sum          = np.copy(ncfile.variables['zp2_sum'        ][:]); zp2_sum         = np.delete(zp2_sum        , ignored_points, axis = 0)
zm2_sum          = np.copy(ncfile.variables['zm2_sum'        ][:]); zm2_sum         = np.delete(zm2_sum        , ignored_points, axis = 0)
smach_rms        = np.copy(ncfile.variables['smach_rms'      ][:]); smach_rms       = np.delete(smach_rms      , ignored_points, axis = 0)
amach_rms        = np.copy(ncfile.variables['amach_rms'      ][:]); amach_rms       = np.delete(amach_rms      , ignored_points, axis = 0)
beta_rms         = np.copy(ncfile.variables['beta_rms'       ][:]); beta_rms        = np.delete(beta_rms       , ignored_points, axis = 0)

# mean magnetic field
b0       = np.copy(ncfile.variables['b0'           ][:]); b0            = np.delete(b0           , ignored_points, axis = 0)

# Load binned spectra
u2_bin           = np.copy(ncfile.variables['u2_bin'         ][:]); u2_bin          = np.delete(u2_bin         , ignored_points, axis = 0)
ux2_bin          = np.copy(ncfile.variables['ux2_bin'        ][:]); ux2_bin         = np.delete(ux2_bin        , ignored_points, axis = 0)
uy2_bin          = np.copy(ncfile.variables['uy2_bin'        ][:]); uy2_bin         = np.delete(uy2_bin        , ignored_points, axis = 0)
uz2_bin          = np.copy(ncfile.variables['uz2_bin'        ][:]); uz2_bin         = np.delete(uz2_bin        , ignored_points, axis = 0)
b2_bin           = np.copy(ncfile.variables['b2_bin'         ][:]); b2_bin          = np.delete(b2_bin         , ignored_points, axis = 0)
bx2_bin          = np.copy(ncfile.variables['bx2_bin'        ][:]); bx2_bin         = np.delete(bx2_bin        , ignored_points, axis = 0)
by2_bin          = np.copy(ncfile.variables['by2_bin'        ][:]); by2_bin         = np.delete(by2_bin        , ignored_points, axis = 0)
bz2_bin          = np.copy(ncfile.variables['bz2_bin'        ][:]); bz2_bin         = np.delete(bz2_bin        , ignored_points, axis = 0)
rho_bin          = np.copy(ncfile.variables['rho_bin'        ][:]); rho_bin         = np.delete(rho_bin        , ignored_points, axis = 0)
zp2_bin          = np.copy(ncfile.variables['zp2_bin'        ][:]); zp2_bin         = np.delete(zp2_bin        , ignored_points, axis = 0)
zm2_bin          = np.copy(ncfile.variables['zm2_bin'        ][:]); zm2_bin         = np.delete(zm2_bin        , ignored_points, axis = 0)

ncfile.close()

tlab  = r'$t$'
xlab  = r'$x$'
ylab  = r'$y$'
zlab  = r'$z$'
kxlab = r'$k_x$'
kylab = r'$k_y$'
kzlab = r'$k_z$'
kplab = r'$k_\+$'

# # Load series modes
# import os

# label_series = {'rho': r'$\rho$', 'mx': r'$m_x$', 'my': r'$m_y$', 'mz': r'$m_z$', 'bx': r'$B_x$', 'by': r'$B_y$', 'bz': r'$B_z$'}

# filename = input_dir+runname+'.modes.out'+restart_num
# if os.path.isfile(filename):
  # data = np.loadtxt(filename, usecols=[0,2,3,4,5,6])
  # name = np.loadtxt(filename, usecols=[1], dtype = "unicode")

  # fldnames_unique = np.unique(name)
  # nfld = fldnames_unique.size

  # tt_unique = np.unique(data.T[0])
  # ntt = tt_unique.size

  # modes_unique = np.unique((data.T[1:4]).T, axis=0)
  # nmodes = modes_unique.shape[0]

  # tt_ = {}
  # f_  = {}
  # for n in fldnames_unique:
    # tt_[n] = np.zeros([ntt, nmodes])
    # f_ [n] = np.zeros([ntt, nmodes], dtype=complex)

  # for n, d in zip(name, data):
    # time = d[0]
    # mode = d[1:4]
    # time_idx = np.argwhere(tt_unique == time)[0][0]
    # mode_idx = np.argwhere([np.array_equal(x, mode) for x in modes_unique])[0][0]

    # tt_[n][time_idx, mode_idx] = time
    # f_ [n][time_idx, mode_idx] = d[4] + 1j*d[5]
