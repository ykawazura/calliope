# -*- coding: utf-8 -*-
import warnings
warnings.filterwarnings('ignore')

import numpy as np
from scipy.io import netcdf

####### time index for final cut #########
input_dir = '../../'
runname = 'calliope'
restart_num = ''
final_idx      = -1
final_fld_idx  = -1
final_kpar_idx = -1
final_SF2_idx  = -1

if final_idx != -1:
  print ('\n!!! CAUTION: final_idx = %d !!!\n' % final_idx)

####### movie or final cut #########
ismovie = False

####### ignore these points #########
ignored_points      = [0]
ignored_points_fld  = [ ]
ignored_points_kpar = [ ]
ignored_points_SF2  = [ ]


####### load netcdf file #########
ncfile = netcdf.netcdf_file(input_dir+runname+'.out.nc'+restart_num, 'r')   

# Load parameters
de = np.copy(ncfile.variables['de'].data)

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
wmag_sum         = np.copy(ncfile.variables['wmag_sum'        ][:]); wmag_sum         = np.delete(wmag_sum        , ignored_points, axis = 0)
wmag_dot_sum     = np.copy(ncfile.variables['wmag_dot_sum'    ][:]); wmag_dot_sum     = np.delete(wmag_dot_sum    , ignored_points, axis = 0)
wmag_dissip_sum  = np.copy(ncfile.variables['wmag_dissip_sum' ][:]); wmag_dissip_sum  = np.delete(wmag_dissip_sum , ignored_points, axis = 0)
p_ext_sum        = np.copy(ncfile.variables['p_ext_sum'       ][:]); p_ext_sum        = np.delete(p_ext_sum       , ignored_points, axis = 0)

# mean magnetic field
b0            = np.copy(ncfile.variables['b0'           ][:]); b0            = np.delete(b0           , ignored_points, axis = 0)

# Load binned spectra
b2_bin        = np.copy(ncfile.variables['b2_bin'       ][:]); b2_bin        = np.delete(b2_bin       , ignored_points, axis = 0)
bx2_bin       = np.copy(ncfile.variables['bx2_bin'      ][:]); bx2_bin       = np.delete(bx2_bin      , ignored_points, axis = 0)
by2_bin       = np.copy(ncfile.variables['by2_bin'      ][:]); by2_bin       = np.delete(by2_bin      , ignored_points, axis = 0)
bz2_bin       = np.copy(ncfile.variables['bz2_bin'      ][:]); bz2_bin       = np.delete(bz2_bin      , ignored_points, axis = 0)

ncfile.close()

# Load kpar
ncfile = netcdf.netcdf_file(input_dir+runname+'.out.kpar.nc'+restart_num, 'r')   
try:
  tt_kpar   = np.copy(ncfile.variables['tt'       ][:]); tt_kpar   = np.delete(tt_kpar  , ignored_points_kpar, axis = 0)
  kpar_b    = np.copy(ncfile.variables['kpar_b'   ][:]); kpar_b    = np.delete(kpar_b   , ignored_points_kpar, axis = 0)
  kpar_u    = np.copy(ncfile.variables['kpar_u'   ][:]); kpar_u    = np.delete(kpar_u   , ignored_points_kpar, axis = 0)
  b1_ovr_b0 = np.copy(ncfile.variables['b1_ovr_b0'][:]); b1_ovr_b0 = np.delete(b1_ovr_b0, ignored_points_kpar, axis = 0)
except KeyError:
  pass
ncfile.close()

# Load SF2
ncfile = netcdf.netcdf_file(input_dir+runname+'.out.SF2.nc'+restart_num, 'r')   
try:
  tt_SF2 = np.copy(ncfile.variables['tt'  ][:]); tt_SF2  = np.delete(tt_SF2 , ignored_points_SF2, axis = 0)
  lpar   = np.copy(ncfile.variables['lpar'][:])
  lper   = np.copy(ncfile.variables['lper'][:])
  nl     = lpar.size
  SF2b   = np.copy(ncfile.variables['SF2b'][:]);  SF2b = np.delete(SF2b, ignored_points_SF2, axis = 0)
  SF2u   = np.copy(ncfile.variables['SF2u'][:]);  SF2u = np.delete(SF2u, ignored_points_SF2, axis = 0)
except KeyError:
  pass


ncfile.close()

tlab  = r'$\Omega_\rmi t$'
xlab  = r'$x/d_\rmi$'
ylab  = r'$y/d_\rmi$'
zlab  = r'$z/d_\rmi$'
kxlab = r'$k_x d_\rmi$'
kylab = r'$k_y d_\rmi$'
kzlab = r'$k_z d_\rmi$'
kplab = r'$k_\+ d_\rmi$'

# # Load series modes
# import os

# label_series = {'bx': r'$B_x$', 'by': r'$B_y$', 'bz': r'$B_z$'}

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

