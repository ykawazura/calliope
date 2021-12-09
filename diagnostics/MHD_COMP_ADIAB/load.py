# -*- coding: utf-8 -*-
import warnings
warnings.filterwarnings('ignore')

import numpy as np
from scipy.io import netcdf

####### time index for final cut #########
input_dir = '../../'
runname = 'calliope'
restart_num = ''
final_idx   = -1

if final_idx != -1:
  print ('\n!!! CAUTION: final_idx = %d !!!\n' % final_idx)

####### movie or final cut #########
ismovie = False

####### ignore these points #########
ignored_points = [0]


####### load netcdf file #########
ncfile = netcdf.netcdf_file(input_dir+runname+'.out.nc'+restart_num, 'r')   

# Load parameters
shear_flg = np.copy(ncfile.variables['shear_flg'].data)
q         = np.copy(ncfile.variables['q'].data)
beta0     = np.copy(ncfile.variables['beta0'].data)
gamma     = np.copy(ncfile.variables['gamma'].data)

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
wprs_sum         = np.copy(ncfile.variables['wprs_sum'       ][:]); wprs_sum        = np.delete(wprs_sum       , ignored_points, axis = 0)
wkin_dot_sum     = np.copy(ncfile.variables['wkin_dot_sum'   ][:]); wkin_dot_sum    = np.delete(wkin_dot_sum   , ignored_points, axis = 0)
wmag_dot_sum     = np.copy(ncfile.variables['wmag_dot_sum'   ][:]); wmag_dot_sum    = np.delete(wmag_dot_sum   , ignored_points, axis = 0)
wprs_dot_sum     = np.copy(ncfile.variables['wprs_dot_sum'   ][:]); wprs_dot_sum    = np.delete(wprs_dot_sum   , ignored_points, axis = 0)
wkin_dissip_sum  = np.copy(ncfile.variables['wkin_dissip_sum'][:]); wkin_dissip_sum = np.delete(wkin_dissip_sum, ignored_points, axis = 0)
wmag_dissip_sum  = np.copy(ncfile.variables['wmag_dissip_sum'][:]); wmag_dissip_sum = np.delete(wmag_dissip_sum, ignored_points, axis = 0)
p_u_sum          = np.copy(ncfile.variables['p_u_sum'        ][:]); p_u_sum         = np.delete(p_u_sum        , ignored_points, axis = 0)
zp2_sum          = np.copy(ncfile.variables['zp2_sum'        ][:]); zp2_sum         = np.delete(zp2_sum        , ignored_points, axis = 0)
zm2_sum          = np.copy(ncfile.variables['zm2_sum'        ][:]); zm2_sum         = np.delete(zm2_sum        , ignored_points, axis = 0)
smach_rms        = np.copy(ncfile.variables['smach_rms'      ][:]); smach_rms       = np.delete(smach_rms      , ignored_points, axis = 0)
amach_rms        = np.copy(ncfile.variables['amach_rms'      ][:]); amach_rms       = np.delete(amach_rms      , ignored_points, axis = 0)
beta_rms         = np.copy(ncfile.variables['beta_rms'       ][:]); beta_rms        = np.delete(beta_rms       , ignored_points, axis = 0)

# mean magnetic field
b0       = np.copy(ncfile.variables['b0'           ][:]); b0            = np.delete(b0           , ignored_points, axis = 0)

# Load binned spectra
u2_bin   = np.copy(ncfile.variables['u2_bin'       ][:]); u2_bin        = np.delete(u2_bin       , ignored_points, axis = 0)
b2_bin   = np.copy(ncfile.variables['b2_bin'       ][:]); b2_bin        = np.delete(b2_bin       , ignored_points, axis = 0)
rho_bin  = np.copy(ncfile.variables['rho_bin'      ][:]); rho_bin       = np.delete(rho_bin      , ignored_points, axis = 0)
sgm_bin  = np.copy(ncfile.variables['sgm_bin'      ][:]); sgm_bin       = np.delete(sgm_bin      , ignored_points, axis = 0)
zp2_bin  = np.copy(ncfile.variables['zp2_bin'      ][:]); zp2_bin       = np.delete(zp2_bin      , ignored_points, axis = 0)
zm2_bin  = np.copy(ncfile.variables['zm2_bin'      ][:]); zm2_bin       = np.delete(zm2_bin      , ignored_points, axis = 0)

# Load 2D cut
rho_r_z0 = np.copy(ncfile.variables['rho_r_z0'][:]); rho_r_z0 = np.delete(rho_r_z0, ignored_points, axis = 0)
rho_r_x0 = np.copy(ncfile.variables['rho_r_x0'][:]); rho_r_x0 = np.delete(rho_r_x0, ignored_points, axis = 0)
rho_r_y0 = np.copy(ncfile.variables['rho_r_y0'][:]); rho_r_y0 = np.delete(rho_r_y0, ignored_points, axis = 0)

mx_r_z0  = np.copy(ncfile.variables['mx_r_z0' ][:]); mx_r_z0  = np.delete(mx_r_z0 , ignored_points, axis = 0)
mx_r_x0  = np.copy(ncfile.variables['mx_r_x0' ][:]); mx_r_x0  = np.delete(mx_r_x0 , ignored_points, axis = 0)
mx_r_y0  = np.copy(ncfile.variables['mx_r_y0' ][:]); mx_r_y0  = np.delete(mx_r_y0 , ignored_points, axis = 0)
                                                                                  
my_r_z0  = np.copy(ncfile.variables['my_r_z0' ][:]); my_r_z0  = np.delete(my_r_z0 , ignored_points, axis = 0)
my_r_x0  = np.copy(ncfile.variables['my_r_x0' ][:]); my_r_x0  = np.delete(my_r_x0 , ignored_points, axis = 0)
my_r_y0  = np.copy(ncfile.variables['my_r_y0' ][:]); my_r_y0  = np.delete(my_r_y0 , ignored_points, axis = 0)
                                                                                  
mz_r_z0  = np.copy(ncfile.variables['mz_r_z0' ][:]); mz_r_z0  = np.delete(mz_r_z0 , ignored_points, axis = 0)
mz_r_x0  = np.copy(ncfile.variables['mz_r_x0' ][:]); mz_r_x0  = np.delete(mz_r_x0 , ignored_points, axis = 0)
mz_r_y0  = np.copy(ncfile.variables['mz_r_y0' ][:]); mz_r_y0  = np.delete(mz_r_y0 , ignored_points, axis = 0)
                                                                                  
wx_r_z0  = np.copy(ncfile.variables['wx_r_z0' ][:]); wx_r_z0  = np.delete(wx_r_z0 , ignored_points, axis = 0)
wx_r_x0  = np.copy(ncfile.variables['wx_r_x0' ][:]); wx_r_x0  = np.delete(wx_r_x0 , ignored_points, axis = 0)
wx_r_y0  = np.copy(ncfile.variables['wx_r_y0' ][:]); wx_r_y0  = np.delete(wx_r_y0 , ignored_points, axis = 0)
                                                                                  
wy_r_z0  = np.copy(ncfile.variables['wy_r_z0' ][:]); wy_r_z0  = np.delete(wy_r_z0 , ignored_points, axis = 0)
wy_r_x0  = np.copy(ncfile.variables['wy_r_x0' ][:]); wy_r_x0  = np.delete(wy_r_x0 , ignored_points, axis = 0)
wy_r_y0  = np.copy(ncfile.variables['wy_r_y0' ][:]); wy_r_y0  = np.delete(wy_r_y0 , ignored_points, axis = 0)
                                                                                  
wz_r_z0  = np.copy(ncfile.variables['wz_r_z0' ][:]); wz_r_z0  = np.delete(wz_r_z0 , ignored_points, axis = 0)
wz_r_x0  = np.copy(ncfile.variables['wz_r_x0' ][:]); wz_r_x0  = np.delete(wz_r_x0 , ignored_points, axis = 0)
wz_r_y0  = np.copy(ncfile.variables['wz_r_y0' ][:]); wz_r_y0  = np.delete(wz_r_y0 , ignored_points, axis = 0)
                                                                                  
bx_r_z0  = np.copy(ncfile.variables['bx_r_z0' ][:]); bx_r_z0  = np.delete(bx_r_z0 , ignored_points, axis = 0)
bx_r_x0  = np.copy(ncfile.variables['bx_r_x0' ][:]); bx_r_x0  = np.delete(bx_r_x0 , ignored_points, axis = 0)
bx_r_y0  = np.copy(ncfile.variables['bx_r_y0' ][:]); bx_r_y0  = np.delete(bx_r_y0 , ignored_points, axis = 0)
                                                                                  
by_r_z0  = np.copy(ncfile.variables['by_r_z0' ][:]); by_r_z0  = np.delete(by_r_z0 , ignored_points, axis = 0)
by_r_x0  = np.copy(ncfile.variables['by_r_x0' ][:]); by_r_x0  = np.delete(by_r_x0 , ignored_points, axis = 0)
by_r_y0  = np.copy(ncfile.variables['by_r_y0' ][:]); by_r_y0  = np.delete(by_r_y0 , ignored_points, axis = 0)
                                                                                  
bz_r_z0  = np.copy(ncfile.variables['bz_r_z0' ][:]); bz_r_z0  = np.delete(bz_r_z0 , ignored_points, axis = 0)
bz_r_x0  = np.copy(ncfile.variables['bz_r_x0' ][:]); bz_r_x0  = np.delete(bz_r_x0 , ignored_points, axis = 0)
bz_r_y0  = np.copy(ncfile.variables['bz_r_y0' ][:]); bz_r_y0  = np.delete(bz_r_y0 , ignored_points, axis = 0)
                                                                                  
jx_r_z0  = np.copy(ncfile.variables['jx_r_z0' ][:]); jx_r_z0  = np.delete(jx_r_z0 , ignored_points, axis = 0)
jx_r_x0  = np.copy(ncfile.variables['jx_r_x0' ][:]); jx_r_x0  = np.delete(jx_r_x0 , ignored_points, axis = 0)
jx_r_y0  = np.copy(ncfile.variables['jx_r_y0' ][:]); jx_r_y0  = np.delete(jx_r_y0 , ignored_points, axis = 0)
                                                                                  
jy_r_z0  = np.copy(ncfile.variables['jy_r_z0' ][:]); jy_r_z0  = np.delete(jy_r_z0 , ignored_points, axis = 0)
jy_r_x0  = np.copy(ncfile.variables['jy_r_x0' ][:]); jy_r_x0  = np.delete(jy_r_x0 , ignored_points, axis = 0)
jy_r_y0  = np.copy(ncfile.variables['jy_r_y0' ][:]); jy_r_y0  = np.delete(jy_r_y0 , ignored_points, axis = 0)
                                                                                  
jz_r_z0  = np.copy(ncfile.variables['jz_r_z0' ][:]); jz_r_z0  = np.delete(jz_r_z0 , ignored_points, axis = 0)
jz_r_x0  = np.copy(ncfile.variables['jz_r_x0' ][:]); jz_r_x0  = np.delete(jz_r_x0 , ignored_points, axis = 0)
jz_r_y0  = np.copy(ncfile.variables['jz_r_y0' ][:]); jz_r_y0  = np.delete(jz_r_y0 , ignored_points, axis = 0)

sgm_r_z0 = np.copy(ncfile.variables['sgm_r_z0'][:]); sgm_r_z0 = np.delete(sgm_r_z0, ignored_points, axis = 0)
sgm_r_x0 = np.copy(ncfile.variables['sgm_r_x0'][:]); sgm_r_x0 = np.delete(sgm_r_x0, ignored_points, axis = 0)
sgm_r_y0 = np.copy(ncfile.variables['sgm_r_y0'][:]); sgm_r_y0 = np.delete(sgm_r_y0, ignored_points, axis = 0)

ux_r_z0  = mx_r_z0/rho_r_z0 
ux_r_x0  = mx_r_x0/rho_r_x0 
ux_r_y0  = mx_r_y0/rho_r_y0 
             
uy_r_z0  = my_r_z0/rho_r_z0
uy_r_x0  = my_r_x0/rho_r_x0 
uy_r_y0  = my_r_y0/rho_r_y0 
              
uz_r_z0  = mz_r_z0/rho_r_z0
uz_r_x0  = mz_r_x0/rho_r_x0 
uz_r_y0  = mz_r_y0/rho_r_y0 

prs_r_z0 = sgm_r_z0/rho_r_z0**(-gamma + 1.0)
prs_r_x0 = sgm_r_x0/rho_r_x0**(-gamma + 1.0) 
prs_r_y0 = sgm_r_y0/rho_r_y0**(-gamma + 1.0) 

ncfile.close()

tlab  = r'$t$'
xlab  = r'$x$'
ylab  = r'$y$'
zlab  = r'$z$'
kxlab = r'$k_x$'
kylab = r'$k_y$'
kzlab = r'$k_z$'
kplab = r'$k_\+$'
