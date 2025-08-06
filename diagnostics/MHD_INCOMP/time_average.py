# -*- coding: utf-8 -*-
from load import *
from fft import *
import sys
sys.path.append('../')
from plots import *
from scipy.integrate import simps, trapz
from scipy import interpolate

avg_start_time = 1.
avg_end_time   = 900.

@jit
def time_average(x, y, axis=0): # x: 1D array, y: any-D array 
  print ('averaging from %.2E to %.2E' % (x[0], x[-1]))
  return trapz(y, x, axis=axis)/(x[-1] - x[0])
  # return trapz(y, x, axis=axis)/(x[-1] - x[0])

@jit
def std_dev(x, y): # x, y: 1D array
  y_intp = interpolate.interp1d(x, y)
  return np.std(y_intp(np.linspace(x[0], x[-1], x.size)), ddof=0)

@jit
def cov(x, y1, y2): # x, y: 1D array
  y1_intp  = interpolate.interp1d(x, y1)
  y2_intp  = interpolate.interp1d(x, y2)
  y1_intp_ = y1_intp(np.linspace(x[0], x[-1], x.size))
  y2_intp_ = y2_intp(np.linspace(x[0], x[-1], x.size))
  
  return np.cov( np.stack((y1_intp_, y2_intp_), axis=0) )

##########################################################
#              average energy time evolution             #
##########################################################
print('\nplotting energy\n')
outdir = './fig_energy/'

W        = u2_sum + b2_sum
W_dot    = u2dot_sum + b2dot_sum
D        = u2dissip_sum + b2dissip_sum
P        = p_ext_sum + p_re_sum + p_ma_sum

avg_start = np.argmin(abs(tt - avg_start_time))
avg_end   = np.argmin(abs(tt - avg_end_time))

if avg_end == -1:
  avg_end = np.argwhere(tt == tt[avg_end])[0][0]

W_avg       = time_average(tt[avg_start:avg_end], W     [avg_start:avg_end], axis=0)
W_dot_avg   = time_average(tt[avg_start:avg_end], W_dot [avg_start:avg_end], axis=0)
D_avg       = time_average(tt[avg_start:avg_end], D     [avg_start:avg_end], axis=0)
P_avg       = time_average(tt[avg_start:avg_end], P     [avg_start:avg_end], axis=0)
u2_sum_avg  = time_average(tt[avg_start:avg_end], u2_sum[avg_start:avg_end], axis=0)
b2_sum_avg  = time_average(tt[avg_start:avg_end], b2_sum[avg_start:avg_end], axis=0)

W_err       = std_dev(tt[avg_start:avg_end], W     [avg_start:avg_end])
W_dot_err   = std_dev(tt[avg_start:avg_end], W_dot [avg_start:avg_end])
D_err       = std_dev(tt[avg_start:avg_end], D     [avg_start:avg_end])
P_err       = std_dev(tt[avg_start:avg_end], P     [avg_start:avg_end])
u2_sum_err  = std_dev(tt[avg_start:avg_end], u2_sum[avg_start:avg_end])
b2_sum_err  = std_dev(tt[avg_start:avg_end], b2_sum[avg_start:avg_end])

s =     'average over t     = [%.3E' % tt[avg_start] + ', %.3E' % tt[avg_end] + ']' + '\n'
s = s + 'average over index = [' + str(avg_start) + ', ' + str(avg_end) + ']' + '\n'
s = s + ' error of ratio is calculated by (a + da)/(b + db) ~ a/b*[1 + sqrt( (da/a)^2 + (db/b)^2 - 2cov/(a*b))]' + '\n\n' 
s = s + '  W     = %.3E \pm %.3E'  % (W_avg    , W_err    ) + '\n'
s = s + '  W_dot = %.3E \pm %.3E'  % (W_dot_avg, W_dot_err) + '\n'
s = s + '  P     = %.3E \pm %.3E'  % (P_avg    , P_err    ) + '\n'
s = s + '  D     = %.3E \pm %.3E'  % (D_avg    , D_err    ) + '\n'
print (s)

f = open('time_average.txt','w') 
f.write(s) 
f.close() 


ys = [ 
       u2dot_sum + b2dot_sum, 
       u2dissip_sum, 
       b2dissip_sum, 
       -p_ext_sum, 
       -p_re_sum, 
       -p_ma_sum, 
       u2dot_sum + b2dot_sum + u2dissip_sum + b2dissip_sum - p_ext_sum - p_re_sum - p_ma_sum,
       np.full([tt[avg_start:avg_end].size], W_dot_avg),
       np.full([tt[avg_start:avg_end].size], D_avg),
       np.full([tt[avg_start:avg_end].size], -P_avg),
     ]
xs = [
       tt,
       tt,
       tt,
       tt,
       tt,
       tt,
       tt,
       tt[avg_start:avg_end],
       tt[avg_start:avg_end],
       tt[avg_start:avg_end],
     ]
ls = [ 
        '', 
        '', 
        '', 
        '', 
        '', 
        '', 
        'k--', 
        '', 
        '', 
        '', 
     ]
legends = [ 
       r'$\mathrm{d}W/\mathrm{d} t$', 
       r'$D_u$', 
       r'$D_B$', 
       r'$P_\mr{ext}$', 
       r'$P_\mr{Re}$', 
       r'$P_\mr{M}$', 
       r'balance', 
       '', 
       '', 
       '', 
     ]
plot_1d_many_average(xs, ys, tt[avg_start], tt[avg_end], xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'balance_all_avg.pdf')


ys = [ 
       u2_sum, 
       b2_sum, 
       np.full([tt[avg_start:avg_end].size], u2_sum_avg),
       np.full([tt[avg_start:avg_end].size], b2_sum_avg),
     ]
ls = [ 
        '', 
        '', 
        '', 
        '', 
     ]
xs = [
       tt,
       tt,
       tt[avg_start:avg_end],
       tt[avg_start:avg_end],
     ]
legends = [ 
       r'$W_{u}$', 
       r'$W_{B}$', 
       '',
       '',
     ]
plot_1d_many_average(xs, ys, tt[avg_start], tt[avg_end], xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'energy_all_avg.pdf')

##########################################################
#                    average kspectrum                   #
##########################################################
print('\nplotting kspectrum\n')
outdir = './fig_kspectrum/'

u2_bin_avg   = time_average(tt[avg_start:avg_end], u2_bin  [avg_start:avg_end], axis=0)
ux2_bin_avg  = time_average(tt[avg_start:avg_end], ux2_bin [avg_start:avg_end], axis=0)
uy2_bin_avg  = time_average(tt[avg_start:avg_end], uy2_bin [avg_start:avg_end], axis=0)
uz2_bin_avg  = time_average(tt[avg_start:avg_end], uz2_bin [avg_start:avg_end], axis=0)
b2_bin_avg   = time_average(tt[avg_start:avg_end], b2_bin  [avg_start:avg_end], axis=0)
bx2_bin_avg  = time_average(tt[avg_start:avg_end], bx2_bin [avg_start:avg_end], axis=0)
by2_bin_avg  = time_average(tt[avg_start:avg_end], by2_bin [avg_start:avg_end], axis=0)
bz2_bin_avg  = time_average(tt[avg_start:avg_end], bz2_bin [avg_start:avg_end], axis=0)
zp2_bin_avg  = time_average(tt[avg_start:avg_end], zp2_bin [avg_start:avg_end], axis=0)
zm2_bin_avg  = time_average(tt[avg_start:avg_end], zm2_bin [avg_start:avg_end], axis=0)
p_re_bin_avg = time_average(tt[avg_start:avg_end], p_re_bin[avg_start:avg_end], axis=0)
p_ma_bin_avg = time_average(tt[avg_start:avg_end], p_ma_bin[avg_start:avg_end], axis=0)

if nlz == nkz:
  kp_end = np.argmin(np.abs(kpbin - kpbin.max()*2./3.))
else:
  kp_end = kpbin.size - 1
kp_log_kpar_end = kpbin_log_kpar.size - 1

# If shear is on, add a vertical line indicating the fastest-MRI mode
if shear_flg == 1:
  b00 = np.max(b0[0])
  k_mri = np.sqrt(15.)/4./b00

def add_vertical_line(xs, ys, ls, legends):
  xs.append([k_mri, k_mri])
  ys.append([np.min([u2_bin_avg[1:kp_end].min(), b2_bin_avg[1:kp_end].min()]), 
             np.max([u2_bin_avg[1:kp_end].max(), b2_bin_avg[1:kp_end].max()])])
  ls.append('k:')
  legends.append(r'$k = k_\mr{MRI} = (\sqrt{15}/4)\varpi_0/v_\rmA$')

  return xs, ys, ls, legends

def fold_in_kz(u):
  unew = [(u[:, i, :] + u[:, -i, :])/2 for i in np.arange(0, int(nkz/2)+1)]
  return np.transpose(np.asarray(unew), axes=[1,0,2])

#--------------------------------------------------------#
#                      plot 1D spectra                   #
#--------------------------------------------------------#
ys = [ 
       u2_bin_avg[1:kp_end], 
       b2_bin_avg[1:kp_end], 
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*b2_bin_avg[1:kp_end][0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*b2_bin_avg[1:kp_end][0]
     ]
xs = [ 
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end], 
      kpbin[1:kp_end]  
     ]
ls = [ 
        '', 
        '', 
        'k--', 
        'k-.', 
     ]
legends = [ 
            r'$E_{u}$', 
            r'$E_{B}$',
            r'-5/3',
            r'-3/2',
          ]

if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'k_spectra_avg.pdf')


# u by components
ys = [ 
       ux2_bin_avg[1:kp_end], 
       uy2_bin_avg[1:kp_end], 
       uz2_bin_avg[1:kp_end], 
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*u2_bin_avg[1:kp_end][0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*u2_bin_avg[1:kp_end][0]
     ]
xs = [ 
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end], 
      kpbin[1:kp_end]  
     ]
ls = [ 
        '', 
        '', 
        '', 
        'k--', 
        'k-.', 
     ]
legends = [ 
            r'$E_{u_x}$', 
            r'$E_{u_y}$', 
            r'$E_{u_z}$', 
            r'-5/3',
            r'-3/2',
          ]

if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'k_spectra_u_avg.pdf')


# B by components
ys = [ 
       bx2_bin_avg[1:kp_end], 
       by2_bin_avg[1:kp_end], 
       bz2_bin_avg[1:kp_end], 
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*b2_bin_avg[1:kp_end][0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*b2_bin_avg[1:kp_end][0]
     ]
xs = [ 
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end], 
      kpbin[1:kp_end]  
     ]
ls = [ 
        '', 
        '', 
        '', 
        'k--', 
        'k-.', 
     ]
legends = [ 
            r'$E_{B_x}$', 
            r'$E_{B_y}$', 
            r'$E_{B_z}$', 
            r'-5/3',
            r'-3/2',
          ]

if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'k_spectra_b_avg.pdf')


# Elsasser fields
ys = [ 
       zp2_bin_avg[1:kp_end], 
       zm2_bin_avg[1:kp_end], 
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*zp2_bin_avg[1:kp_end][0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*zp2_bin_avg[1:kp_end][0],
     ]
xs = [ 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
     ]
ls = [ 
        '', 
        '', 
        'k--', 
        'k-.', 
     ]
legends = [ 
            r'$E_{Z^+}$', 
            r'$E_{Z^-}$',
            r'-5/3',
            r'-3/2',
          ]

if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'k_spectra_ELS_avg.pdf')


# MRI Injection
ys = [ 
       p_re_bin_avg[1:kp_end], 
       p_ma_bin_avg[1:kp_end], 
     ]
xs = [ 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
     ]
ls = [ 
        '', 
        '', 
     ]
legends = [ 
            r'$P_\mr{Re}$', 
            r'$P_\mr{M}$',
          ]

if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'k_MRI_avg.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ek_avg.txt'  , np.column_stack((kpbin[:kp_end], 
                                                  u2_bin_avg [:kp_end],
                                                 ux2_bin_avg [:kp_end],
                                                 uy2_bin_avg [:kp_end],
                                                 uz2_bin_avg [:kp_end],
                                                  b2_bin_avg [:kp_end],
                                                 bx2_bin_avg [:kp_end],
                                                 by2_bin_avg [:kp_end],
                                                 bz2_bin_avg [:kp_end],
                                                 zp2_bin_avg [:kp_end],
                                                 zm2_bin_avg [:kp_end],
                                                 p_re_bin_avg[:kp_end],
                                                 p_ma_bin_avg[:kp_end],
                                                 )), fmt='%E')

# kpar
if tt_kpar.size > 0:
  avg_start = np.argmin(abs(tt_kpar - avg_start_time))
  avg_end   = np.argmin(abs(tt_kpar - avg_end_time))

  if avg_end == -1:
    avg_end = np.argwhere(tt == tt_kpar[-1])[0][0]

  kpar_b_avg    = time_average(tt_kpar[avg_start:avg_end], kpar_b   [avg_start:avg_end], axis=0)
  kpar_u_avg    = time_average(tt_kpar[avg_start:avg_end], kpar_u   [avg_start:avg_end], axis=0)
  b1_ovr_b0_avg = time_average(tt_kpar[avg_start:avg_end], b1_ovr_b0[avg_start:avg_end], axis=0)

  b1par2_avg    = time_average(tt_kpar[avg_start:avg_end], b1par2   [avg_start:avg_end], axis=0)
  b1prp2_avg    = time_average(tt_kpar[avg_start:avg_end], b1prp2   [avg_start:avg_end], axis=0)
  u1par2_avg    = time_average(tt_kpar[avg_start:avg_end], u1par2   [avg_start:avg_end], axis=0)
  u1prp2_avg    = time_average(tt_kpar[avg_start:avg_end], u1prp2   [avg_start:avg_end], axis=0)

  ys = [ 
         kpar_b_avg[1:kp_log_kpar_end], 
         kpar_u_avg[1:kp_log_kpar_end], 
         kpbin_log_kpar[1:kp_log_kpar_end]**(2./3.)/kpbin_log_kpar[1]**(2./3.)*kpar_b_avg[1:kp_log_kpar_end][0],
       ]
  xs = [ 
        kpbin_log_kpar[1:kp_log_kpar_end], 
        kpbin_log_kpar[1:kp_log_kpar_end], 
        kpbin_log_kpar[1:kp_log_kpar_end], 
       ]
  ls = [ 
          '', 
          '', 
          'k--', 
       ]
  legends = [ 
              r'$k_\|^b$', 
              r'$k_\|^u$', 
              r'2/3',
            ]

  if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

  plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'kpar_avg.pdf')


  # delta b/b0
  ys = [ 
         b1_ovr_b0_avg[1:kp_log_kpar_end], 
       ]
  xs = [ 
        kpbin_log_kpar[1:kp_log_kpar_end], 
       ]
  ls = [ 
          '', 
       ]
  legends = [ 
              r'$\delta b/b_0$', 
            ]

  if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

  plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'b1_ovr_b0_avg.pdf')


  # shear AWs vs pseudo AWs
  ys = [ 
         b1par2_avg[1:kp_log_kpar_end], 
         b1prp2_avg[1:kp_log_kpar_end], 
         (b1par2_avg/b1prp2_avg)[1:kp_log_kpar_end], 
         u1par2_avg[1:kp_log_kpar_end], 
         u1prp2_avg[1:kp_log_kpar_end], 
         (u1par2_avg/u1prp2_avg)[1:kp_log_kpar_end], 
       ]
  xs = [ 
        kpbin_log_kpar[1:kp_log_kpar_end], 
        kpbin_log_kpar[1:kp_log_kpar_end], 
        kpbin_log_kpar[1:kp_log_kpar_end], 
        kpbin_log_kpar[1:kp_log_kpar_end], 
        kpbin_log_kpar[1:kp_log_kpar_end], 
        kpbin_log_kpar[1:kp_log_kpar_end], 
       ]
  ls = [ 
          'r-', 
          'r--', 
          'r:', 
          'b-', 
          'b--', 
          'b:', 
       ]
  legends = [ 
              r'$\delta B_\|^2$', 
              r'$\delta B_\+^2$', 
              r'$\delta B_\|^2/\delta B_\+^2$', 
              r'$u_\|^2$', 
              r'$u_\+^2$', 
              r'$u_\|^2/u_\+^2$', 
            ]

  if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

  plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'shearAW_vs_pseudoAW_avg.pdf')

  #------------------#
  #   output ascii   #
  #------------------#
  np.savetxt(outdir + 'shearAW_vs_pseudoAW_avg.txt'  , np.column_stack((kpbin_log_kpar[:kp_log_kpar_end], 
                                                       b1par2_avg[:kp_log_kpar_end],
                                                       b1prp2_avg[:kp_log_kpar_end],
                                                       u1par2_avg[:kp_log_kpar_end],
                                                       u1prp2_avg[:kp_log_kpar_end],
                                                       )), fmt='%E')

  np.savetxt(outdir + 'kpar_avg.txt'  , np.column_stack((kpbin_log_kpar[:kp_log_kpar_end], 
                                                   kpar_b_avg   [:kp_log_kpar_end],
                                                   kpar_u_avg   [:kp_log_kpar_end],
                                                   b1_ovr_b0_avg[:kp_log_kpar_end],
                                                   )), fmt='%E')

#--------------------------------------------------------#
#                       2D spectra                       #
#--------------------------------------------------------#
# summed over kx or ky or kz
tt_2d = np.transpose(np.loadtxt(input_dir+'out2d'+restart_num+'/time.dat'))[0]
nt_2d = tt_2d.size
if nt_2d == 1 : tt_2d = [tt_2d]

avg_start = np.argmin(abs(tt_2d - avg_start_time))
avg_end   = np.argmin(abs(tt_2d - avg_end_time))

nkx_mid = np.where(np.sign(kx) == -1)[0][0]
nkz_mid = np.where(np.sign(kz) == -1)[0][0]
kx_ = np.concatenate([kx[nkx_mid:], kx[0:nkx_mid]])
kz_ = np.concatenate([kz[nkz_mid:], kz[0:nkz_mid]])

u2_kxy_sum_kz = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/u2_kxy_sum_kz.dat').reshape(nt_2d, nky, nkx)); u2_kxy_sum_kz = np.transpose(u2_kxy_sum_kz, axes=(2, 1, 0));  u2_kxy_sum_kz = np.delete(u2_kxy_sum_kz, ignored_points_2d, axis = 0); u2_kxy_sum_kz = np.concatenate([u2_kxy_sum_kz[:, :, nkx_mid:], u2_kxy_sum_kz[:, :, 0:nkx_mid]], axis=2)
u2_kyz_sum_kx = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/u2_kyz_sum_kx.dat').reshape(nt_2d, nkz, nky)); u2_kyz_sum_kx = np.transpose(u2_kyz_sum_kx, axes=(2, 1, 0));  u2_kyz_sum_kx = np.delete(u2_kyz_sum_kx, ignored_points_2d, axis = 0); u2_kyz_sum_kx = np.concatenate([u2_kyz_sum_kx[:, nkz_mid:, :], u2_kyz_sum_kx[:, 0:nkz_mid, :]], axis=1)
u2_kxz_sum_ky = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/u2_kxz_sum_ky.dat').reshape(nt_2d, nkz, nkx)); u2_kxz_sum_ky = np.transpose(u2_kxz_sum_ky, axes=(2, 1, 0));  u2_kxz_sum_ky = np.delete(u2_kxz_sum_ky, ignored_points_2d, axis = 0); u2_kxz_sum_ky = np.concatenate([u2_kxz_sum_ky[:, nkz_mid:, :], u2_kxz_sum_ky[:, 0:nkz_mid, :]], axis=1); u2_kxz_sum_ky = np.concatenate([u2_kxz_sum_ky[:, :, nkx_mid:], u2_kxz_sum_ky[:, :, 0:nkx_mid]], axis=2)
b2_kxy_sum_kz = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/b2_kxy_sum_kz.dat').reshape(nt_2d, nky, nkx)); b2_kxy_sum_kz = np.transpose(b2_kxy_sum_kz, axes=(2, 1, 0));  b2_kxy_sum_kz = np.delete(b2_kxy_sum_kz, ignored_points_2d, axis = 0); b2_kxy_sum_kz = np.concatenate([b2_kxy_sum_kz[:, :, nkx_mid:], b2_kxy_sum_kz[:, :, 0:nkx_mid]], axis=2)
b2_kyz_sum_kx = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/b2_kyz_sum_kx.dat').reshape(nt_2d, nkz, nky)); b2_kyz_sum_kx = np.transpose(b2_kyz_sum_kx, axes=(2, 1, 0));  b2_kyz_sum_kx = np.delete(b2_kyz_sum_kx, ignored_points_2d, axis = 0); b2_kyz_sum_kx = np.concatenate([b2_kyz_sum_kx[:, nkz_mid:, :], b2_kyz_sum_kx[:, 0:nkz_mid, :]], axis=1)
b2_kxz_sum_ky = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/b2_kxz_sum_ky.dat').reshape(nt_2d, nkz, nkx)); b2_kxz_sum_ky = np.transpose(b2_kxz_sum_ky, axes=(2, 1, 0));  b2_kxz_sum_ky = np.delete(b2_kxz_sum_ky, ignored_points_2d, axis = 0); b2_kxz_sum_ky = np.concatenate([b2_kxz_sum_ky[:, nkz_mid:, :], b2_kxz_sum_ky[:, 0:nkz_mid, :]], axis=1); b2_kxz_sum_ky = np.concatenate([b2_kxz_sum_ky[:, :, nkx_mid:], b2_kxz_sum_ky[:, :, 0:nkx_mid]], axis=2)

u2_kxy_sum_kz_avg = time_average(tt_2d[avg_start:avg_end], u2_kxy_sum_kz[avg_start:avg_end,:,:], axis=0)
u2_kyz_sum_kx_avg = time_average(tt_2d[avg_start:avg_end], u2_kyz_sum_kx[avg_start:avg_end,:,:], axis=0)
u2_kxz_sum_ky_avg = time_average(tt_2d[avg_start:avg_end], u2_kxz_sum_ky[avg_start:avg_end,:,:], axis=0)
b2_kxy_sum_kz_avg = time_average(tt_2d[avg_start:avg_end], b2_kxy_sum_kz[avg_start:avg_end,:,:], axis=0)
b2_kyz_sum_kx_avg = time_average(tt_2d[avg_start:avg_end], b2_kyz_sum_kx[avg_start:avg_end,:,:], axis=0)
b2_kxz_sum_ky_avg = time_average(tt_2d[avg_start:avg_end], b2_kxz_sum_ky[avg_start:avg_end,:,:], axis=0)

plot_2d(np.log10(u2_kxy_sum_kz_avg[1:-1,1:-1]), kx_[1:-1], ky [1:-1], xlab=r'$k_x$', ylab=r'$k_y$', title=r'$E_u$', contour=True, aspect='auto', save=outdir+'u2_kxy_avg.pdf')
plot_2d(np.log10(u2_kyz_sum_kx_avg[1:-1,1:-1]), ky [1:-1], kz_[1:-1], xlab=r'$k_y$', ylab=r'$k_z$', title=r'$E_u$', contour=True, aspect='auto', save=outdir+'u2_kyz_avg.pdf')
plot_2d(np.log10(u2_kxz_sum_ky_avg[1:-1,1:-1]), kx_[1:-1], kz_[1:-1], xlab=r'$k_x$', ylab=r'$k_z$', title=r'$E_u$', contour=True, aspect='auto', save=outdir+'u2_kxz_avg.pdf')
plot_2d(np.log10(b2_kxy_sum_kz_avg[1:-1,1:-1]), kx_[1:-1], ky [1:-1], xlab=r'$k_x$', ylab=r'$k_y$', title=r'$E_B$', contour=True, aspect='auto', save=outdir+'b2_kxy_avg.pdf')
plot_2d(np.log10(b2_kyz_sum_kx_avg[1:-1,1:-1]), ky [1:-1], kz_[1:-1], xlab=r'$k_y$', ylab=r'$k_z$', title=r'$E_B$', contour=True, aspect='auto', save=outdir+'b2_kyz_avg.pdf')
plot_2d(np.log10(b2_kxz_sum_ky_avg[1:-1,1:-1]), kx_[1:-1], kz_[1:-1], xlab=r'$k_x$', ylab=r'$k_z$', title=r'$E_B$', contour=True, aspect='auto', save=outdir+'b2_kxz_avg.pdf')


# bin over (kx, ky) leaving kz or (kx, kz) leaving ky
u2_kxy_vs_kz = np.fromfile(input_dir+'out2d'+restart_num+'/u2_kxy_vs_kz.dat').reshape(nt_2d, nkz, nkpolar); u2_kxy_vs_kz = fold_in_kz(u2_kxy_vs_kz)
u2_kxz_vs_ky = np.fromfile(input_dir+'out2d'+restart_num+'/u2_kxz_vs_ky.dat').reshape(nt_2d, nky, nkpolar)
b2_kxy_vs_kz = np.fromfile(input_dir+'out2d'+restart_num+'/b2_kxy_vs_kz.dat').reshape(nt_2d, nkz, nkpolar); b2_kxy_vs_kz = fold_in_kz(b2_kxy_vs_kz)
b2_kxz_vs_ky = np.fromfile(input_dir+'out2d'+restart_num+'/b2_kxz_vs_ky.dat').reshape(nt_2d, nky, nkpolar)

u2_kxy_vs_kz_avg = time_average(tt_2d[avg_start:avg_end], u2_kxy_vs_kz[avg_start:avg_end,:,:], axis=0)
u2_kxz_vs_ky_avg = time_average(tt_2d[avg_start:avg_end], u2_kxz_vs_ky[avg_start:avg_end,:,:], axis=0)
b2_kxy_vs_kz_avg = time_average(tt_2d[avg_start:avg_end], b2_kxy_vs_kz[avg_start:avg_end,:,:], axis=0)
b2_kxz_vs_ky_avg = time_average(tt_2d[avg_start:avg_end], b2_kxz_vs_ky[avg_start:avg_end,:,:], axis=0)

plot_log2d(u2_kxy_vs_kz_avg[1:int(nkz/2),1:kp_end], kpbin[1:kp_end], kz[1:int(nkz/2)] , xlab=r'$k_\+$', ylab=r'$k_z$', title=r'$E_u$', contour=True, save=outdir+'u2_kxy_vs_kz_avg.pdf')
plot_log2d(b2_kxy_vs_kz_avg[1:int(nkz/2),1:kp_end], kpbin[1:kp_end], kz[1:int(nkz/2)] , xlab=r'$k_\+$', ylab=r'$k_z$', title=r'$E_u$', contour=True, save=outdir+'b2_kxy_vs_kz_avg.pdf')

plot_log2d(u2_kxz_vs_ky_avg[1:nky - 1   ,1:kp_end], kpbin[1:kp_end], ky[1:nky - 1   ] , xlab=r'$k_\+$', ylab=r'$k_y$', title=r'$E_B$', contour=True, save=outdir+'u2_kxz_vs_ky_avg.pdf')
plot_log2d(b2_kxz_vs_ky_avg[1:nky - 1   ,1:kp_end], kpbin[1:kp_end], ky[1:nky - 1   ] , xlab=r'$k_\+$', ylab=r'$k_y$', title=r'$E_B$', contour=True, save=outdir+'b2_kxz_vs_ky_avg.pdf')

#------------------#
#   output data    #
#------------------#
from scipy.io import savemat

savemat(outdir + '2d_bin_in_3d_grid' , {'kpbin':kpbin, 'ky':ky, 'kz':kz[:int(nkz/2)+1]})
savemat(outdir + '2d_bin_in_3d_avg' , {
                            'u2_kxy_vs_kz' :u2_kxy_vs_kz_avg[:,:],
                            'u2_kxz_vs_ky' :u2_kxz_vs_ky_avg[:,:],
                            'b2_kxy_vs_kz' :b2_kxy_vs_kz_avg[:,:],
                            'b2_kxz_vs_ky' :b2_kxz_vs_ky_avg[:,:],
                          })


##########################################################
#                    average energy SF2                  #
##########################################################
# if tt_SF2.size > 0:
  # print('\nplotting SF2\n')
  # outdir = './fig_SF2/'

  # avg_start = np.argmin(abs(tt_SF2 - avg_start_time))
  # avg_end   = np.argmin(abs(tt_SF2 - avg_end_time))

  # if avg_end == -1:
    # avg_end = np.argwhere(tt == tt_SF2[-1])[0][0]

  # SF2b_avg = time_average(tt_SF2[avg_start:avg_end], SF2b[avg_start:avg_end,:,:], axis=0)
  # SF2u_avg = time_average(tt_SF2[avg_start:avg_end], SF2u[avg_start:avg_end,:,:], axis=0)

  # plot_SF2(SF2b_avg[:,:], SF2u_avg[:,:], lpar, lper, xlab=r'$\ell_\|$', ylab=r'$\ell_\+$', title=r'', cmp=parula_map, save=outdir+'SF2_avg.pdf')

  # ys = [ 
         # SF2b_avg[1:, 0 ], 
         # SF2u_avg[1:, 0 ], 
         # SF2b_avg[0 , 1:], 
         # SF2u_avg[0 , 1:], 
         # lper[1:]**(2./3.)/lper[1]**(2./3.)*SF2b_avg[1, 0],
         # lpar[1:]**(3./3.)/lpar[1]**(3./3.)*SF2b_avg[0, 1],
       # ]
  # xs = [ 
        # lpar[1:], 
        # lpar[1:], 
        # lpar[1:], 
        # lpar[1:], 
        # lper[1:], 
        # lper[1:], 
       # ]
  # ls = [ 
          # '', 
          # '', 
          # '', 
          # '', 
          # 'k--', 
          # 'k--', 
       # ]
  # legends = [ 
              # r'$\mr{SF}^2_b(0,\ell_\+)$',
              # r'$\mr{SF}^2_u(0,\ell_\+)$',
              # r'$\mr{SF}^2_b(\ell_\|,0)$', 
              # r'$\mr{SF}^2_u(\ell_\|,0)$', 
              # r'2/3',
              # r'1',
            # ]
  # plot_log1d_many(xs, ys, xlab=r'$\ell_\|$ or $\ell_\+$', legends=legends, ls=ls, legendloc='upper left', ylab='', term=True, save=outdir+'SF2_0_avg.pdf')

  # # kpar(kper)
  # SF2b_0_lper = SF2b_avg[: , 0] # SF2(0, lper)
  # SF2u_0_lper = SF2u_avg[: , 0] # SF2(0, lper)
  # SF2b_lpar_0 = SF2b_avg[0 , :] # SF2(lpar, 0)
  # SF2u_lpar_0 = SF2u_avg[0 , :] # SF2(lpar, 0)

  # idx_b = [np.argmin(abs(x - SF2b_lpar_0[1:])) for x in SF2b_0_lper[1:]] # where SF2(lpar, 0) = SF2(0, lper)
  # idx_u = [np.argmin(abs(x - SF2u_lpar_0[1:])) for x in SF2u_0_lper[1:]] # where SF2(lpar, 0) = SF2(0, lper)

  # kper   = 2*np.pi/lper[1:]
  # kpar_b = 2*np.pi/lpar[idx_b]
  # kpar_u = 2*np.pi/lpar[idx_u]

  # ys = [ 
        # kpar_b, 
        # kpar_b, 
        # kper**(2./3.), 
        # kper**(1./2.), 
       # ]
  # xs = [ 
        # kper, 
        # kper, 
        # kper, 
        # kper, 
       # ]
  # ls = [ 
          # '', 
          # '', 
          # 'k--', 
          # 'k:', 
       # ]
  # legends = [ 
              # '$b$',
              # '$u$',
              # r'2/3',
              # r'1/2',
            # ]
  # plot_log1d_many(xs, ys, xlab=r'$k_\+$', ylab=r'$k_\|$', legends=legends, ls=ls, legendloc='upper left', title='', term=True, save=outdir+'kpar_avg.pdf')

  # #------------------#
  # #    output data   #
  # #------------------#
  # np.savetxt(outdir + 'kpar_avg.txt'  , np.column_stack((kper, kpar_b, kpar_u)), fmt='%E')

  # from scipy.io import savemat
  # savemat(outdir + 'grid', {'lper':lper, 'lpar':lpar})
  # savemat(outdir + 'SF2_avg' , {
                             # 'SF2b' :SF2b_avg[:,:],
                             # 'SF2u' :SF2u_avg[:,:],
                           # })

##########################################################
#             average shell-to-shell transfer            #
##########################################################
if tt_nltrans.size > 0:
  print('\nplotting shell-to-shell transfer function\n')
  outdir = './fig_nltrans/'

  k0idx = 3

  kp       = kpbin_log_nltrans[k0idx:]
  trans_uu = trans_uu[:, k0idx:, k0idx:]
  trans_bb = trans_bb[:, k0idx:, k0idx:]
  trans_ub = trans_ub[:, k0idx:, k0idx:]
  trans_bu = trans_bu[:, k0idx:, k0idx:]

  avg_start = np.argmin(abs(tt_nltrans - avg_start_time))
  avg_end   = np.argmin(abs(tt_nltrans - avg_end_time))

  if avg_end == -1:
    avg_end = np.argwhere(tt == tt_nltrans[-1])[0][0]

  trans_uu_avg = time_average(tt_nltrans[avg_start:avg_end], trans_uu[avg_start:avg_end,:,:], axis=0)
  trans_bb_avg = time_average(tt_nltrans[avg_start:avg_end], trans_bb[avg_start:avg_end,:,:], axis=0)
  trans_ub_avg = time_average(tt_nltrans[avg_start:avg_end], trans_ub[avg_start:avg_end,:,:], axis=0)
  trans_bu_avg = time_average(tt_nltrans[avg_start:avg_end], trans_bu[avg_start:avg_end,:,:], axis=0)

  #--------------------------------------------------------#
  #                     2D transfer map                    #
  #--------------------------------------------------------#
  fig, axes = plt.subplots(2, 2, figsize=(15, 12))

  transes = [trans_uu_avg[:,:], trans_bb_avg[:,:], trans_ub_avg[:,:], trans_bu_avg[:,:]]
  names   = [r'$\calT_{uu}$', r'$\calT_{BB}$', r'$\calT_{uB}$', r'$\calT_{Bu}$']

  for i, trans in enumerate(transes):
    col = i%2
    row = int(i/2)

    umax = np.max(abs(trans))

    ax = axes[row][col]
    X, Y = np.meshgrid(kp, kp)
    im = ax.pcolormesh(X, Y, trans.T, cmap='RdBu_r', vmin = -umax, vmax = umax)

    ax.set_title(names[i])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(im, cax=cax)
    cb.formatter.set_powerlimits((0, 0))
    cb.ax.yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
    cb.ax.yaxis.set_offset_position('left')  
    cb.locator = ticker.MaxNLocator(nbins=9)
    cb.update_ticks()

  for ax in axes.flatten():
    ax.set_aspect('equal')
    ax.minorticks_on()

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(kp[0], kp[-1])
    ax.set_ylim(kp[0], kp[-1])

  axes[0][0].tick_params(labelbottom=False)
  axes[0][1].tick_params(labelbottom=False)
  axes[0][1].tick_params(labelleft  =False)
  axes[1][1].tick_params(labelleft  =False)

  axes[0][0].set_ylabel(r'$K$')
  axes[1][0].set_ylabel(r'$K$')
  axes[1][0].set_xlabel(r'$Q$')
  axes[1][1].set_xlabel(r'$Q$')

  plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
  plt.savefig(outdir+'trans_avg.pdf')

  #------------------#
  #   output data    #
  #------------------#
  from scipy.io import savemat

  savemat(outdir + 'grid' , {'kpbin_log':kp})
  savemat(outdir + 'trans_avg' , {
                              'nltrans_uu' :trans_uu_avg[:,:],
                              'nltrans_bb' :trans_bb_avg[:,:],
                              'nltrans_ub' :trans_ub_avg[:,:],
                              'nltrans_bu' :trans_bu_avg[:,:],
                            })


  #--------------------------------------------------------#
  #                         Fluxes                         #
  #--------------------------------------------------------#
  flx_uu = np.asarray([np.sum(trans_uu_avg.T[i:, :i]) for i in np.arange(0, kp.size)])
  flx_bb = np.asarray([np.sum(trans_bb_avg.T[i:, :i]) for i in np.arange(0, kp.size)])
  flx_ub = np.asarray([np.sum(trans_ub_avg.T[i:, :i]) for i in np.arange(0, kp.size)])
  flx_bu = np.asarray([np.sum(trans_bu_avg.T[i:, :i]) for i in np.arange(0, kp.size)])

  ys = [ 
        flx_uu[:],
        flx_bb[:],
        flx_ub[:],
        flx_bu[:],
        (flx_uu[:] + flx_bb[:] + flx_bu[:] + flx_ub[:]),
       ]
  xs = [ 
        kp[:], 
        kp[:], 
        kp[:], 
        kp[:], 
        kp[:], 
       ]
  ls = [ 
          '', 
          '', 
          '', 
          '', 
          '', 
       ]
  legends = [ 
              r'$\Pi^{u^<}_{u^>}$', 
              r'$\Pi^{B^<}_{B^>}$', 
              r'$\Pi^{u^<}_{B^>}$', 
              r'$\Pi^{B^<}_{u^>}$', 
              r'$\Pi_\mr{tot}$', 
            ]

  plot_semilogx1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'flux_avg.pdf')
