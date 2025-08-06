# -*- coding: utf-8 -*-
from load import *
from fft import *
import sys
sys.path.append('../')
from plots import *
from scipy.integrate import simps, trapz
from scipy import interpolate

avg_start = 0
avg_end   = -1

avg_start_kpar    = 1
avg_end_kpar      = -1

@jit
def time_average(x, y, axis=0): # x: 1D array, y: any-D array 
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

if avg_end == -1:
  avg_end = np.argwhere(tt == tt[avg_end])[0][0]

##########################################################
#              average energy time evolution             #
##########################################################
print('\nplotting energy\n')
outdir = './fig_energy/'

W        = wmag_sum
W_dot    = wmag_dot_sum
D        = wmag_dissip_sum
P        = p_ext_sum

W_avg        = time_average(tt[avg_start:avg_end], W       [avg_start:avg_end], axis=0)
W_dot_avg    = time_average(tt[avg_start:avg_end], W_dot   [avg_start:avg_end], axis=0)
D_avg        = time_average(tt[avg_start:avg_end], D       [avg_start:avg_end], axis=0)
P_avg        = time_average(tt[avg_start:avg_end], P       [avg_start:avg_end], axis=0)
wmag_sum_avg = time_average(tt[avg_start:avg_end], wmag_sum[avg_start:avg_end], axis=0)
             
W_err        = std_dev(tt[avg_start:avg_end], W       [avg_start:avg_end])
W_dot_err    = std_dev(tt[avg_start:avg_end], W_dot   [avg_start:avg_end])
D_err        = std_dev(tt[avg_start:avg_end], D       [avg_start:avg_end])
P_err        = std_dev(tt[avg_start:avg_end], P       [avg_start:avg_end])
wmag_sum_err = std_dev(tt[avg_start:avg_end], wmag_sum[avg_start:avg_end])

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
       wmag_dot_sum, 
       wmag_dissip_sum, 
       -p_ext_sum, 
       wmag_dot_sum + wmag_dissip_sum - p_ext_sum,
       np.full([tt[avg_start:avg_end].size], W_dot_avg),
       np.full([tt[avg_start:avg_end].size], D_avg),
       np.full([tt[avg_start:avg_end].size], -P_avg),
     ]
xs = [
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
        'k--', 
        '', 
        '', 
        '', 
     ]
legends = [ 
       r'$\mathrm{d}W/\mathrm{d} t$', 
       r'$D_B$', 
       r'$-P$', 
       r'balance', 
       '', 
       '', 
       '', 
     ]
plot_1d_many_average(xs, ys, tt[avg_start], tt[avg_end], xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'balance_all_avg.pdf')


ys = [ 
       wmag_sum, 
       np.full([tt[avg_start:avg_end].size], wmag_sum_avg),
     ]
ls = [ 
        '', 
        '', 
     ]
xs = [
       tt,
       tt[avg_start:avg_end],
     ]
legends = [ 
       r'$W_{B}$', 
       '',
     ]
plot_1d_many_average(xs, ys, tt[avg_start], tt[avg_end], xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'energy_all_avg.pdf')

##########################################################
#                    average kspectrum                   #
##########################################################
print('\nplotting kspectrum\n')
outdir = './fig_kspectrum/'

b2_bin_avg  = time_average(tt[avg_start:avg_end], b2_bin [avg_start:avg_end], axis=0)
bx2_bin_avg = time_average(tt[avg_start:avg_end], bx2_bin[avg_start:avg_end], axis=0)
by2_bin_avg = time_average(tt[avg_start:avg_end], by2_bin[avg_start:avg_end], axis=0)
bz2_bin_avg = time_average(tt[avg_start:avg_end], bz2_bin[avg_start:avg_end], axis=0)

if nlz == nkz:
  kp_end = np.argmin(np.abs(kpbin - kpbin.max()*2./3.))
else:
  kp_end = kpbin.size - 1
  kp_log_end = kpbin_log.size - 1

#--------------------------------------------------------#
#                      plot 1D spectra                   #
#--------------------------------------------------------#
ys = [ 
       b2_bin_avg[1:kp_end], 
       kpbin[1:kp_end]**(-7./3.)/kpbin[1]**(-7./3.)*b2_bin_avg[1:kp_end][0],
     ]
xs = [ 
      kpbin[1:kp_end], 
      kpbin[1:kp_end]  
     ]
ls = [ 
        '', 
        'k-.', 
     ]
legends = [ 
            r'$E_{B}$',
            r'-7/3',
          ]

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'k_spectra_avg.pdf')


# B by components
ys = [ 
       bx2_bin_avg[1:kp_end], 
       by2_bin_avg[1:kp_end], 
       bz2_bin_avg[1:kp_end], 
       kpbin[1:kp_end]**(-7./3.)/kpbin[1]**(-7./3.)*b2_bin_avg[1:kp_end][0],
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
        '', 
        'k-.', 
     ]
legends = [ 
            r'$E_{B_x}$', 
            r'$E_{B_y}$', 
            r'$E_{B_z}$', 
            r'-7/3',
          ]

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'k_spectra_b_avg.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ek_avg.txt'  , np.column_stack((kpbin[:kp_end], 
                                                  b2_bin_avg[:kp_end],
                                                 bx2_bin_avg[:kp_end],
                                                 by2_bin_avg[:kp_end],
                                                 bz2_bin_avg[:kp_end],
                                                 )), fmt='%E')
# kpar
if tt_kpar.size > 0:
  avg_start = avg_start_kpar
  avg_end   = avg_end_kpar
  if avg_end == -1:
    avg_end = np.argwhere(tt == tt_kpar[-1])[0][0]

  kpar_b_avg    = time_average(tt_kpar[avg_start:avg_end], kpar_b   [avg_start:avg_end], axis=0)
  b1_ovr_b0_avg = time_average(tt_kpar[avg_start:avg_end], b1_ovr_b0[avg_start:avg_end], axis=0)

  b1par2_avg    = time_average(tt_kpar[avg_start:avg_end], b1par2   [avg_start:avg_end], axis=0)
  b1prp2_avg    = time_average(tt_kpar[avg_start:avg_end], b1prp2   [avg_start:avg_end], axis=0)

  ys = [ 
         kpar_b_avg[1:kp_log_end], 
         kpbin_log[1:kp_log_end]**(1./3.)/kpbin_log[1]**(1./3.)*kpar_b_avg[1:kp_log_end][0],
       ]
  xs = [ 
        kpbin_log[1:kp_log_end], 
        kpbin_log[1:kp_log_end], 
       ]
  ls = [ 
          '', 
          'k--', 
       ]
  legends = [ 
              r'$k_\|^b$', 
              r'1/3',
            ]

  plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'kpar_avg.pdf')


  # delta b/b0
  ys = [ 
         b1_ovr_b0_avg[1:kp_log_end], 
       ]
  xs = [ 
        kpbin_log[1:kp_log_end], 
       ]
  ls = [ 
          '', 
       ]
  legends = [ 
              r'$\delta b/b_0$', 
            ]

  plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'b1_ovr_b0_avg.pdf')


  # shear AWs vs pseudo AWs
  ys = [ 
         b1par2_avg[1:kp_log_end], 
         b1prp2_avg[1:kp_log_end], 
         (b1par2_avg/b1prp2_avg)[1:kp_log_end], 
       ]
  xs = [ 
        kpbin_log[1:kp_log_end], 
        kpbin_log[1:kp_log_end], 
        kpbin_log[1:kp_log_end], 
       ]
  ls = [ 
          'r-', 
          'r--', 
          'r:', 
       ]
  legends = [ 
              r'$\delta B_\|^2$', 
              r'$\delta B_\+^2$', 
              r'$\delta B_\|^2/\delta B_\+^2$', 
            ]

  plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'shearAW_vs_pseudoAW_avg.pdf')

  #------------------#
  #   output ascii   #
  #------------------#
  np.savetxt(outdir + 'shearAW_vs_pseudoAW_avg.txt'  , np.column_stack((kpbin_log[:kp_log_end], 
                                                       b1par2_avg[:kp_log_end],
                                                       b1prp2_avg[:kp_log_end],
                                                       )), fmt='%E')

  np.savetxt(outdir + 'kpar_avg.txt'  , np.column_stack((kpbin_log[:kp_log_end], 
                                                   kpar_b_avg   [:kp_log_end],
                                                   b1_ovr_b0_avg[:kp_log_end],
                                                   )), fmt='%E')


##########################################################
#                    average energy SF2                  #
##########################################################
print('\nplotting SF2\n')
outdir = './fig_SF2/'

SF2b_avg = time_average(tt_SF2[avg_start:avg_end], SF2b[avg_start:avg_end,:,:], axis=0)
SF2u_avg = time_average(tt_SF2[avg_start:avg_end], SF2u[avg_start:avg_end,:,:], axis=0)

plot_SF2(SF2b_avg[:,:], SF2u_avg[:,:], lpar, lper, xlab=r'$\ell_\|$', ylab=r'$\ell_\+$', title=r'', save=outdir+'SF2_avg.pdf')

ys = [ 
       SF2b_avg[1:, 0 ], 
       SF2u_avg[1:, 0 ], 
       SF2b_avg[0 , 1:], 
       SF2u_avg[0 , 1:], 
       lper[1:]**(2./3.)/lper[1]**(2./3.)*SF2b_avg[1, 0],
       lpar[1:]**(3./3.)/lpar[1]**(3./3.)*SF2b_avg[0, 1],
     ]
xs = [ 
      lpar[1:], 
      lpar[1:], 
      lpar[1:], 
      lpar[1:], 
      lper[1:], 
      lper[1:], 
     ]
ls = [ 
        '', 
        '', 
        '', 
        '', 
        'k--', 
        'k--', 
     ]
legends = [ 
            r'$\mr{SF}^2_b(0,\ell_\+)$',
            r'$\mr{SF}^2_u(0,\ell_\+)$',
            r'$\mr{SF}^2_b(\ell_\|,0)$', 
            r'$\mr{SF}^2_u(\ell_\|,0)$', 
            r'2/3',
            r'1',
          ]
plot_log1d_many(xs, ys, xlab=r'$\ell_\|$ or $\ell_\+$', legends=legends, ls=ls, legendloc='upper left', ylab='', term=True, save=outdir+'SF2_0_avg.pdf')

# kpar(kper)
SF2b_0_lper = SF2b_avg[: , 0] # SF2(0, lper)
SF2u_0_lper = SF2u_avg[: , 0] # SF2(0, lper)
SF2b_lpar_0 = SF2b_avg[0 , :] # SF2(lpar, 0)
SF2u_lpar_0 = SF2u_avg[0 , :] # SF2(lpar, 0)

idx_b = [np.argmin(abs(x - SF2b_lpar_0[1:])) for x in SF2b_0_lper[1:]] # where SF2(lpar, 0) = SF2(0, lper)
idx_u = [np.argmin(abs(x - SF2u_lpar_0[1:])) for x in SF2u_0_lper[1:]] # where SF2(lpar, 0) = SF2(0, lper)

kper   = 2*np.pi/lper[1:]
kpar_b = 2*np.pi/lpar[idx_b]
kpar_u = 2*np.pi/lpar[idx_u]

ys = [ 
      kpar_b, 
      kpar_b, 
      kper**(2./3.), 
      kper**(1./2.), 
     ]
xs = [ 
      kper, 
      kper, 
      kper, 
      kper, 
     ]
ls = [ 
        '', 
        '', 
        'k--', 
        'k:', 
     ]
legends = [ 
            '$b$',
            '$u$',
            r'2/3',
            r'1/2',
          ]
plot_log1d_many(xs, ys, xlab=r'$k_\+$', ylab=r'$k_\|$', legends=legends, ls=ls, legendloc='upper left', title='', term=True, save=outdir+'kpar_avg.pdf')

#------------------#
#    output data   #
#------------------#
np.savetxt(outdir + 'kpar_avg.txt'  , np.column_stack((kper, kpar_b, kpar_u)), fmt='%E')

from scipy.io import savemat
savemat(outdir + 'grid', {'lper':lper, 'lpar':lpar})
savemat(outdir + 'SF2_avg' , {
                           'SF2b' :SF2b_avg[:,:],
                           'SF2u' :SF2u_avg[:,:],
                         })
