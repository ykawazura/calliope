# -*- coding: utf-8 -*-
from load import *
from fft import *
from plots import *
from scipy.integrate import simps, trapz
from scipy import interpolate

avg_start = 0
avg_end   = -1

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

Waw        = upe2_sum + bpe2_sum
Wtot       = Waw
Waw_dot    = upe2dot_sum + bpe2dot_sum
Wtot_dot   = Waw_dot
Daw        = upe2dissip_sum + bpe2dissip_sum
Dtot       = Daw
Paw        = p_omg_sum + p_psi_sum
Ptot       = Paw

Waw_avg        = time_average(tt[avg_start:avg_end], Waw        [avg_start:avg_end], axis=0)
Wtot_avg       = time_average(tt[avg_start:avg_end], Wtot       [avg_start:avg_end], axis=0)
Waw_dot_avg    = time_average(tt[avg_start:avg_end], Waw_dot    [avg_start:avg_end], axis=0)
Wtot_dot_avg   = time_average(tt[avg_start:avg_end], Wtot_dot   [avg_start:avg_end], axis=0)
Daw_avg        = time_average(tt[avg_start:avg_end], Daw        [avg_start:avg_end], axis=0)
Dtot_avg       = time_average(tt[avg_start:avg_end], Dtot       [avg_start:avg_end], axis=0)
Paw_avg        = time_average(tt[avg_start:avg_end], Paw        [avg_start:avg_end], axis=0)
Ptot_avg       = time_average(tt[avg_start:avg_end], Ptot       [avg_start:avg_end], axis=0)
upe2_sum_avg   = time_average(tt[avg_start:avg_end], upe2_sum   [avg_start:avg_end], axis=0)
bpe2_sum_avg   = time_average(tt[avg_start:avg_end], bpe2_sum   [avg_start:avg_end], axis=0)

Waw_err        = std_dev(tt[avg_start:avg_end], Waw        [avg_start:avg_end])
Wtot_err       = std_dev(tt[avg_start:avg_end], Wtot       [avg_start:avg_end])
Waw_dot_err    = std_dev(tt[avg_start:avg_end], Waw_dot    [avg_start:avg_end])
Wtot_dot_err   = std_dev(tt[avg_start:avg_end], Wtot_dot   [avg_start:avg_end])
Daw_err        = std_dev(tt[avg_start:avg_end], Daw        [avg_start:avg_end])
Dtot_err       = std_dev(tt[avg_start:avg_end], Dtot       [avg_start:avg_end])
Paw_err        = std_dev(tt[avg_start:avg_end], Paw        [avg_start:avg_end])
Ptot_err       = std_dev(tt[avg_start:avg_end], Ptot       [avg_start:avg_end])
upe2_sum_err   = std_dev(tt[avg_start:avg_end], upe2_sum   [avg_start:avg_end])
bpe2_sum_err   = std_dev(tt[avg_start:avg_end], bpe2_sum   [avg_start:avg_end])

s =     'average over t     = [%.3E' % tt[avg_start] + ', %.3E' % tt[avg_end] + ']' + '\n'
s = s + 'average over index = [' + str(avg_start) + ', ' + str(avg_end) + ']' + '\n'
s = s + ' error of ratio is calculated by (a + da)/(b + db) ~ a/b*[1 + sqrt( (da/a)^2 + (db/b)^2 - 2cov/(a*b))]' + '\n\n' 
s = s + '  Wtot             = %.3E \pm %.3E'  % (Wtot_avg      , Wtot_err      ) + '\n'
s = s + '    (Waw           = %.3E \pm %.3E)' % (Waw_avg       , Waw_err       ) + '\n'
s = s + '  Wtot_dot         = %.3E \pm %.3E'  % (Wtot_dot_avg  , Wtot_dot_err  ) + '\n'
s = s + '    (Waw_dot       = %.3E \pm %.3E)' % (Waw_dot_avg   , Waw_dot_err   ) + '\n'
s = s + '  Ptot             = %.3E \pm %.3E'  % (Ptot_avg      , Ptot_err      ) + '\n'
s = s + '    Paw            = %.3E \pm %.3E'  % (Paw_avg       , Paw_err       ) + '\n'
s = s + '  Dtot             = %.3E \pm %.3E'  % (Dtot_avg      , Dtot_err      ) + '\n'
s = s + '    Daw            = %.3E \pm %.3E'  % (Daw_avg       , Daw_err       ) + '\n'
print (s)

f = open('time_average.txt','w') 
f.write(s) 
f.close() 


ys = [ 
			 upe2dot_sum + bpe2dot_sum, 
			 upe2dissip_sum + bpe2dissip_sum, 
       -p_omg_sum, 
       -p_psi_sum, 
			 upe2dot_sum + bpe2dot_sum + upe2dissip_sum + bpe2dissip_sum - p_omg_sum - p_psi_sum,
			 np.full([tt[avg_start:avg_end].size], Waw_dot_avg),
			 np.full([tt[avg_start:avg_end].size], Daw_avg),
			 np.full([tt[avg_start:avg_end].size], -Paw_avg),
		 ]
xs = [
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
				'k--', 
				'', 
				'', 
				'', 
		 ]
legends = [ 
			 r'$\mathrm{d}W/\mathrm{d} t$', 
			 r'$D_\mathrm{AW}$', 
       r'$-P_{\nbl_\+^2\Phi}$', 
       r'$-P_{\Psi}$', 
			 r'balance', 
			 '', 
			 '', 
			 '', 
		 ]
plot_1d_many_average(xs, ys, tt[avg_start], tt[avg_end], xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'balance_all_avg.pdf')


ys = [ 
       upe2_sum, 
       bpe2_sum, 
			 np.full([tt[avg_start:avg_end].size], upe2_sum_avg),
			 np.full([tt[avg_start:avg_end].size], bpe2_sum_avg),
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
       r'$W_{u_\perp}$', 
       r'$W_{\delta B_\perp}$', 
			 '',
			 '',
     ]
plot_1d_many_average(xs, ys, tt[avg_start], tt[avg_end], xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'energy_all_avg.pdf')

##########################################################
#                    average kspectrum                   #
##########################################################
print('\nplotting kspectrum\n')
outdir = './fig_kspectrum/'

upe2_bin = sum_negative_kz2d(upe2_bin)
bpe2_bin = sum_negative_kz2d(bpe2_bin)

upe2_bin_avg      = time_average(tt[avg_start:avg_end], upe2_bin     [avg_start:avg_end], axis=0)
bpe2_bin_avg      = time_average(tt[avg_start:avg_end], bpe2_bin     [avg_start:avg_end], axis=0)

kp23i = np.argmin(np.abs(kpbin - kpbin.max()*2./3.))
kz23i = np.argmin(np.abs(kz[1:int(nkz/2)] - kz[1:int(nkz/2)].max()*2./3.))

# kperp spectrum
ys = [ 
       np.sum(upe2_bin_avg   [:, 1:kp23i], axis=0), 
       np.sum(bpe2_bin_avg   [:, 1:kp23i], axis=0), 
			 kpbin[1:kp23i]**(-5./3.)/kpbin[1]**(-5./3.)*np.sum(bpe2_bin_avg[:,1:kp23i], axis=0)[0]
     ]
xs = [ 
			kpbin[1:kp23i],
			kpbin[1:kp23i],
			kpbin[1:kp23i]  
     ]
ls = [ 
        '', 
				'', 
        'k--', 
     ]
legends = [ 
            r'$E_{u_{\perp}}$', 
            r'$E_{\delta B_{\perp}}$',
            r'-5/3',
          ]
plot_log1d_many(xs, ys, xlab='$k_\+ L_\+$', legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'kperp_spectra_avg.pdf')

# kz spectrum
ys = [ 
       np.sum(upe2_bin_avg   [1:kz23i, :], axis=1), 
       np.sum(bpe2_bin_avg   [1:kz23i, :], axis=1), 
     ]
xs = [ 
        kz[1:kz23i], 
        kz[1:kz23i], 
     ]
ls = [ 
        '', 
        '', 
     ]
legends = [ 
            r'$E_{u_{\perp}}$', 
            r'$E_{\delta B_{\perp}}$',
          ]
plot_log1d_many(xs, ys, xlab='$'+kzlab+'$', legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'kz_spectra_avg.pdf')

#--------------------------------------------------------#
#                      plot 2D spectra                   #
#--------------------------------------------------------#
plot_log2d(upe2_bin_avg[1:kz23i, 1:kp23i], kpbin[1:kp23i], kz[1:kz23i], xlab='$k_\+ L_\+$', ylab='$'+kzlab+'$', 
		title=r'$E_{u_{\perp}}$', save=outdir + 'upe2_avg.pdf')
plot_log2d(bpe2_bin_avg[1:kz23i, 1:kp23i], kpbin[1:kp23i], kz[1:kz23i], xlab='$k_\+ L_\+$', ylab='$'+kzlab+'$', 
		title=r'$E_{\delta B_{\perp}}$', save=outdir + 'bpe2_avg.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ekperp_avg.txt'  , np.column_stack((kpbin[:kp23i], 
                                                       np.sum(upe2_bin_avg   [:,:kp23i], axis=0),
                                                       np.sum(bpe2_bin_avg   [:,:kp23i], axis=0),
                                                     )), fmt='%E')
np.savetxt(outdir + 'Ekz_avg.txt'  , np.column_stack((kz[:kz23i], 
                                                       np.sum(upe2_bin_avg   [:kz23i,:], axis=1),
                                                       np.sum(bpe2_bin_avg   [:kz23i,:], axis=1),
                                                     )), fmt='%E')


del upe2_bin
del bpe2_bin
