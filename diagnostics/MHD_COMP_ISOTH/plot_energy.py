# -*- coding: utf-8 -*-
from load import *
from fft import *
from plots import *

print('\nplotting energy\n')
outdir = './fig_energy/'

# plot energy balance
ys = [ 
       wkin_dot_sum, 
       wmag_dot_sum, 
       wrho_dot_sum, 
       wkin_dissip_sum, 
       wmag_dissip_sum, 
       wrho_dissip_sum, 
       -p_ext_sum, 
       -p_re_sum, 
       -p_ma_sum, 
			 wkin_dot_sum + wmag_dot_sum + wrho_dot_sum + wkin_dissip_sum + wmag_dissip_sum + wrho_dissip_sum - p_ext_sum - p_re_sum - p_ma_sum,
     ]
xs = [
       tt,
       tt,
       tt,
       tt,
       tt,
       tt,
       tt,
       tt,
       tt,
       tt,
     ]
ls = [ 
        '', 
        '', 
        '', 
        '', 
        '', 
        '', 
        '', 
        '', 
        '', 
        'k--', 
     ]
legends = [ 
       r'$\rmd W_\mr{kin}/\rmd t$', 
       r'$\rmd W_\mr{mag}/\rmd t$', 
       r'$\rmd W_\rho/\rmd t$', 
       r'$D_\mr{kin}$', 
       r'$D_\mr{mag}$', 
       r'$D_\rho$', 
       r'$P_\mr{ext}$', 
       r'$P_\mr{Re}$', 
       r'$P_\mr{M}$', 
			 r'balance', 
     ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'balance_all.pdf')


# plot energy change
ys = [ 
       wkin_sum, 
       wmag_sum, 
       wrho_sum, 
     ]
ls = [ 
       '', 
       '', 
       '', 
     ]
xs = [
       tt,
       tt,
       tt,
     ]
legends = [ 
       r'$W_\mr{kin}$', 
       r'$W_\mr{mag}$', 
       r'$W_\rho$', 
     ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'energy_all.pdf')


# plot helicity change
plot_1d(tt, (zp2_sum - zm2_sum)/(zp2_sum + zm2_sum), xlab='$'+tlab+'$', ylab=r'$\sigma := \f{\int\rmd^3\bm{x}[(Z^+)^2 - (Z^-)^2]}{\int\rmd^3\bm{x}[(Z^+)^2 + (Z^-)^2]}$', ymin=-1.0, ymax=1.0, save=outdir+'sigma.pdf')


# plot mach number & beta
plot_1d(tt, smach_rms, xlab=tlab, ylab=r'$M_\mr{S,rms}$'  , save=outdir+'smach_rms.pdf')
plot_1d(tt, amach_rms, xlab=tlab, ylab=r'$M_\mr{A,rms}$'  , save=outdir+'amach_rms.pdf')
plot_1d(tt, beta_rms , xlab=tlab, ylab=r'$\beta_\mr{rms}$', save=outdir+'beta_rms.pdf' )


# ascii output
np.savetxt(outdir + 'energies.txt' , np.column_stack((tt, wkin_sum, wmag_sum, wrho_sum, 
	                                                        wkin_dot_sum, wmag_dot_sum, wrho_dot_sum,
	                                                        wkin_dissip_sum, wmag_dissip_sum, wrho_dissip_sum,
	                                                        p_ext_sum, p_re_sum, p_ma_sum
																													)), fmt='%E')
np.savetxt(outdir + 'balance.txt' , np.column_stack((tt, wkin_dot_sum + wmag_dot_sum + wrho_dot_sum + wkin_dissip_sum + wmag_dissip_sum + wrho_dissip_sum - p_ext_sum - p_re_sum - p_ma_sum
																													)), fmt='%E')
np.savetxt(outdir + 'energy_dot.txt' , np.column_stack((tt, wkin_dot_sum + wmag_dot_sum + wrho_dot_sum
																													)), fmt='%E')


# calculate balance
balance = wkin_dot_sum + wmag_dot_sum + wrho_dot_sum + wkin_dissip_sum + wmag_dissip_sum + wrho_dissip_sum - p_ext_sum - p_re_sum - p_ma_sum

print ('|balance| > 1e1 at')
print (np.where(abs(balance) > 1e1)[0])

