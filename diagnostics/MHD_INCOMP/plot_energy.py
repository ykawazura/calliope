# -*- coding: utf-8 -*-
from load import *
from fft import *
from plots import *

print('\nplotting energy\n')
outdir = './fig_energy/'

# plot energy balance
ys = [ 
       u2dot_sum + b2dot_sum, 
       u2dissip_sum, 
       b2dissip_sum, 
       -p_u_sum, 
			 u2dot_sum + b2dot_sum + u2dissip_sum + b2dissip_sum - p_u_sum,
     ]
xs = [
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
        'k--', 
     ]
legends = [ 
       r'$\rmd W/\rmd t$', 
       r'$D_u$', 
       r'$D_B$', 
       r'$P$', 
			 r'balance', 
     ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'balance_all.pdf')


# plot energy change
ys = [ 
       u2_sum, b2_sum, 
     ]
ls = [ 
       '', '', 
     ]
xs = [
       tt, tt,
     ]
legends = [ 
       r'$W_u$', r'$W_B$', 
     ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'energy_all.pdf')


# plot helicity change
plot_1d(tt, (zp2_sum - zm2_sum)/(zp2_sum + zm2_sum), xlab='$'+tlab+'$', ylab=r'$\sigma := \f{\int\rmd^3\bm{x}[(Z^+)^2 - (Z^-)^2]}{\int\rmd^3\bm{x}[(Z^+)^2 + (Z^-)^2]}$', ymin=-1.0, ymax=1.0, save=outdir+'sigma.pdf')


# ascii output
np.savetxt(outdir + 'energies.txt' , np.column_stack((tt, u2_sum, b2_sum, 
	                                                        u2dot_sum, b2dot_sum,
	                                                        u2dissip_sum, b2dissip_sum,
	                                                        p_u_sum,
	                                                        zp2_sum, zm2_sum,
																													)), fmt='%E')
np.savetxt(outdir + 'balance.txt' , np.column_stack((tt, u2dot_sum + b2dot_sum + u2dissip_sum + b2dissip_sum - p_u_sum
																													)), fmt='%E')
np.savetxt(outdir + 'energy_dot.txt' , np.column_stack((tt, u2dot_sum + b2dot_sum
																													)), fmt='%E')


# calculate balance
balance = u2dot_sum + b2dot_sum + u2dissip_sum + b2dissip_sum - p_u_sum

print ('|balance| > 1e1 at')
print (np.where(abs(balance) > 1e1)[0])

