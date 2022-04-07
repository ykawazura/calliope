# -*- coding: utf-8 -*-
from load import *
from fft import *
from plots import *

print('\nplotting energy\n')
outdir = './fig_energy/'

# plot energy balance
ys = [ 
       wmag_dot_sum, 
       wmag_dissip_sum, 
       -p_ext_sum, 
			 wmag_dot_sum + wmag_dissip_sum - p_ext_sum,
     ]
xs = [
       tt,
       tt,
       tt,
       tt,
     ]
ls = [ 
        '', 
        '', 
        '', 
        'k--', 
     ]
legends = [ 
       r'$\rmd W/\rmd t$', 
       r'$D_B$', 
       r'$P_\mr{ext}$', 
			 r'balance', 
     ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'balance_all.pdf')


# plot energy change
ys = [ 
       wmag_sum, 
     ]
ls = [ 
       '', 
     ]
xs = [
       tt,
     ]
legends = [ 
       r'$W$', 
     ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'energy_all.pdf')



# ascii output
np.savetxt(outdir + 'energies.txt' , np.column_stack((tt, wmag_sum, 
	                                                        wmag_dot_sum,
	                                                        wmag_dissip_sum,
	                                                        p_ext_sum,
																													)), fmt='%E')
np.savetxt(outdir + 'balance.txt' , np.column_stack((tt, wmag_dot_sum + wmag_dissip_sum - p_ext_sum
																													)), fmt='%E')
np.savetxt(outdir + 'energy_dot.txt' , np.column_stack((tt, wmag_dot_sum
																													)), fmt='%E')


# calculate balance
balance = wmag_dot_sum + wmag_dissip_sum - p_ext_sum

print ('|balance| > 1e1 at')
print (np.where(abs(balance) > 1e1)[0])

