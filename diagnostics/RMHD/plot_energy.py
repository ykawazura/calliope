# -*- coding: utf-8 -*-
from load import *
from fft import *
from plots import *

print('\nplotting energy\n')
outdir = './fig_energy/'

# plot energy balance
ys = [ 
       upe2dot_sum + bpe2dot_sum , 
       upe2dissip_sum + bpe2dissip_sum, 
       -p_omg_sum, 
       -p_psi_sum, 
			 upe2dot_sum + bpe2dot_sum + upe2dissip_sum + bpe2dissip_sum - p_omg_sum - p_psi_sum,
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
       r'$D_\mr{AW}$', 
       r'$-P_{\nbl_\+^2\Phi}$', 
       r'$-P_{\Psi}$', 
			 r'balance', 
     ]
plot_1d_many(xs, ys, xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'balance_all.pdf')


# plot energy change
ys = [ 
       upe2_sum, 
       bpe2_sum, 
     ]
ls = [ 
        '', 
        '', 
     ]
xs = [
       tt,
       tt,
     ]
legends = [ 
       r'$W_{u_\+}$', 
       r'$W_{\delta B_\+}$', 
     ]
plot_1d_many(xs, ys, xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'energy_all.pdf')


# ascii output
np.savetxt(outdir + 'energies.txt' , np.column_stack((tt, upe2_sum, bpe2_sum, 
	                                                        upe2dot_sum, bpe2dot_sum,
	                                                        upe2dissip_sum, bpe2dissip_sum,
	                                                        p_omg_sum, p_psi_sum
																													)), fmt='%E')
np.savetxt(outdir + 'balance.txt' , np.column_stack((tt, upe2dot_sum + bpe2dot_sum + upe2dissip_sum + bpe2dissip_sum - p_omg_sum - p_psi_sum
																													)), fmt='%E')
np.savetxt(outdir + 'energy_dot.txt' , np.column_stack((tt, upe2dot_sum + bpe2dot_sum
																													)), fmt='%E')



# calculate balance
balance = upe2dot_sum + bpe2dot_sum + upe2dissip_sum + bpe2dissip_sum - p_omg_sum - p_psi_sum

print ('|balance| > 1e1 at')
print (np.where(abs(balance) > 1e1)[0])
