# -*- coding: utf-8 -*-
from load import *
from fft import *
import sys
sys.path.append('../')
from plots import *

print('\nplotting energy\n')
outdir = './fig_energy/'

# plot energy balance
ys = [ 
       upe2dot_sum + bpe2dot_sum + upa2dot_sum + bpa2dot_sum, 
       upe2dissip_sum + bpe2dissip_sum, 
       upa2dissip_sum + bpa2dissip_sum, 
       -p_aw_sum, 
       -p_compr_sum, 
			 upe2dot_sum + bpe2dot_sum + upa2dot_sum + bpa2dot_sum + upe2dissip_sum + bpe2dissip_sum + upa2dissip_sum + bpa2dissip_sum - p_aw_sum - p_compr_sum,
     ]
xs = [
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
        'k--', 
     ]
legends = [ 
       r'$\rmd W/\rmd t$', 
       r'$D_\mr{AW}$', 
       r'$D_\mr{compr}$', 
       r'$-P_\mr{AW}$', 
       r'$-P_\mr{compr}$', 
			 r'balance', 
     ]
plot_1d_many(xs, ys, xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'balance_all.pdf')


# plot energy change
ys = [ 
       upe2_sum, 
       bpe2_sum, 
       upa2_sum, 
       bpa2_sum, 
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
       tt,
       tt,
     ]
legends = [ 
       r'$W_{u_\+}$', 
       r'$W_{\delta B_\+}$', 
       r'$W_{u_\|}$', 
       r'$W_{\delta B_\|}$', 
     ]
plot_1d_many(xs, ys, xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'energy_all.pdf')


# ascii output
np.savetxt(outdir + 'energies.txt' , np.column_stack((tt, upe2_sum, bpe2_sum, upa2_sum, bpa2_sum, 
	                                                        upe2dot_sum, bpe2dot_sum, upa2dot_sum, bpa2dot_sum,
	                                                        upe2dissip_sum, bpe2dissip_sum, upa2dissip_sum, bpa2dissip_sum,
	                                                        p_aw_sum, p_compr_sum, 
                                                          zpep2_sum, zpem2_sum, zpap2_sum, zpam2_sum
																													)), fmt='%E')
np.savetxt(outdir + 'balance.txt' , np.column_stack((tt, upe2dot_sum + bpe2dot_sum + upa2dot_sum + bpa2dot_sum + upe2dissip_sum + bpe2dissip_sum + upa2dissip_sum + bpa2dissip_sum - p_aw_sum - p_compr_sum
																													)), fmt='%E')
np.savetxt(outdir + 'energy_dot.txt' , np.column_stack((tt, upe2dot_sum + bpe2dot_sum + upa2dot_sum + bpa2dot_sum
																													)), fmt='%E')



# calculate balance
balance = upe2dot_sum + bpe2dot_sum + upa2dot_sum + bpa2dot_sum + upe2dissip_sum + bpe2dissip_sum + upa2dissip_sum + bpa2dissip_sum - p_aw_sum - p_compr_sum

print ('|balance| > 1e1 at')
print (np.where(abs(balance) > 1e1)[0])
