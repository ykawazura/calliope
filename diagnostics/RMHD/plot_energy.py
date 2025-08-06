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
       upe2dot_sum + bpe2dot_sum , 
       upe2dissip_sum + bpe2dissip_sum, 
       -p_phi_sum, 
       -p_psi_sum, 
       -p_phi_sum - p_psi_sum, 
			 upe2dot_sum + bpe2dot_sum + upe2dissip_sum + bpe2dissip_sum - p_phi_sum - p_psi_sum,
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
       r'$-P_{\Phi}$', 
       r'$-P_{\Psi}$', 
       r'$-P_\mr{total}$', 
			 r'balance', 
     ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'balance_all.pdf')


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
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'energy_all.pdf')

# plot helicity change
plot_1d(tt, (zppe2_sum - zmpe2_sum)/(zppe2_sum + zmpe2_sum), xlab=tlab, ylab=r'$H := \f{\int\rmd^3\bm{x}[(Z_\+^+)^2 - (Z_\+^-)^2]}{\int\rmd^3\bm{x}[(Z_\+^+)^2 + (Z_\+^-)^2]}$', ymin=-1.0, ymax=1.0, save=outdir+'helicity.pdf')

ys = [ 
       p_phi_sum + p_psi_sum, p_xhl_sum, 
     ]
ls = [ 
       '', '', 
     ]
xs = [
       tt, tt,
     ]
legends = [ 
       r'$\int \rmd^3\bm{x}(-\zeta_\+^+\nbl^2 f_{\zeta_\+^+} - \zeta_\+^-\nbl^2 f_{\zeta_\+^-})$', r'$\int \rmd^3\bm{x}(-\zeta_\+^+\nbl^2 f_{\zeta_\+^+} + \zeta_\+^-\nbl^2 f_{\zeta_\+^-})$', 
     ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, ymin=0.0, save=outdir + 'helicity_injection.pdf')


# ascii output
np.savetxt(outdir + 'energies.txt' , np.column_stack((tt, upe2_sum, bpe2_sum, 
	                                                        upe2dot_sum, bpe2dot_sum,
	                                                        upe2dissip_sum, bpe2dissip_sum,
	                                                        p_phi_sum, p_psi_sum,
	                                                        zppe2_sum, zmpe2_sum,
																													)), fmt='%E')
np.savetxt(outdir + 'balance.txt' , np.column_stack((tt, upe2dot_sum + bpe2dot_sum + upe2dissip_sum + bpe2dissip_sum - p_phi_sum - p_psi_sum
																													)), fmt='%E')
np.savetxt(outdir + 'energy_dot.txt' , np.column_stack((tt, upe2dot_sum + bpe2dot_sum
																													)), fmt='%E')



# calculate balance
balance = upe2dot_sum + bpe2dot_sum + upe2dissip_sum + bpe2dissip_sum - p_phi_sum - p_psi_sum

print ('|balance| > 1e1 at')
print (np.where(abs(balance) > 1e1)[0])
