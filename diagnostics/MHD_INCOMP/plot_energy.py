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
       u2dot_sum + b2dot_sum, 
       u2dissip_sum, 
       b2dissip_sum, 
       -p_ext_ene_sum, 
       -p_re_sum, 
       -p_ma_sum, 
			 u2dot_sum + b2dot_sum + u2dissip_sum + b2dissip_sum - p_ext_ene_sum - p_re_sum - p_ma_sum,
     ]
xs = [
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
        'k--', 
     ]
legends = [ 
       r'$\rmd W/\rmd t$', 
       r'$D_u$', 
       r'$D_B$', 
       r'$P_\mr{ext}$', 
       r'$P_\mr{Re}$', 
       r'$P_\mr{M}$', 
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
plot_1d(tt, (zp2_sum - zm2_sum)/(zp2_sum + zm2_sum), xlab='$'+tlab+'$', ylab=r'$H := \f{\int\rmd^3\bm{x}[(Z^+)^2 - (Z^-)^2]}{\int\rmd^3\bm{x}[(Z^+)^2 + (Z^-)^2]}$', ymin=-1.0, ymax=1.0, save=outdir+'helicity.pdf')
ys = [ 
       p_ext_ene_sum, p_ext_xhl_sum, 
     ]
ls = [ 
       '', '', 
     ]
xs = [
       tt, tt,
     ]
legends = [ 
       r'$\int \rmd^3\bm{x}(\bm{u}\cdot\bm{f}_u + \bm{b}\cdot\bm{f}_b)$', r'$\int \rmd^3\bm{x}(\bm{u}\cdot\bm{f}_b + \bm{b}\cdot\bm{f}_u)$', 
     ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', ymin=0.0, term=True, save=outdir + 'helicity_injection.pdf')

# plot rms
ys = [ 
       np.sqrt(ux_rms**2 + uy_rms**2 + uz_rms**2),
       np.sqrt((bx_rms - b0[0][0])**2 + (by_rms - b0[0][1])**2 + (bz_rms - b0[0][2])**2)
     ]
xs = [
       tt,
       tt,
     ]
ls = [ 
        '', 
        '', 
     ]
legends = [ 
       r'$|\bm{u}|_\mr{rms}$', 
       r'$|\delta \bm{B}|_\mr{rms}$', 
     ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'rms.pdf')

# plot stress
ys = [ 
       p_re_sum/zz.max()**2, 
       p_ma_sum/zz.max()**2, 
       (p_re_sum + p_ma_sum)/zz.max()**2, 
     ]
xs = [
       tt,
       tt,
       tt,
     ]
ls = [ 
        '', 
        '', 
        '', 
     ]
legends = [ 
       r'$\alpha_\mr{Re}$', 
       r'$\alpha_\mr{M}$', 
       r'$\alpha_\mr{Re} + \alpha_\mr{M}$', 
     ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'stress.pdf')


# ascii output
np.savetxt(outdir + 'energies.txt' , np.column_stack((tt, u2_sum, b2_sum, 
	                                                        u2dot_sum, b2dot_sum,
	                                                        u2dissip_sum, b2dissip_sum,
	                                                        p_ext_ene_sum, p_re_sum, p_ma_sum,
	                                                        zp2_sum, zm2_sum,
																													)), fmt='%E')
np.savetxt(outdir + 'balance.txt' , np.column_stack((tt, u2dot_sum + b2dot_sum + u2dissip_sum + b2dissip_sum - p_ext_ene_sum - p_re_sum - p_ma_sum
																													)), fmt='%E')
np.savetxt(outdir + 'energy_dot.txt' , np.column_stack((tt, u2dot_sum + b2dot_sum
																													)), fmt='%E')
np.savetxt(outdir + 'rms.txt' , np.column_stack((tt, ux_rms, uy_rms, uz_rms, 
                                                     np.sqrt((bx_rms - b0[0][0])**2),
                                                     np.sqrt((by_rms - b0[0][1])**2),
                                                     np.sqrt((bz_rms - b0[0][2])**2),
																													)), fmt='%E')

np.savetxt(outdir + 'stress.txt' , np.column_stack((tt, p_re_sum/zz.max()**2, p_ma_sum/zz.max()**2)), fmt='%E')


# calculate balance
balance = u2dot_sum + b2dot_sum + u2dissip_sum + b2dissip_sum - p_ext_ene_sum - p_re_sum - p_ma_sum

print ('|balance| > 1e1 at')
print (np.where(abs(balance) > 1e1)[0])

