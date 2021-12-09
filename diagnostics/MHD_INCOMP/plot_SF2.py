# -*- coding: utf-8 -*-
from load import *
from fft import *
from plots import *

print('\nplotting SF2\n')
outdir = './fig_SF2/'

plot_SF2(SF2b[final_SF2_idx, :,:], SF2u[final_SF2_idx, :,:], lpar, lper, xlab=r'$\ell_\|$', ylab=r'$\ell_\+$', title=r'$t = %.2E$' % tt_SF2[final_SF2_idx], cmp=parula_map, save=outdir+'SF2.pdf')

ys = [ 
       SF2b[final_SF2_idx, 1:, 0 ], 
       SF2u[final_SF2_idx, 1:, 0 ], 
       SF2b[final_SF2_idx, 0 , 1:], 
       SF2u[final_SF2_idx, 0 , 1:], 
       lper[1:]**(2./3.)/lper[1]**(2./3.)*SF2b[final_SF2_idx, 1, 0],
       lpar[1:]**(3./3.)/lpar[1]**(3./3.)*SF2b[final_SF2_idx, 0, 1],
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
plot_log1d_many(xs, ys, xlab=r'$\ell_\|$ or $\ell_\+$', legends=legends, ls=ls, legendloc='upper left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'SF2_0.pdf')

# kpar(kper)
SF2b_0_lper = SF2b[final_SF2_idx, : , 0] # SF2b(0, lper)
SF2b_lpar_0 = SF2b[final_SF2_idx, 0 , :] # SF2b(lpar, 0)
SF2u_0_lper = SF2u[final_SF2_idx, : , 0] # SF2u(0, lper)
SF2u_lpar_0 = SF2u[final_SF2_idx, 0 , :] # SF2u(lpar, 0)

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
plot_log1d_many(xs, ys, xlab=r'$k_\+$', ylab=r'$k_\|$', legends=legends, ls=ls, legendloc='upper left', title=r'$t = %.2E $' % tt[final_idx], term=True, save=outdir+'kpar.pdf')

#------------------#
#    output data   #
#------------------#
np.savetxt(outdir + 'kpar.txt'  , np.column_stack((kper, kpar_b, kpar_u)), fmt='%E')

from scipy.io import savemat
savemat(outdir + 'grid', {'lper':lper, 'lpar':lpar})
savemat(outdir + 'SF2' , {
                           'SF2b' :SF2b[final_SF2_idx, :,:],
                           'SF2u' :SF2u[final_SF2_idx, :,:],
                         })
