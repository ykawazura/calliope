# -*- coding: utf-8 -*-
from load import *
from fft import *
import sys
sys.path.append('../')
from plots import *

print('\nplotting kspectrum\n')
outdir = './fig_kspectrum/'

upe2_bin    = sum_negative_kz2d(upe2_bin)
bpe2_bin    = sum_negative_kz2d(bpe2_bin)
zppe2_bin   = sum_negative_kz2d(zppe2_bin)
zmpe2_bin   = sum_negative_kz2d(zmpe2_bin)
hp_bin      = sum_negative_kz2d(hp_bin   )
hm_bin      = sum_negative_kz2d(hm_bin   )

if nlz == nkz:
  kp_end = np.argmin(np.abs(kpbin - kpbin.max()*2./3.))
  if not is2D:
    kz_end = np.argmin(np.abs(kz[1:int(nkz/2)] - kz[1:int(nkz/2)].max()*2./3.))
else:
  kp_end = kpbin.size - 1
  kz_end = int(nkz/2)

#--------------------------------------------------------#
#                      plot 1D spectra                   #
#--------------------------------------------------------#
# kprp spectrum
ys = [ 
			 np.sum(upe2_bin[final_idx, :, 1:kp_end], axis=0), 
			 np.sum(bpe2_bin[final_idx, :, 1:kp_end], axis=0), 
			 kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*np.sum(bpe2_bin[final_idx,:,1:kp_end], axis=0)[0],
			 kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*np.sum(bpe2_bin[final_idx,:,1:kp_end], axis=0)[0]
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
				'k--', 
				'k-.', 
		 ]
legends = [ 
						r'$E_{u_\+}$', 
						r'$E_{\delta B_\+}$',
						r'-5/3',
						r'-3/2',
					]
plot_log1d_many(xs, ys, xlab='$k_\+ L_\+$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kprp_spectra.pdf')

## Elsasser fields
ys = [ 
			 np.sum(zppe2_bin[final_idx, :, 1:kp_end], axis=0), 
			 np.sum(zmpe2_bin[final_idx, :, 1:kp_end], axis=0), 
			 kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*np.sum(bpe2_bin[final_idx,:,1:kp_end], axis=0)[0],
			 kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*np.sum(bpe2_bin[final_idx,:,1:kp_end], axis=0)[0]
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
				'k--', 
				'k-.', 
		 ]
legends = [ 
            r'$E_{Z_\+^+}$', 
            r'$E_{Z_\+^-}$',
						r'-5/3',
						r'-3/2',
					]
plot_log1d_many(xs, ys, xlab='$k_\+ L_\+$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kprp_spectra_ELS.pdf')

# kz spectrum
if not is2D:
  ys = [ 
         np.sum(upe2_bin[final_idx, 1:kz_end, :kp_end], axis=1), 
         np.sum(bpe2_bin[final_idx, 1:kz_end, :kp_end], axis=1), 
       ]
  xs = [ 
          kz[1:kz_end], 
          kz[1:kz_end], 
       ]
  ls = [ 
          '', 
          '', 
       ]
  legends = [ 
              r'$E_{u_\+}$', 
              r'$E_{\delta B_\+}$',
            ]
  plot_log1d_many(xs, ys, xlab=kzlab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kz_spectra.pdf')

## Elsasser fields
if not is2D:
  ys = [ 
         np.sum(zppe2_bin[final_idx, 1:kz_end, :kp_end], axis=1), 
         np.sum(zmpe2_bin[final_idx, 1:kz_end, :kp_end], axis=1), 
       ]
  xs = [ 
          kz[1:kz_end], 
          kz[1:kz_end], 
       ]
  ls = [ 
          '', 
          '', 
       ]
  legends = [ 
              r'$E_{Z_\+^+}$', 
              r'$E_{Z_\+^-}$',
            ]
  plot_log1d_many(xs, ys, xlab=kzlab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kz_spectra_ELS.pdf')

#--------------------------------------------------------#
#                      plot 2D spectra                   #
#--------------------------------------------------------#
if not is2D:
  plot_log2d(upe2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$k_\+ L_\+$', ylab=kzlab, 
      title=r'$E_{u_{\+}}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'upe2.pdf')
  plot_log2d(bpe2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$k_\+ L_\+$', ylab=kzlab, 
      title=r'$E_{\delta B_\+}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'bpe2.pdf')
  plot_log2d(zppe2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$k_\+ L_\+$', ylab=kzlab, 
      title=r'$E_{Z_{\+}^+}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'zppe2.pdf')
  plot_log2d(zmpe2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$k_\+ L_\+$', ylab=kzlab, 
      title=r'$E_{Z_{\+}^-}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'zmpe2.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ekprp.txt'  , np.column_stack((kpbin[:kp_end], 
																											 np.sum(upe2_bin [final_idx,:kz_end,:kp_end], axis=0),
																											 np.sum(bpe2_bin [final_idx,:kz_end,:kp_end], axis=0),
																											 np.sum(zppe2_bin[final_idx,:kz_end,:kp_end], axis=0),
																											 np.sum(zmpe2_bin[final_idx,:kz_end,:kp_end], axis=0),
																										 )), fmt='%E')
if not is2D:
  np.savetxt(outdir + 'Ekz.txt'  , np.column_stack((kz[:kz_end], 
                                                         np.sum(upe2_bin [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(bpe2_bin [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(zppe2_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(zmpe2_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                       )), fmt='%E')

del upe2_bin
del bpe2_bin


#---------------------------------------#
#   time-correlations                   #
#   10.1103/PhysRevResearch.2.023357    #
#---------------------------------------#
L = xx.max() - xx.min()
urms = np.sqrt(ux_rms[-1]**2 + uy_rms[-1]**2)
tau = L/(2.0*np.pi)/urms

from scipy.optimize import curve_fit
def gaussian_func(x, sigma):
    return np.exp(- (x**2) / (2 * sigma**2))
initial_guess = [1.0]

line_num = 5
fig, axes = plt.subplots(1, 2, figsize=(16, 8))

for k in kpbin[0:-1:kpbin.size//line_num][1:]:
    i = np.argmin(abs(kpbin - k))

    # Fig 1
    ax = axes[0]
    x = (tt - tt[0])/tau
    y = np.sum(hp_bin[:, :, i], axis=1)/np.sum(zppe2_bin[0, :, i], axis=0)

    ax.plot(x, y, lw=3, label=r'$k_\+ = %2d$' % k)

    # fitting with Gaussian
    popt, _ = curve_fit(gaussian_func, x, y, p0=initial_guess)
    ax.plot(x, gaussian_func(x, popt[0]), 'k--', lw=2)

    ax.set_ylabel(r'$\Gamma^+(k_\+, \, \tau)$')


    # Fig 2
    ax = axes[1]
    x = (tt - tt[0])/tau
    y = np.sum(hp_bin[:, :, i], axis=1)/np.sum(zppe2_bin[0, :, i], axis=0)

    ax.plot(x, y, lw=3, label=r'$k_\+ = %2d$' % k)

    # fitting with Gaussian
    popt, _ = curve_fit(gaussian_func, x, y, p0=initial_guess)
    ax.plot(x, gaussian_func(x, popt[0]), 'k--', lw=2)

    ax.set_ylabel(r'$\Gamma^-(k_\+, \, \tau)$')

for ax in axes:
    ax.set_xlabel(r'$\tau$')
    ax.axhline(y=0.0, color='k', linestyle='-', lw=1)
    leg = ax.legend(frameon=False, loc='upper right')


plt.savefig(outdir+'time-correlations.pdf')
#---------------------------------------------
