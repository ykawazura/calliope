# -*- coding: utf-8 -*-
from load import *
from fft import *
from plots import *

print('\nplotting kspectrum\n')
outdir = './fig_kspectrum/'

upe2_bin    = sum_negative_kz2d(upe2_bin)
bpe2_bin    = sum_negative_kz2d(bpe2_bin)
aaa2_bin    = sum_negative_kz2d(aaa2_bin)
ux2_bin     = sum_negative_kz2d(ux2_bin)
uy2_bin     = sum_negative_kz2d(uy2_bin)
bx2_bin     = sum_negative_kz2d(bx2_bin)
by2_bin     = sum_negative_kz2d(by2_bin)
zpep2_bin   = sum_negative_kz2d(zpep2_bin)
zpem2_bin   = sum_negative_kz2d(zpem2_bin)
p_bin       = sum_negative_kz2d(p_bin)

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
             np.sum(aaa2_bin[final_idx, :, 1:kp_end], axis=0), 
			 kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*np.sum(bpe2_bin[final_idx,:,1:kp_end], axis=0)[0]
		 ]
xs = [ 
			kpbin[1:kp_end], 
			kpbin[1:kp_end], 
			kpbin[1:kp_end], 
			kpbin[1:kp_end], 
		 ]
ls = [ 
				'', 
				'', 
				'', 
				'k--', 
		 ]
legends = [ 
            r'$E_{u_\+}$', 
            r'$E_{\delta B_\+}$',
            r'$E_{a}$', 
            r'-5/3',
					]
plot_log1d_many(xs, ys, xlab='$k_\+ L_\+$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kprp_spectra.pdf')

# kprp spectrum by components
ys = [ 
       np.sum(upe2_bin[final_idx, :, 1:kp_end], axis=0), 
       np.sum(ux2_bin [final_idx, :, 1:kp_end], axis=0), 
       np.sum(uy2_bin [final_idx, :, 1:kp_end], axis=0), 
       np.sum(bpe2_bin[final_idx, :, 1:kp_end], axis=0), 
       np.sum(bx2_bin [final_idx, :, 1:kp_end], axis=0), 
       np.sum(by2_bin [final_idx, :, 1:kp_end], axis=0), 
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*np.sum(bpe2_bin[final_idx,:,1:kp_end], axis=0)[0]
     ]
xs = [ 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end]  
     ]
ls = [ 
        'r-', 
        'r--', 
        'r:', 
        'b-', 
        'b--', 
        'b:', 
        'k--', 
     ]
legends = [ 
            r'$E_{u_\+}$', 
            r'$E_{u_x}$', 
            r'$E_{u_y}$', 
            r'$E_{\delta B_\+}$',
            r'$E_{\delta B_x}$',
            r'$E_{\delta B_y}$',
            r'-5/3',
          ]
plot_log1d_many(xs, ys, xlab='$k_\+ L_\+$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kprp_spectra_components.pdf')

# Elsasser fields
ys = [ 
       np.sum(zpep2_bin   [final_idx, :, 1:kp_end], axis=0), 
       np.sum(zpem2_bin   [final_idx, :, 1:kp_end], axis=0), 
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*np.sum(zpep2_bin[final_idx,:,1:kp_end], axis=0)[0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*np.sum(zpep2_bin[final_idx,:,1:kp_end], axis=0)[0]
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
        'k--', 
     ]
legends = [ 
            r'$E_{Z^+_\+}$', 
            r'$E_{Z^-_\+}$', 
            r'-5/3',
            r'-3/2',
          ]
plot_log1d_many(xs, ys, xlab='$k_\+ L_\+$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kprp_spectra_ELS.pdf')

# kz spectrum
if not is2D:
  ys = [ 
         np.sum(upe2_bin[final_idx, 1:kz_end, :kp_end], axis=1), 
         np.sum(bpe2_bin[final_idx, 1:kz_end, :kp_end], axis=1), 
         np.sum(aaa2_bin[final_idx, 1:kz_end, :kp_end], axis=1), 
       ]
  xs = [ 
          kz[1:kz_end], 
          kz[1:kz_end], 
          kz[1:kz_end], 
       ]
  ls = [ 
          '', 
          '', 
          '', 
       ]
  legends = [ 
              r'$E_{u_\+}$', 
              r'$E_{\delta B_\+}$',
              r'$E_{a}$', 
            ]
  plot_log1d_many(xs, ys, xlab='$'+kzlab+'$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kz_spectra.pdf')

  # Elsasser fields
  ys = [ 
         np.sum(zpep2_bin[final_idx, 1:kz_end, :kp_end], axis=1), 
         np.sum(zpem2_bin[final_idx, 1:kz_end, :kp_end], axis=1), 
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
              r'$E_{Z^+_\+}$', 
              r'$E_{Z^-_\+}$', 
            ]
  plot_log1d_many(xs, ys, xlab='$'+kzlab+'$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kz_spectra_ELS.pdf')

#--------------------------------------------------------#
#                      plot 2D spectra                   #
#--------------------------------------------------------#
if not is2D:
  plot_log2d(upe2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$k_\+ L_\+$', ylab='$'+kzlab+'$', 
      title=r'$E_{u_{\+}}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'upe2.pdf')
  plot_log2d(bpe2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$k_\+ L_\+$', ylab='$'+kzlab+'$', 
      title=r'$E_{\delta B_\+}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'bpe2.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ekprp.txt'  , np.column_stack((kpbin[:kp_end], 
                                                       np.sum(upe2_bin            [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(bpe2_bin            [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(aaa2_bin            [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ux2_bin             [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(uy2_bin             [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(bx2_bin             [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(by2_bin             [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(zpep2_bin           [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(zpem2_bin           [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(p_bin               [final_idx,:kz_end,:kp_end], axis=0),
                                                     )), fmt='%E')
if not is2D:
  np.savetxt(outdir + 'Ekz.txt'  , np.column_stack((kz[:kz_end], 
                                                         np.sum(upe2_bin            [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(bpe2_bin            [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(aaa2_bin            [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ux2_bin             [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(uy2_bin             [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(bx2_bin             [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(by2_bin             [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(zpep2_bin           [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(zpem2_bin           [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(p_bin               [final_idx,:kz_end,:kp_end], axis=1),
                                                       )), fmt='%E')

del upe2_bin
del bpe2_bin
del aaa2_bin
