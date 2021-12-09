# -*- coding: utf-8 -*-
from load import *
from fft import *
from plots import *

print('\nplotting kspectrum\n')
outdir = './fig_kspectrum/'

# If shear is on, add a vertical line indicating the fastest-MRI mode
if shear_flg == 1:
  b00 = np.max(b0[0])
  k_mri = np.sqrt(15.)/4./b00

def add_vertical_line(xs, ys, ls, legends):
  xs.append([k_mri, k_mri])
  ys.append([np.min([u2_bin[final_idx, 1:-1].min(), b2_bin[final_idx, 1:-1].min()]), 
             np.max([u2_bin[final_idx, 1:-1].max(), b2_bin[final_idx, 1:-1].max()])])
  ls.append('k:')
  legends.append(r'$k = k_\mr{MRI} = (\sqrt{15}/4)\varpi_0/v_\rmA$')

  return xs, ys, ls, legends

#--------------------------------------------------------#
#                      plot 1D spectra                   #
#--------------------------------------------------------#
ys = [ 
			 u2_bin[final_idx, 1:-1], 
			 b2_bin[final_idx, 1:-1], 
			 kpbin[1:-1]**(-5./3.)/kpbin[1]**(-5./3.)*b2_bin[final_idx,1:-1][0]
		 ]
xs = [ 
			kpbin[1:-1], 
			kpbin[1:-1], 
			kpbin[1:-1]  
		 ]
ls = [ 
				'', 
				'', 
				'k--', 
		 ]
legends = [ 
						r'$E_{u_\+}$', 
						r'$E_{\delta B_\+}$',
						r'-5/3',
					]

if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra.pdf')


# Elsasser fields
ys = [ 
			 zp2_bin[final_idx, 1:-1], 
			 zm2_bin[final_idx, 1:-1], 
			 kpbin[1:-1]**(-5./3.)/kpbin[1]**(-5./3.)*zp2_bin[final_idx,1:-1][0]
		 ]
xs = [ 
			kpbin[1:-1], 
			kpbin[1:-1], 
			kpbin[1:-1]  
		 ]
ls = [ 
				'', 
				'', 
				'k--', 
		 ]
legends = [ 
						r'$E_{Z^+}$', 
						r'$E_{Z^-}$',
						r'-5/3',
					]

if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra_ELS.pdf')


# Thermodynamic fields
ys = [ 
			 rho_bin[final_idx, 1:-1], 
			 sgm_bin[final_idx, 1:-1], 
			 kpbin[1:-1]**(-5./3.)/kpbin[1]**(-5./3.)*rho_bin[final_idx,1:-1][0]
		 ]
xs = [ 
			kpbin[1:-1], 
			kpbin[1:-1], 
			kpbin[1:-1]  
		 ]
ls = [ 
				'', 
				'', 
				'k--', 
		 ]
legends = [ 
						r'$E_{\rho}$', 
						r'$E_{\sigma}$',
						r'-5/3',
					]

if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra_therm.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ek.txt'  , np.column_stack((kpbin[:], 
																								 u2_bin [final_idx,:],
																								 b2_bin [final_idx,:],
                                                 rho_bin[final_idx,:],
                                                 sgm_bin[final_idx,:],
																								 )), fmt='%E')

#--------------------------------------------------------#
#                       plot movie                       #
#--------------------------------------------------------#
if ismovie:
  ys = []
  ys.append(u2_bin[:, 1:-1])
  ys.append(b2_bin[:, 1:-1])
  ys.append(np.tile(kpbin[1:-1], (nt, 1))**(-5./3.)/kpbin[1]*np.tile(u2_bin[:,1], (kpbin.size-2, 1)).T)
  ys = np.transpose(ys, (1, 0, 2))

  xs = [ 
        kpbin[1:-1],
        kpbin[1:-1],
        kpbin[1:-1]  
       ]
  ls = [ 
          '', 
          '', 
          'k--', 
       ]
  legends = [ 
              r'$E_{u}$', 
              r'$E_{b}$', 
              r'-5/3',
            ]
  movie_log1d_many(tt, xs, ys, legends=legends, ls=ls, legendloc='lower left', xlab='$k_\+ L_\+$', save=outdir+'k_spectra.gif')


  # Elsasser fields
  ys = []
  ys.append(zp2_bin[:, 1:-1])
  ys.append(zm2_bin[:, 1:-1])
  ys.append(np.tile(kpbin[1:-1], (nt, 1))**(-5./3.)/kpbin[1]*np.tile(zp2_bin[:,1], (kpbin.size-2, 1)).T)
  ys = np.transpose(ys, (1, 0, 2))

  xs = [ 
        kpbin[1:-1],
        kpbin[1:-1],
        kpbin[1:-1]  
       ]
  ls = [ 
          '', 
          '', 
          'k--', 
       ]
  legends = [ 
              r'$E_{Z^+}$', 
              r'$E_{Z^-}$', 
              r'-5/3',
            ]
  movie_log1d_many(tt, xs, ys, legends=legends, ls=ls, legendloc='lower left', xlab='$k_\+ L_\+$', save=outdir+'k_spectra_ELS.gif')


  # Thermodynamic fields
  ys = []
  ys.append(rho_bin[:, 1:-1])
  ys.append(sgm_bin[:, 1:-1])
  ys.append(np.tile(kpbin[1:-1], (nt, 1))**(-5./3.)/kpbin[1]*np.tile(rho_bin[:,1], (kpbin.size-2, 1)).T)
  ys = np.transpose(ys, (1, 0, 2))

  xs = [ 
        kpbin[1:-1],
        kpbin[1:-1]  
       ]
  ls = [ 
          '', 
          '', 
          'k--', 
       ]
  legends = [ 
              r'$E_{\rho}$', 
              r'$E_{\sigma}$', 
              r'-5/3',
            ]
  movie_log1d_many(tt, xs, ys, legends=legends, ls=ls, legendloc='lower left', xlab='$k_\+ L_\+$', save=outdir+'k_spectra_therm.gif')
