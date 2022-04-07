# -*- coding: utf-8 -*-
from load import *
from fft import *
from plots import *

print('\nplotting kspectrum\n')
outdir = './fig_kspectrum/'

if nlz == nkz:
  kp_end = np.argmin(np.abs(kpbin - kpbin.max()*2./3.))
else:
  kp_end = kpbin.size - 1

def add_vertical_line(xs, ys, ls, legends):
  xs.append([1./de, 1./de])
  ys.append([b2_bin[final_idx, 1:kp_end].min(), 
             b2_bin[final_idx, 1:kp_end].max()])
  ls.append('k:')
  legends.append(r'$kd_\rme = 1$')

  return xs, ys, ls, legends

#--------------------------------------------------------#
#                      plot 1D spectra                   #
#--------------------------------------------------------#
ys = [ 
       b2_bin[final_idx, 1:kp_end], 
       kpbin[1:kp_end]**(-7./3.)/kpbin[1]**(-7./3.)*b2_bin[final_idx,1:kp_end][0],
     ]
xs = [ 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
     ]
ls = [ 
        '', 
        'k--', 
     ]
legends = [ 
            r'$E$',
            r'-7/3',
          ]
if de > 0.0: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra.pdf')


# B by components
ys = [ 
       bx2_bin[final_idx, 1:kp_end], 
       by2_bin[final_idx, 1:kp_end], 
       bz2_bin[final_idx, 1:kp_end], 
       kpbin[1:kp_end]**(-7./3.)/kpbin[1]**(-7./3.)*b2_bin[final_idx,1:kp_end][0],
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
            r'$E_{B_x}$', 
            r'$E_{B_y}$', 
            r'$E_{B_z}$', 
            r'-7/3',
          ]

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra_b.pdf')

# kpar
ys = [ 
       kpar_b[final_kpar_idx, 1:kp_end], 
       kpbin[1:kp_end]**(2./3.)/kpbin[1]**(2./3.)*kpar_b[final_kpar_idx,1:kp_end][0],
       kpbin[1:kp_end]**(1./3.)/kpbin[1]**(1./3.)*kpar_b[final_kpar_idx,1:kp_end][0],
     ]
xs = [ 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
     ]
ls = [ 
        '', 
        'k--', 
        'k-.', 
     ]
legends = [ 
            r'$k_\|$', 
            r'2/3',
            r'1/3',
          ]

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt_kpar[final_kpar_idx], ylab='', term=True, save=outdir+'kpar.pdf')


# # delta b/b0
# ys = [ 
       # b1_ovr_b0[final_kpar_idx, 1:kp_end], 
     # ]
# xs = [ 
      # kpbin[1:kp_end], 
     # ]
# ls = [ 
        # '', 
     # ]
# legends = [ 
            # r'$\delta b/b_0$', 
          # ]

# plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt_kpar[final_kpar_idx], ylab='', term=True, save=outdir+'b1_ovr_b0.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ek.txt'  , np.column_stack((kpbin    [:kp_end], 
                                                  b2_bin  [final_idx,:kp_end],
                                                 bx2_bin  [final_idx,:kp_end],
                                                 by2_bin  [final_idx,:kp_end],
                                                 bz2_bin  [final_idx,:kp_end],
                                                 )), fmt='%E')
np.savetxt(outdir + 'kpar.txt'  , np.column_stack((kpbin  [:kp_end], 
                                                 kpar_b   [final_kpar_idx,:kp_end],
                                                 # b1_ovr_b0[final_kpar_idx,:kp_end],
                                                 )), fmt='%E')

#--------------------------------------------------------#
#                       plot movie                       #
#--------------------------------------------------------#
if ismovie:
  ys = []
  ys.append(b2_bin[:, 1:kp_end])
  ys.append(np.tile(kpbin[1:kp_end], (nt, 1))**(-7./3.)/kpbin[1]*np.tile(b2_bin[:,1], (kp_end-1, 1)).T)
  ys = np.transpose(ys, (1, 0, 2))

  xs = [ 
        kpbin[1:kp_end],
        kpbin[1:kp_end],
        kpbin[1:kp_end]  
       ]
  ls = [ 
          '', 
          'k--', 
          'k-.', 
       ]
  legends = [ 
              r'$E_{b}$', 
              r'-7/3',
            ]
  movie_log1d_many(tt, xs, ys, legends=legends, ls=ls, legendloc='lower left', xlab='$k_\+ L_\+$', save=outdir+'Ek.gif')

#--------------------------------------------------------#
#                       2D spectra                       #
#--------------------------------------------------------#
ncfile = netcdf.netcdf_file(input_dir+runname+'.out.2D.nc'+restart_num, 'r')   

tt_fld  = np.copy(ncfile.variables['tt' ][:]); tt_fld  = np.delete(tt_fld , ignored_points_fld, axis = 0)
kx_fld  = np.copy(ncfile.variables['kx' ][:])
ky_fld  = np.copy(ncfile.variables['ky' ][:])
kz_fld  = np.copy(ncfile.variables['kz' ][:])

nkx_mid = np.where(np.sign(kx_fld) == -1)[0][0]
nkz_mid = np.where(np.sign(kz_fld) == -1)[0][0]
kx_fld_ = np.concatenate([kx_fld[nkx_mid:], kx_fld[0:nkx_mid]])
kz_fld_ = np.concatenate([kz_fld[nkz_mid:], kz_fld[0:nkz_mid]])

b2_kxy = np.copy(ncfile.variables['b2_kxy'][:]);  b2_kxy = np.delete(b2_kxy, ignored_points_fld, axis = 0); b2_kxy = np.concatenate([b2_kxy[:, :, nkx_mid:], b2_kxy[:, :, 0:nkx_mid]], axis=2)
b2_kyz = np.copy(ncfile.variables['b2_kyz'][:]);  b2_kyz = np.delete(b2_kyz, ignored_points_fld, axis = 0); b2_kyz = np.concatenate([b2_kyz[:, nkz_mid:, :], b2_kyz[:, 0:nkz_mid, :]], axis=1)
b2_kxz = np.copy(ncfile.variables['b2_kxz'][:]);  b2_kxz = np.delete(b2_kxz, ignored_points_fld, axis = 0); b2_kxz = np.concatenate([b2_kxz[:, nkz_mid:, :], b2_kxz[:, 0:nkz_mid, :]], axis=1); b2_kxz = np.concatenate([b2_kxz[:, :, nkx_mid:], b2_kxz[:, :, 0:nkx_mid]], axis=2)

plot_2d(np.log10(b2_kxy[final_idx,:,:]), kx_fld_, ky_fld , xlab=r'$k_x$', ylab=r'$k_y$', title=r'$E_{b} (t = %.2E) $' % tt[final_idx], save=outdir+'b2_kxy.pdf')
plot_2d(np.log10(b2_kyz[final_idx,:,:]), ky_fld , kz_fld_, xlab=r'$k_y$', ylab=r'$k_z$', title=r'$E_{b} (t = %.2E) $' % tt[final_idx], save=outdir+'b2_kyz.pdf')
plot_2d(np.log10(b2_kxz[final_idx,:,:]), kx_fld_, kz_fld_, xlab=r'$k_x$', ylab=r'$k_z$', title=r'$E_{b} (t = %.2E) $' % tt[final_idx], save=outdir+'b2_kxz.pdf')
