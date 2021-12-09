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

# If shear is on, add a vertical line indicating the fastest-MRI mode
if shear_flg == 1:
  b00 = np.max(b0[0])
  k_mri = np.sqrt(15.)/4./b00

def add_vertical_line(xs, ys, ls, legends):
  xs.append([k_mri, k_mri])
  ys.append([np.min([u2_bin[final_idx, 1:kp_end].min(), b2_bin[final_idx, 1:kp_end].min()]), 
             np.max([u2_bin[final_idx, 1:kp_end].max(), b2_bin[final_idx, 1:kp_end].max()])])
  ls.append('k:')
  legends.append(r'$k = k_\mr{MRI} = (\sqrt{15}/4)\varpi_0/v_\rmA$')

  return xs, ys, ls, legends

#--------------------------------------------------------#
#                      plot 1D spectra                   #
#--------------------------------------------------------#
ys = [ 
       u2_bin[final_idx, 1:kp_end], 
       b2_bin[final_idx, 1:kp_end], 
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*b2_bin[final_idx,1:kp_end][0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*b2_bin[final_idx,1:kp_end][0],
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
        'k--', 
        'k-.', 
     ]
legends = [ 
            r'$E_u$', 
            r'$E_B$',
            r'-5/3',
            r'-3/2',
          ]

if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra.pdf')


# u by components
ys = [ 
       ux2_bin[final_idx, 1:kp_end], 
       uy2_bin[final_idx, 1:kp_end], 
       uz2_bin[final_idx, 1:kp_end], 
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*u2_bin[final_idx,1:kp_end][0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*u2_bin[final_idx,1:kp_end][0],
     ]
xs = [ 
      kpbin[1:kp_end], 
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
        'k-.', 
     ]
legends = [ 
            r'$E_{u_x}$', 
            r'$E_{u_y}$', 
            r'$E_{u_z}$', 
            r'-5/3',
            r'-3/2',
          ]

if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra_u.pdf')


# B by components
ys = [ 
       bx2_bin[final_idx, 1:kp_end], 
       by2_bin[final_idx, 1:kp_end], 
       bz2_bin[final_idx, 1:kp_end], 
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*b2_bin[final_idx,1:kp_end][0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*b2_bin[final_idx,1:kp_end][0],
     ]
xs = [ 
      kpbin[1:kp_end], 
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
        'k-.', 
     ]
legends = [ 
            r'$E_{B_x}$', 
            r'$E_{B_y}$', 
            r'$E_{B_z}$', 
            r'-5/3',
            r'-3/2',
          ]

if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra_b.pdf')


# Elsasser fields
ys = [ 
       zp2_bin[final_idx, 1:kp_end], 
       zm2_bin[final_idx, 1:kp_end], 
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*zp2_bin[final_idx,1:kp_end][0],
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*zp2_bin[final_idx,1:kp_end][0],
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
        'k--', 
        'k-.', 
     ]
legends = [ 
            r'$E_{Z^+}$', 
            r'$E_{Z^-}$',
            r'-5/3',
            r'-3/2',
          ]

if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra_ELS.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ek.txt'  , np.column_stack((kpbin[:kp_end], 
                                                  u2_bin[final_idx,:kp_end],
                                                 ux2_bin[final_idx,:kp_end],
                                                 uy2_bin[final_idx,:kp_end],
                                                 uz2_bin[final_idx,:kp_end],
                                                  b2_bin[final_idx,:kp_end],
                                                 bx2_bin[final_idx,:kp_end],
                                                 by2_bin[final_idx,:kp_end],
                                                 bz2_bin[final_idx,:kp_end],
                                                 zp2_bin[final_idx,:kp_end],
                                                 zm2_bin[final_idx,:kp_end],
                                                 )), fmt='%E')

#--------------------------------------------------------#
#                       plot movie                       #
#--------------------------------------------------------#
if ismovie:
  ys = []
  ys.append(u2_bin[:, 1:kp_end])
  ys.append(b2_bin[:, 1:kp_end])
  ys.append(np.tile(kpbin[1:kp_end], (nt, 1))**(-5./3.)/kpbin[1]*np.tile(u2_bin[:,1], (kp_end-1, 1)).T)
  ys.append(np.tile(kpbin[1:kp_end], (nt, 1))**(-3./2.)/kpbin[1]*np.tile(u2_bin[:,1], (kp_end-1, 1)).T)
  ys = np.transpose(ys, (1, 0, 2))

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
              r'$E_{u}$', 
              r'$E_{b}$', 
              r'-5/3',
              r'-3/2',
            ]
  movie_log1d_many(tt, xs, ys, legends=legends, ls=ls, legendloc='lower left', xlab='$k_\+ L_\+$', save=outdir+'Ek.gif')

  # Elsasser fields
  ys = []
  ys.append(zp2_bin[:, 1:kp_end])
  ys.append(zm2_bin[:, 1:kp_end])
  ys.append(np.tile(kpbin[1:kp_end], (nt, 1))**(-5./3.)/kpbin[1]*np.tile(zp2_bin[:,1], (kp_end-1, 1)).T)
  ys.append(np.tile(kpbin[1:kp_end], (nt, 1))**(-3./2.)/kpbin[1]*np.tile(zp2_bin[:,1], (kp_end-1, 1)).T)
  ys = np.transpose(ys, (1, 0, 2))

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
              r'$E_{Z^+}$', 
              r'$E_{Z^-}$',
              r'-5/3',
              r'-3/2',
            ]
  movie_log1d_many(tt, xs, ys, legends=legends, ls=ls, legendloc='lower left', xlab='$k_\+ L_\+$', save=outdir+'Ek_ELS.gif')

#--------------------------------------------------------#
#                       2D spectra                       #
#--------------------------------------------------------#
ncfile = netcdf.netcdf_file(input_dir+runname+'.out.fields_section.nc'+restart_num, 'r')   

tt_fld  = np.copy(ncfile.variables['tt' ][:]); tt_fld  = np.delete(tt_fld , ignored_points_fld, axis = 0)
kx_fld  = np.copy(ncfile.variables['kx' ][:])
ky_fld  = np.copy(ncfile.variables['ky' ][:])
kz_fld  = np.copy(ncfile.variables['kz' ][:])

nkx_mid = np.where(np.sign(kx_fld) == -1)[0][0]
nkz_mid = np.where(np.sign(kz_fld) == -1)[0][0]
kx_fld_ = np.concatenate([kx_fld[nkx_mid:], kx_fld[0:nkx_mid]])
kz_fld_ = np.concatenate([kz_fld[nkz_mid:], kz_fld[0:nkz_mid]])

u2_kxy = np.copy(ncfile.variables['u2_kxy'][:]);  u2_kxy = np.delete(u2_kxy, ignored_points_fld, axis = 0); u2_kxy = np.concatenate([u2_kxy[:, :, nkx_mid:], u2_kxy[:, :, 0:nkx_mid]], axis=2)
u2_kyz = np.copy(ncfile.variables['u2_kyz'][:]);  u2_kyz = np.delete(u2_kyz, ignored_points_fld, axis = 0); u2_kyz = np.concatenate([u2_kyz[:, nkz_mid:, :], u2_kyz[:, 0:nkz_mid, :]], axis=1)
u2_kxz = np.copy(ncfile.variables['u2_kxz'][:]);  u2_kxz = np.delete(u2_kxz, ignored_points_fld, axis = 0); u2_kxz = np.concatenate([u2_kxz[:, nkz_mid:, :], u2_kxz[:, 0:nkz_mid, :]], axis=1); u2_kxz = np.concatenate([u2_kxz[:, :, nkx_mid:], u2_kxz[:, :, 0:nkx_mid]], axis=2)
b2_kxy = np.copy(ncfile.variables['b2_kxy'][:]);  b2_kxy = np.delete(b2_kxy, ignored_points_fld, axis = 0); b2_kxy = np.concatenate([b2_kxy[:, :, nkx_mid:], b2_kxy[:, :, 0:nkx_mid]], axis=2)
b2_kyz = np.copy(ncfile.variables['b2_kyz'][:]);  b2_kyz = np.delete(b2_kyz, ignored_points_fld, axis = 0); b2_kyz = np.concatenate([b2_kyz[:, nkz_mid:, :], b2_kyz[:, 0:nkz_mid, :]], axis=1)
b2_kxz = np.copy(ncfile.variables['b2_kxz'][:]);  b2_kxz = np.delete(b2_kxz, ignored_points_fld, axis = 0); b2_kxz = np.concatenate([b2_kxz[:, nkz_mid:, :], b2_kxz[:, 0:nkz_mid, :]], axis=1); b2_kxz = np.concatenate([b2_kxz[:, :, nkx_mid:], b2_kxz[:, :, 0:nkx_mid]], axis=2)

plot_2d(np.log10(u2_kxy[final_idx,:,:]), kx_fld_, ky_fld , xlab=r'$k_x$', ylab=r'$k_y$', title=r'$E_{u} (t = %.2E) $' % tt[final_idx], save=outdir+'u2_kxy.pdf')
plot_2d(np.log10(u2_kyz[final_idx,:,:]), ky_fld , kz_fld_, xlab=r'$k_y$', ylab=r'$k_z$', title=r'$E_{u} (t = %.2E) $' % tt[final_idx], save=outdir+'u2_kyz.pdf')
plot_2d(np.log10(u2_kxz[final_idx,:,:]), kx_fld_, kz_fld_, xlab=r'$k_x$', ylab=r'$k_z$', title=r'$E_{u} (t = %.2E) $' % tt[final_idx], save=outdir+'u2_kxz.pdf')
plot_2d(np.log10(b2_kxy[final_idx,:,:]), kx_fld_, ky_fld , xlab=r'$k_x$', ylab=r'$k_y$', title=r'$E_{b} (t = %.2E) $' % tt[final_idx], save=outdir+'b2_kxy.pdf')
plot_2d(np.log10(b2_kyz[final_idx,:,:]), ky_fld , kz_fld_, xlab=r'$k_y$', ylab=r'$k_z$', title=r'$E_{b} (t = %.2E) $' % tt[final_idx], save=outdir+'b2_kyz.pdf')
plot_2d(np.log10(b2_kxz[final_idx,:,:]), kx_fld_, kz_fld_, xlab=r'$k_x$', ylab=r'$k_z$', title=r'$E_{b} (t = %.2E) $' % tt[final_idx], save=outdir+'b2_kxz.pdf')
