# -*- coding: utf-8 -*-
from load import *
from fft import *
import sys
sys.path.append('../')
from plots import *

print('\nplotting kspectrum\n')
outdir = './fig_kspectrum/'

if nlz == nkz:
  kp_end = np.argmin(np.abs(kpbin - kpbin.max()*2./3.))
else:
  kp_end = kpbin.size - 1
  kp_log_end = kpbin_log.size - 1

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

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ek.txt'  , np.column_stack((kpbin    [:kp_end], 
                                                  b2_bin  [final_idx,:kp_end],
                                                 bx2_bin  [final_idx,:kp_end],
                                                 by2_bin  [final_idx,:kp_end],
                                                 bz2_bin  [final_idx,:kp_end],
                                                 )), fmt='%E')

# kpar
ys = [ 
       kpar_b[final_kpar_idx, 1:kp_log_end], 
       kpbin_log[1:kp_log_end]**(1./3.)/kpbin_log[1]**(1./3.)*kpar_b[final_kpar_idx,1:kp_log_end][0],
     ]
xs = [ 
      kpbin_log[1:kp_log_end], 
      kpbin_log[1:kp_log_end], 
     ]
ls = [ 
        '', 
        'k--', 
     ]
legends = [ 
            r'$k_\|^b$', 
            r'1/3',
          ]

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt_kpar[final_kpar_idx], ylab='', term=True, save=outdir+'kpar.pdf')


# delta b/b0
ys = [ 
       b1_ovr_b0[final_kpar_idx, 1:kp_log_end], 
     ]
xs = [ 
      kpbin_log[1:kp_log_end], 
     ]
ls = [ 
        '', 
     ]
legends = [ 
            r'$\delta b/b_0$', 
          ]

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt_kpar[final_kpar_idx], ylab='', term=True, save=outdir+'b1_ovr_b0.pdf')


# shear AWs vs pseudo AWs
ys = [ 
       b1par2[final_kpar_idx, 1:kp_log_end], 
       b1prp2[final_kpar_idx, 1:kp_log_end], 
       (b1par2/b1prp2)[final_kpar_idx, 1:kp_log_end], 
     ]
xs = [ 
      kpbin_log[1:kp_log_end], 
      kpbin_log[1:kp_log_end], 
      kpbin_log[1:kp_log_end], 
     ]
ls = [ 
        'r-', 
        'r--', 
        'r:', 
     ]
legends = [ 
            r'$\delta B_\|^2$', 
            r'$\delta B_\+^2$', 
            r'$\delta B_\|^2/\delta B_\+^2$', 
          ]

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt_kpar[final_kpar_idx], ylab='', term=True, save=outdir+'shearAW_vs_pseudoAW.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'kpar.txt'  , np.column_stack((kpbin_log[:kp_log_end], 
                                                 kpar_b      [final_kpar_idx,:kp_log_end],
                                                 b1_ovr_b0   [final_kpar_idx,:kp_log_end],
                                                 )), fmt='%E')

np.savetxt(outdir + 'shearAW_vs_pseudoAW.txt'  , np.column_stack((kpbin_log[:kp_log_end], 
                                                 b1par2[final_kpar_idx,:kp_log_end],
                                                 b1prp2[final_kpar_idx,:kp_log_end],
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
tt_fld = np.loadtxt(input_dir+'out2d/time.dat'+restart_num)
nt_fld = tt_fld.size
if nt_fld == 1 : tt_fld = [tt_fld]

nkx_mid = np.where(np.sign(kx) == -1)[0][0]
nkz_mid = np.where(np.sign(kz) == -1)[0][0]
kx_ = np.concatenate([kx[nkx_mid:], kx[0:nkx_mid]])
kz_ = np.concatenate([kz[nkz_mid:], kz[0:nkz_mid]])

b2_kxy = np.transpose(np.fromfile(input_dir+'out2d/b2_kxy_sum_kz.dat'+restart_num).reshape(nt_fld, nky, nkx)); b2_kxy = np.transpose(b2_kxy, axes=(2, 1, 0));  b2_kxy = np.delete(b2_kxy, ignored_points_fld, axis = 0); b2_kxy = np.concatenate([b2_kxy[:, :, nkx_mid:], b2_kxy[:, :, 0:nkx_mid]], axis=2)
b2_kyz = np.transpose(np.fromfile(input_dir+'out2d/b2_kyz_sum_kx.dat'+restart_num).reshape(nt_fld, nkz, nky)); b2_kyz = np.transpose(b2_kyz, axes=(2, 1, 0));  b2_kyz = np.delete(b2_kyz, ignored_points_fld, axis = 0); b2_kyz = np.concatenate([b2_kyz[:, nkz_mid:, :], b2_kyz[:, 0:nkz_mid, :]], axis=1)
b2_kxz = np.transpose(np.fromfile(input_dir+'out2d/b2_kxz_sum_ky.dat'+restart_num).reshape(nt_fld, nkz, nkx)); b2_kxz = np.transpose(b2_kxz, axes=(2, 1, 0));  b2_kxz = np.delete(b2_kxz, ignored_points_fld, axis = 0); b2_kxz = np.concatenate([b2_kxz[:, nkz_mid:, :], b2_kxz[:, 0:nkz_mid, :]], axis=1); b2_kxz = np.concatenate([b2_kxz[:, :, nkx_mid:], b2_kxz[:, :, 0:nkx_mid]], axis=2)

plot_2d(np.log10(b2_kxy[final_idx,1:-1,1:-1]), kx_[1:-1], ky [1:-1], xlab=r'$k_x$', ylab=r'$k_y$', title=r'$E_{b} (t = %.2E) $' % tt[final_idx], save=outdir+'b2_kxy.pdf')
plot_2d(np.log10(b2_kyz[final_idx,1:-1,1:-1]), ky [1:-1], kz_[1:-1], xlab=r'$k_y$', ylab=r'$k_z$', title=r'$E_{b} (t = %.2E) $' % tt[final_idx], save=outdir+'b2_kyz.pdf')
plot_2d(np.log10(b2_kxz[final_idx,1:-1,1:-1]), kx_[1:-1], kz_[1:-1], xlab=r'$k_x$', ylab=r'$k_z$', title=r'$E_{b} (t = %.2E) $' % tt[final_idx], save=outdir+'b2_kxz.pdf')


#--------------------------------------------------------#
#                      plot 1D spectra                   #
#--------------------------------------------------------#
f  = np.sum(b2_kyz[final_idx,:,:], axis=1)[nkz_mid:]
kz = kz_[nkz_mid:]
ys = [ 
       f, 
       kz**(-3.)/kz[0]**(-3.)*f[0],
       kz**(-4.)/kz[0]**(-4.)*f[0],
     ]
xs = [ 
      kz, 
      kz, 
      kz, 
     ]
ls = [ 
        '', 
        'k--', 
        'k:', 
     ]
legends = [ 
            r'$E$',
            r'-3',
            r'-4',
          ]
# if de > 0.0: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=r"$k_z d_\rmi$", legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kz_spectra.pdf')
