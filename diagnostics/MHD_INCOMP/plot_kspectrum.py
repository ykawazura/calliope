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
kp_log_kpar_end = kpbin_log_kpar.size - 1

# If shear is on, add a vertical line indicating the fastest-MRI mode
if shear_flg == 1 and np.max(b0[0]) > 1e-10:
  b00 = np.max(b0[0])
  k_mri = np.sqrt(15.)/4./b00

def add_vertical_line(xs, ys, ls, legends):
  xs.append([k_mri, k_mri])
  ys.append([np.min([u2_bin[final_idx, 1:kp_end].min(), b2_bin[final_idx, 1:kp_end].min()]), 
             np.max([u2_bin[final_idx, 1:kp_end].max(), b2_bin[final_idx, 1:kp_end].max()])])
  ls.append('k:')
  legends.append(r'$k = k_\mr{MRI} = (\sqrt{15}/4)\Omega/v_\mathrm{A}$')

  return xs, ys, ls, legends

def fold_in_kz(u):
  unew = [(u[:, i, :] + u[:, -i, :])/2 for i in np.arange(0, int(nkz/2)+1)]
  return np.transpose(np.asarray(unew), axes=[1,0,2])

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

if shear_flg == 1 and np.max(b0[0]) > 1e-10: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

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

if shear_flg == 1 and np.max(b0[0]) > 1e-10: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

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

if shear_flg == 1 and np.max(b0[0]) > 1e-10: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra_b.pdf')


# Elsasser fields
ys = [ 
       zp2_bin[final_idx, 1:kp_end], 
       zm2_bin[final_idx, 1:kp_end], 
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*zp2_bin[final_idx,1:kp_end][0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*zp2_bin[final_idx,1:kp_end][0],
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

if shear_flg == 1 and np.max(b0[0]) > 1e-10: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra_ELS.pdf')


# Dissipation
ys = [ 
       u2dissip_bin[final_idx, 1:kp_end], 
       b2dissip_bin[final_idx, 1:kp_end], 
     ]
xs = [ 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
     ]
ls = [ 
        '', 
        '', 
     ]
legends = [ 
            r'$D_\mr{u}$', 
            r'$D_\mr{B}$',
          ]

if shear_flg == 1 and np.max(b0[0]) > 1e-10: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_dissip.pdf')


# MRI Injection
ys = [ 
       p_re_bin[final_idx, 1:kp_end], 
       p_ma_bin[final_idx, 1:kp_end], 
     ]
xs = [ 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
     ]
ls = [ 
        '', 
        '', 
     ]
legends = [ 
            r'$P_\mr{Re}$', 
            r'$P_\mr{M}$',
          ]

if shear_flg == 1 and np.max(b0[0]) > 1e-10: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_MRI.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ek.txt'  , np.column_stack((kpbin       [:kp_end], 
                                                  u2_bin     [final_idx,:kp_end],
                                                 ux2_bin     [final_idx,:kp_end],
                                                 uy2_bin     [final_idx,:kp_end],
                                                 uz2_bin     [final_idx,:kp_end],
                                                  b2_bin     [final_idx,:kp_end],
                                                 bx2_bin     [final_idx,:kp_end],
                                                 by2_bin     [final_idx,:kp_end],
                                                 bz2_bin     [final_idx,:kp_end],
                                                 zp2_bin     [final_idx,:kp_end],
                                                 zm2_bin     [final_idx,:kp_end],
                                                 u2dissip_bin[final_idx,:kp_end],
                                                 b2dissip_bin[final_idx,:kp_end],
                                                 p_re_bin    [final_idx,:kp_end],
                                                 p_ma_bin    [final_idx,:kp_end],
                                                 )), fmt='%E')

# kpar
try:
  ys = [ 
         kpar_b[final_kpar_idx, 1:kp_log_kpar_end], 
         kpar_u[final_kpar_idx, 1:kp_log_kpar_end], 
         kpbin_log_kpar[1:kp_log_kpar_end]**(2./3.)/kpbin_log_kpar[1]**(2./3.)*kpar_b[final_kpar_idx,1:kp_log_kpar_end][0],
       ]
  xs = [ 
        kpbin_log_kpar[1:kp_log_kpar_end], 
        kpbin_log_kpar[1:kp_log_kpar_end], 
        kpbin_log_kpar[1:kp_log_kpar_end], 
       ]
  ls = [ 
          '', 
          '', 
          'k--', 
       ]
  legends = [ 
              r'$k_\|^b$', 
              r'$k_\|^u$', 
              r'2/3',
            ]

  if shear_flg == 1 and np.max(b0[0]) > 1e-10: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

  plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt_kpar[final_kpar_idx], ylab='', term=True, save=outdir+'kpar.pdf')


  # delta b/b0
  ys = [ 
         b1_ovr_b0[final_kpar_idx, 1:kp_log_kpar_end], 
       ]
  xs = [ 
        kpbin_log_kpar[1:kp_log_kpar_end], 
       ]
  ls = [ 
          '', 
       ]
  legends = [ 
              r'$\delta b/b_0$', 
            ]

  if shear_flg == 1 and np.max(b0[0]) > 1e-10: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

  plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt_kpar[final_kpar_idx], ylab='', term=True, save=outdir+'b1_ovr_b0.pdf')


  # shear AWs vs pseudo AWs
  ys = [ 
         b1par2[final_kpar_idx, 1:kp_log_kpar_end], 
         b1prp2[final_kpar_idx, 1:kp_log_kpar_end], 
         (b1par2/b1prp2)[final_kpar_idx, 1:kp_log_kpar_end], 
         u1par2[final_kpar_idx, 1:kp_log_kpar_end], 
         u1prp2[final_kpar_idx, 1:kp_log_kpar_end], 
         (u1par2/u1prp2)[final_kpar_idx, 1:kp_log_kpar_end], 
       ]
  xs = [ 
        kpbin_log_kpar[1:kp_log_kpar_end], 
        kpbin_log_kpar[1:kp_log_kpar_end], 
        kpbin_log_kpar[1:kp_log_kpar_end], 
        kpbin_log_kpar[1:kp_log_kpar_end], 
        kpbin_log_kpar[1:kp_log_kpar_end], 
        kpbin_log_kpar[1:kp_log_kpar_end], 
       ]
  ls = [ 
          'r-', 
          'r--', 
          'r:', 
          'b-', 
          'b--', 
          'b:', 
       ]
  legends = [ 
              r'$\delta B_\|^2$', 
              r'$\delta B_\+^2$', 
              r'$\delta B_\|^2/\delta B_\+^2$', 
              r'$u_\|^2$', 
              r'$u_\+^2$', 
              r'$u_\|^2/u_\+^2$', 
            ]

  if shear_flg == 1 and np.max(b0[0]) > 1e-10: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

  plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt_kpar[final_kpar_idx], ylab='', term=True, save=outdir+'shearAW_vs_pseudoAW.pdf')

  #------------------#
  #   output ascii   #
  #------------------#
  np.savetxt(outdir + 'tt_kpar.txt'  , tt_kpar, fmt='%E')
  np.savetxt(outdir + 'kpar.txt'  , np.column_stack((kpbin_log_kpar[:kp_log_kpar_end], 
                                                   kpar_b      [final_kpar_idx,:kp_log_kpar_end],
                                                   kpar_u      [final_kpar_idx,:kp_log_kpar_end],
                                                   b1_ovr_b0   [final_kpar_idx,:kp_log_kpar_end],
                                                   )), fmt='%E')

  np.savetxt(outdir + 'shearAW_vs_pseudoAW.txt'  , np.column_stack((kpbin_log_kpar[:kp_log_kpar_end], 
                                                   b1par2[final_kpar_idx,:kp_log_kpar_end],
                                                   b1prp2[final_kpar_idx,:kp_log_kpar_end],
                                                   u1par2[final_kpar_idx,:kp_log_kpar_end],
                                                   u1prp2[final_kpar_idx,:kp_log_kpar_end],
                                                   )), fmt='%E')
except:
  pass

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
# summed over kx or ky or kz
tt_2d = np.transpose(np.loadtxt(input_dir+'out2d'+restart_num+'/time.dat'))[0]
nt_2d = tt_2d.size
if nt_2d == 1 : tt_2d = [tt_2d]

nkx_mid = np.where(np.sign(kx) == -1)[0][0]
nkz_mid = np.where(np.sign(kz) == -1)[0][0]
kx_ = np.concatenate([kx[nkx_mid:], kx[0:nkx_mid]])
kz_ = np.concatenate([kz[nkz_mid:], kz[0:nkz_mid]])

u2_kxy_sum_kz = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/u2_kxy_sum_kz.dat').reshape(nt_2d, nky, nkx)); u2_kxy_sum_kz = np.transpose(u2_kxy_sum_kz, axes=(2, 1, 0));  u2_kxy_sum_kz = np.delete(u2_kxy_sum_kz, ignored_points_2d, axis = 0); u2_kxy_sum_kz = np.concatenate([u2_kxy_sum_kz[:, :, nkx_mid:], u2_kxy_sum_kz[:, :, 0:nkx_mid]], axis=2)
u2_kyz_sum_kx = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/u2_kyz_sum_kx.dat').reshape(nt_2d, nkz, nky)); u2_kyz_sum_kx = np.transpose(u2_kyz_sum_kx, axes=(2, 1, 0));  u2_kyz_sum_kx = np.delete(u2_kyz_sum_kx, ignored_points_2d, axis = 0); u2_kyz_sum_kx = np.concatenate([u2_kyz_sum_kx[:, nkz_mid:, :], u2_kyz_sum_kx[:, 0:nkz_mid, :]], axis=1)
u2_kxz_sum_ky = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/u2_kxz_sum_ky.dat').reshape(nt_2d, nkz, nkx)); u2_kxz_sum_ky = np.transpose(u2_kxz_sum_ky, axes=(2, 1, 0));  u2_kxz_sum_ky = np.delete(u2_kxz_sum_ky, ignored_points_2d, axis = 0); u2_kxz_sum_ky = np.concatenate([u2_kxz_sum_ky[:, nkz_mid:, :], u2_kxz_sum_ky[:, 0:nkz_mid, :]], axis=1); u2_kxz_sum_ky = np.concatenate([u2_kxz_sum_ky[:, :, nkx_mid:], u2_kxz_sum_ky[:, :, 0:nkx_mid]], axis=2)
b2_kxy_sum_kz = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/b2_kxy_sum_kz.dat').reshape(nt_2d, nky, nkx)); b2_kxy_sum_kz = np.transpose(b2_kxy_sum_kz, axes=(2, 1, 0));  b2_kxy_sum_kz = np.delete(b2_kxy_sum_kz, ignored_points_2d, axis = 0); b2_kxy_sum_kz = np.concatenate([b2_kxy_sum_kz[:, :, nkx_mid:], b2_kxy_sum_kz[:, :, 0:nkx_mid]], axis=2)
b2_kyz_sum_kx = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/b2_kyz_sum_kx.dat').reshape(nt_2d, nkz, nky)); b2_kyz_sum_kx = np.transpose(b2_kyz_sum_kx, axes=(2, 1, 0));  b2_kyz_sum_kx = np.delete(b2_kyz_sum_kx, ignored_points_2d, axis = 0); b2_kyz_sum_kx = np.concatenate([b2_kyz_sum_kx[:, nkz_mid:, :], b2_kyz_sum_kx[:, 0:nkz_mid, :]], axis=1)
b2_kxz_sum_ky = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/b2_kxz_sum_ky.dat').reshape(nt_2d, nkz, nkx)); b2_kxz_sum_ky = np.transpose(b2_kxz_sum_ky, axes=(2, 1, 0));  b2_kxz_sum_ky = np.delete(b2_kxz_sum_ky, ignored_points_2d, axis = 0); b2_kxz_sum_ky = np.concatenate([b2_kxz_sum_ky[:, nkz_mid:, :], b2_kxz_sum_ky[:, 0:nkz_mid, :]], axis=1); b2_kxz_sum_ky = np.concatenate([b2_kxz_sum_ky[:, :, nkx_mid:], b2_kxz_sum_ky[:, :, 0:nkx_mid]], axis=2)

plot_2d(np.log10(u2_kxy_sum_kz[final_idx,1:-1,1:-1]), kx_[1:-1], ky [1:-1], xlab=r'$k_x$', ylab=r'$k_y$', title=r'$E_u (t = %.2E) $' % tt_2d[final_2d_idx], contour=True, aspect='auto', save=outdir+'u2_kxy.pdf')
plot_2d(np.log10(u2_kyz_sum_kx[final_idx,1:-1,1:-1]), ky [1:-1], kz_[1:-1], xlab=r'$k_y$', ylab=r'$k_z$', title=r'$E_u (t = %.2E) $' % tt_2d[final_2d_idx], contour=True, aspect='auto', save=outdir+'u2_kyz.pdf')
plot_2d(np.log10(u2_kxz_sum_ky[final_idx,1:-1,1:-1]), kx_[1:-1], kz_[1:-1], xlab=r'$k_x$', ylab=r'$k_z$', title=r'$E_u (t = %.2E) $' % tt_2d[final_2d_idx], contour=True, aspect='auto', save=outdir+'u2_kxz.pdf')
plot_2d(np.log10(b2_kxy_sum_kz[final_idx,1:-1,1:-1]), kx_[1:-1], ky [1:-1] , xlab=r'$k_x$', ylab=r'$k_y$', title=r'$E_B (t = %.2E) $' % tt_2d[final_2d_idx], contour=True, aspect='auto', save=outdir+'b2_kxy.pdf')
plot_2d(np.log10(b2_kyz_sum_kx[final_idx,1:-1,1:-1]), ky [1:-1], kz_[1:-1], xlab=r'$k_y$', ylab=r'$k_z$', title=r'$E_B (t = %.2E) $' % tt_2d[final_2d_idx], contour=True, aspect='auto', save=outdir+'b2_kyz.pdf')
plot_2d(np.log10(b2_kxz_sum_ky[final_idx,1:-1,1:-1]), kx_[1:-1], kz_[1:-1], xlab=r'$k_x$', ylab=r'$k_z$', title=r'$E_B (t = %.2E) $' % tt_2d[final_2d_idx], contour=True, aspect='auto', save=outdir+'b2_kxz.pdf')


# bin over (kx, ky) leaving kz or (kx, kz) leaving ky
u2_kxy_vs_kz = np.fromfile(input_dir+'out2d'+restart_num+'/u2_kxy_vs_kz.dat').reshape(nt_2d, nkz, nkpolar); u2_kxy_vs_kz = fold_in_kz(u2_kxy_vs_kz)
u2_kxz_vs_ky = np.fromfile(input_dir+'out2d'+restart_num+'/u2_kxz_vs_ky.dat').reshape(nt_2d, nky, nkpolar)
b2_kxy_vs_kz = np.fromfile(input_dir+'out2d'+restart_num+'/b2_kxy_vs_kz.dat').reshape(nt_2d, nkz, nkpolar); b2_kxy_vs_kz = fold_in_kz(b2_kxy_vs_kz)
b2_kxz_vs_ky = np.fromfile(input_dir+'out2d'+restart_num+'/b2_kxz_vs_ky.dat').reshape(nt_2d, nky, nkpolar)

plot_log2d(u2_kxy_vs_kz[final_idx,1:int(nkz/2),1:kp_end], kpbin[1:kp_end], kz[1:int(nkz/2)] , xlab=r'$k_\+$', ylab=r'$k_z$', title=r'$E_u (t = %.2E) $' % tt_2d[final_2d_idx], contour=True, save=outdir+'u2_kxy_vs_kz.pdf')
plot_log2d(b2_kxy_vs_kz[final_idx,1:int(nkz/2),1:kp_end], kpbin[1:kp_end], kz[1:int(nkz/2)] , xlab=r'$k_\+$', ylab=r'$k_z$', title=r'$E_u (t = %.2E) $' % tt_2d[final_2d_idx], contour=True, save=outdir+'b2_kxy_vs_kz.pdf')

plot_log2d(u2_kxz_vs_ky[final_idx,1:nky - 1   ,1:kp_end], kpbin[1:kp_end], ky[1:nky - 1   ] , xlab=r'$k_\+$', ylab=r'$k_y$', title=r'$E_B (t = %.2E) $' % tt_2d[final_2d_idx], contour=True, save=outdir+'u2_kxz_vs_ky.pdf')
plot_log2d(b2_kxz_vs_ky[final_idx,1:nky - 1   ,1:kp_end], kpbin[1:kp_end], ky[1:nky - 1   ] , xlab=r'$k_\+$', ylab=r'$k_y$', title=r'$E_B (t = %.2E) $' % tt_2d[final_2d_idx], contour=True, save=outdir+'b2_kxz_vs_ky.pdf')


ys = [ 
       u2_bin[final_idx, 1:kp_end], 
       np.sum(u2_kxy_vs_kz[final_idx,:,1:kp_end], axis=0), 
       np.sum(u2_kxz_vs_ky[final_idx,:,1:kp_end], axis=0), 
       b2_bin[final_idx, 1:kp_end], 
       np.sum(b2_kxy_vs_kz[final_idx,:,1:kp_end], axis=0), 
       np.sum(b2_kxz_vs_ky[final_idx,:,1:kp_end], axis=0), 
     ]
xs = [ 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
     ]
ls = [ 
        'r-', 
        'r--', 
        'r.', 
        'b-', 
        'b--', 
        'b.', 
     ]
legends = [ 
            r'$E_u$', 
            r'$\sum_{k_z}E_u(k_z, \sqrt{k_x^2 + k_y^2})$',
            r'$\sum_{k_y}E_u(k_y, \sqrt{k_x^2 + k_z^2})$',
            r'$E_B$', 
            r'$\sum_{k_z}E_B(k_z, \sqrt{k_x^2 + k_y^2})$',
            r'$\sum_{k_y}E_B(k_y, \sqrt{k_x^2 + k_z^2})$',
          ]

if shear_flg == 1 and np.max(b0[0]) > 1e-10: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt_2d[final_2d_idx], ylab='', term=True, save=outdir+'2d_3d_bin_check.pdf')

#------------------#
#   output data    #
#------------------#
from scipy.io import savemat

savemat(outdir + '2d_bin_in_3d_grid' , {'kpbin':kpbin, 'ky':ky, 'kz':kz[:int(nkz/2)+1]})
savemat(outdir + '2d_bin_in_3d' , {
                            'tt'           :tt_2d[final_2d_idx],
                            'u2_kxy_vs_kz' :u2_kxy_vs_kz[final_2d_idx,:,:],
                            'u2_kxz_vs_ky' :u2_kxz_vs_ky[final_2d_idx,:,:],
                            'b2_kxy_vs_kz' :b2_kxy_vs_kz[final_2d_idx,:,:],
                            'b2_kxz_vs_ky' :b2_kxz_vs_ky[final_2d_idx,:,:],
                          })
