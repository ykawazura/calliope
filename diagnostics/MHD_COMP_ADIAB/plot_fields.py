# -*- coding: utf-8 -*-
from load import *
from fft import *
from plots import *

print('\nplotting fields\n')
outdir = './fig_fields/'

#--------------------------------------------------------#
#                   plot final snapshot                  #
#--------------------------------------------------------#
u_r_z0 = np.sqrt(ux_r_z0**2 + uy_r_z0**2 + uz_r_z0**2)
u_r_y0 = np.sqrt(ux_r_y0**2 + uy_r_y0**2 + uz_r_y0**2)
u_r_x0 = np.sqrt(ux_r_x0**2 + uy_r_x0**2 + uz_r_x0**2)
w_r_z0 = np.sqrt(wx_r_z0**2 + wy_r_z0**2 + wz_r_z0**2)
w_r_y0 = np.sqrt(wx_r_y0**2 + wy_r_y0**2 + wz_r_y0**2)
w_r_x0 = np.sqrt(wx_r_x0**2 + wy_r_x0**2 + wz_r_x0**2)
b_r_z0 = np.sqrt(bx_r_z0**2 + by_r_z0**2 + bz_r_z0**2)
b_r_y0 = np.sqrt(bx_r_y0**2 + by_r_y0**2 + bz_r_y0**2)
b_r_x0 = np.sqrt(bx_r_x0**2 + by_r_x0**2 + bz_r_x0**2)
j_r_z0 = np.sqrt(jx_r_z0**2 + jy_r_z0**2 + jz_r_z0**2)
j_r_y0 = np.sqrt(jx_r_y0**2 + jy_r_y0**2 + jz_r_y0**2)
j_r_x0 = np.sqrt(jx_r_x0**2 + jy_r_x0**2 + jz_r_x0**2)
plot_3d(rho_r_z0[final_idx,:,:], rho_r_y0[final_idx,:,:], rho_r_x0[final_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$\rho (t = %.2E)$'     % tt[final_idx], cmp=parula_map, save=outdir+'rho.pdf')
plot_3d(u_r_z0  [final_idx,:,:], u_r_y0  [final_idx,:,:], u_r_x0  [final_idx,:,:], xx, yy, zz, (ux_r_z0[final_idx,:,:], uy_r_z0[final_idx,:,:]), (ux_r_y0[final_idx,:,:], uz_r_y0[final_idx,:,:]), (uy_r_x0[final_idx,:,:], uz_r_x0[final_idx,:,:]), xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\bm{u}| (t = %.2E)$' % tt[final_idx], cmp=parula_map, streamline_density=1, streamline_width=0.5, streamline_color='w', save=outdir+'u.pdf')
plot_3d(b_r_z0  [final_idx,:,:], b_r_y0  [final_idx,:,:], b_r_x0  [final_idx,:,:], xx, yy, zz, (bx_r_z0[final_idx,:,:], by_r_z0[final_idx,:,:]), (bx_r_y0[final_idx,:,:], bz_r_y0[final_idx,:,:]), (by_r_x0[final_idx,:,:], bz_r_x0[final_idx,:,:]), xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\bm{B}| (t = %.2E)$' % tt[final_idx], cmp=parula_map, streamline_density=1, streamline_width=0.5, streamline_color='w', save=outdir+'b.pdf')
plot_3d(prs_r_z0[final_idx,:,:], prs_r_y0[final_idx,:,:], prs_r_x0[final_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$p (t = %.2E)$'        % tt[final_idx], cmp=parula_map, save=outdir+'prs.pdf')
if is2D:
  plot_3d(wz_r_z0[final_idx,:,:], wz_r_y0[final_idx,:,:], wz_r_x0[final_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$(\Curl\bm{u})_z (t = %.2E)$' % tt[final_idx], cmp='RdBu_r', save=outdir+'w.pdf')
  plot_3d(jz_r_z0[final_idx,:,:], jz_r_y0[final_idx,:,:], jz_r_x0[final_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$(\Curl\bm{B})_z (t = %.2E)$' % tt[final_idx], cmp='RdBu_r', save=outdir+'j.pdf')
else:
  plot_3d( w_r_z0[final_idx,:,:],  w_r_y0[final_idx,:,:],  w_r_x0[final_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\Curl\bm{u}| (t = %.2E)$' % tt[final_idx], cmp=parula_map, save=outdir+'w.pdf')
  plot_3d( j_r_z0[final_idx,:,:],  j_r_y0[final_idx,:,:],  j_r_x0[final_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\Curl\bm{B}| (t = %.2E)$' % tt[final_idx], cmp=parula_map, save=outdir+'j.pdf')


#--------------------------------------------------------#
#                       plot movie                       #
#--------------------------------------------------------#
if ismovie:
  # evenly space time
  # nframe = min(200, int(tt[:final_idx].size))
  nframe = tt[:final_idx].size
  idx = np.unique([np.argmin(abs(tt - np.linspace(tt[0], tt[-1], nframe)[i])) for i in range(0, nframe)])
  tt  = tt [idx]

  rho_r_z0 = np.take(rho_r_z0, idx, axis=0)
  rho_r_y0 = np.take(rho_r_y0, idx, axis=0)
  rho_r_x0 = np.take(rho_r_x0, idx, axis=0)

  ux_r_z0  = np.take(ux_r_z0 , idx, axis=0)
  ux_r_y0  = np.take(ux_r_y0 , idx, axis=0)
  ux_r_x0  = np.take(ux_r_x0 , idx, axis=0)
  uy_r_z0  = np.take(uy_r_z0 , idx, axis=0)
  uy_r_y0  = np.take(uy_r_y0 , idx, axis=0)
  uy_r_x0  = np.take(uy_r_x0 , idx, axis=0)
  uz_r_z0  = np.take(uz_r_z0 , idx, axis=0)
  uz_r_y0  = np.take(uz_r_y0 , idx, axis=0)
  uz_r_x0  = np.take(uz_r_x0 , idx, axis=0)
                             
  wx_r_z0  = np.take(wx_r_z0 , idx, axis=0)
  wx_r_y0  = np.take(wx_r_y0 , idx, axis=0)
  wx_r_x0  = np.take(wx_r_x0 , idx, axis=0)
  wy_r_z0  = np.take(wy_r_z0 , idx, axis=0)
  wy_r_y0  = np.take(wy_r_y0 , idx, axis=0)
  wy_r_x0  = np.take(wy_r_x0 , idx, axis=0)
  wz_r_z0  = np.take(wz_r_z0 , idx, axis=0)
  wz_r_y0  = np.take(wz_r_y0 , idx, axis=0)
  wz_r_x0  = np.take(wz_r_x0 , idx, axis=0)
                             
  bx_r_z0  = np.take(bx_r_z0 , idx, axis=0)
  bx_r_y0  = np.take(bx_r_y0 , idx, axis=0)
  bx_r_x0  = np.take(bx_r_x0 , idx, axis=0)
  by_r_z0  = np.take(by_r_z0 , idx, axis=0)
  by_r_y0  = np.take(by_r_y0 , idx, axis=0)
  by_r_x0  = np.take(by_r_x0 , idx, axis=0)
  bz_r_z0  = np.take(bz_r_z0 , idx, axis=0)
  bz_r_y0  = np.take(bz_r_y0 , idx, axis=0)
  bz_r_x0  = np.take(bz_r_x0 , idx, axis=0)
                             
  jx_r_z0  = np.take(jx_r_z0 , idx, axis=0)
  jx_r_y0  = np.take(jx_r_y0 , idx, axis=0)
  jx_r_x0  = np.take(jx_r_x0 , idx, axis=0)
  jy_r_z0  = np.take(jy_r_z0 , idx, axis=0)
  jy_r_y0  = np.take(jy_r_y0 , idx, axis=0)
  jy_r_x0  = np.take(jy_r_x0 , idx, axis=0)
  jz_r_z0  = np.take(jz_r_z0 , idx, axis=0)
  jz_r_y0  = np.take(jz_r_y0 , idx, axis=0)
  jz_r_x0  = np.take(jz_r_x0 , idx, axis=0)

  u_r_z0 = np.sqrt(ux_r_z0**2 + uy_r_z0**2 + uz_r_z0**2)
  u_r_y0 = np.sqrt(ux_r_y0**2 + uy_r_y0**2 + uz_r_y0**2)
  u_r_x0 = np.sqrt(ux_r_x0**2 + uy_r_x0**2 + uz_r_x0**2)
  w_r_z0 = np.sqrt(wx_r_z0**2 + wy_r_z0**2 + wz_r_z0**2)
  w_r_y0 = np.sqrt(wx_r_y0**2 + wy_r_y0**2 + wz_r_y0**2)
  w_r_x0 = np.sqrt(wx_r_x0**2 + wy_r_x0**2 + wz_r_x0**2)
  b_r_z0 = np.sqrt(bx_r_z0**2 + by_r_z0**2 + bz_r_z0**2)
  b_r_y0 = np.sqrt(bx_r_y0**2 + by_r_y0**2 + bz_r_y0**2)
  b_r_x0 = np.sqrt(bx_r_x0**2 + by_r_x0**2 + bz_r_x0**2)
  j_r_z0 = np.sqrt(jx_r_z0**2 + jy_r_z0**2 + jz_r_z0**2)
  j_r_y0 = np.sqrt(jx_r_y0**2 + jy_r_y0**2 + jz_r_y0**2)
  j_r_x0 = np.sqrt(jx_r_x0**2 + jy_r_x0**2 + jz_r_x0**2)
  movie_3d(tt, rho_r_z0, rho_r_y0, rho_r_x0, xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$\rho$'    , cmp=parula_map, save=outdir+'rho_anim.gif')
  movie_3d(tt, u_r_z0  , u_r_y0, u_r_x0    , xx, yy, zz, (ux_r_z0, uy_r_z0), (ux_r_y0, uz_r_y0), (uy_r_x0, uz_r_x0), xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\bm{u}|$', cmp=parula_map, save=outdir+'u_anim.gif')
  movie_3d(tt, b_r_z0  , b_r_y0, b_r_x0    , xx, yy, zz, (bx_r_z0, by_r_z0), (bx_r_y0, bz_r_y0), (by_r_x0, bz_r_x0), xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\bm{B}|$', cmp=parula_map, save=outdir+'b_anim.gif')
  movie_3d(tt, prs_r_z0, prs_r_y0, prs_r_x0, xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$p$'       , cmp=parula_map, save=outdir+'prs_anim.gif')
  if is2D:
    movie_3d(tt, wz_r_z0, wz_r_y0, wz_r_x0, xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$(\Curl\bm{u})_z$', cmp='RdBu_r', save=outdir+'w_anim.gif')
    movie_3d(tt, jz_r_z0, jz_r_y0, jz_r_x0, xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$(\Curl\bm{B})_z$', cmp='RdBu_r', save=outdir+'j_anim.gif')
  else:
    movie_3d(tt,  w_r_z0,  w_r_y0,  w_r_x0, xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\Curl\bm{u}|$', cmp=parula_map, save=outdir+'w_anim.gif')
    movie_3d(tt,  j_r_z0,  j_r_y0,  j_r_x0, xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\Curl\bm{B}|$', cmp=parula_map, save=outdir+'j_anim.gif')

