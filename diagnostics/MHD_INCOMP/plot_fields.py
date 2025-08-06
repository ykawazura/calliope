# -*- coding: utf-8 -*-
from load import *
from fft import *
import sys
sys.path.append('../')
from plots import *

print('\nplotting fields\n')
outdir = './fig_fields/'

tt_2d = np.transpose(np.loadtxt(input_dir+'out2d'+restart_num+'/time.dat'))[0]
nt_2d = tt_2d.size
if nt_2d == 1 : tt_2d = [tt_2d]

ux_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/ux_r_z0.dat').reshape(nt_2d, nly, nlx)); ux_r_z0 = np.transpose(ux_r_z0, axes=(2, 1, 0))
ux_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/ux_r_x0.dat').reshape(nt_2d, nlz, nly)); ux_r_x0 = np.transpose(ux_r_x0, axes=(2, 1, 0))
ux_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/ux_r_y0.dat').reshape(nt_2d, nlz, nlx)); ux_r_y0 = np.transpose(ux_r_y0, axes=(2, 1, 0))
                                                                  
uy_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/uy_r_z0.dat').reshape(nt_2d, nly, nlx)); uy_r_z0 = np.transpose(uy_r_z0, axes=(2, 1, 0))
uy_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/uy_r_x0.dat').reshape(nt_2d, nlz, nly)); uy_r_x0 = np.transpose(uy_r_x0, axes=(2, 1, 0))
uy_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/uy_r_y0.dat').reshape(nt_2d, nlz, nlx)); uy_r_y0 = np.transpose(uy_r_y0, axes=(2, 1, 0))
                                                                  
uz_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/uz_r_z0.dat').reshape(nt_2d, nly, nlx)); uz_r_z0 = np.transpose(uz_r_z0, axes=(2, 1, 0))
uz_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/uz_r_x0.dat').reshape(nt_2d, nlz, nly)); uz_r_x0 = np.transpose(uz_r_x0, axes=(2, 1, 0))
uz_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/uz_r_y0.dat').reshape(nt_2d, nlz, nlx)); uz_r_y0 = np.transpose(uz_r_y0, axes=(2, 1, 0))
                                                                  
wx_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/wx_r_z0.dat').reshape(nt_2d, nly, nlx)); wx_r_z0 = np.transpose(wx_r_z0, axes=(2, 1, 0))
wx_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/wx_r_x0.dat').reshape(nt_2d, nlz, nly)); wx_r_x0 = np.transpose(wx_r_x0, axes=(2, 1, 0))
wx_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/wx_r_y0.dat').reshape(nt_2d, nlz, nlx)); wx_r_y0 = np.transpose(wx_r_y0, axes=(2, 1, 0))
                                                                                                                                                            
wy_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/wy_r_z0.dat').reshape(nt_2d, nly, nlx)); wy_r_z0 = np.transpose(wy_r_z0, axes=(2, 1, 0))
wy_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/wy_r_x0.dat').reshape(nt_2d, nlz, nly)); wy_r_x0 = np.transpose(wy_r_x0, axes=(2, 1, 0))
wy_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/wy_r_y0.dat').reshape(nt_2d, nlz, nlx)); wy_r_y0 = np.transpose(wy_r_y0, axes=(2, 1, 0))
                                                                                                                                                            
wz_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/wz_r_z0.dat').reshape(nt_2d, nly, nlx)); wz_r_z0 = np.transpose(wz_r_z0, axes=(2, 1, 0))
wz_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/wz_r_x0.dat').reshape(nt_2d, nlz, nly)); wz_r_x0 = np.transpose(wz_r_x0, axes=(2, 1, 0))
wz_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/wz_r_y0.dat').reshape(nt_2d, nlz, nlx)); wz_r_y0 = np.transpose(wz_r_y0, axes=(2, 1, 0))
                                                                  
bx_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/bx_r_z0.dat').reshape(nt_2d, nly, nlx)); bx_r_z0 = np.transpose(bx_r_z0, axes=(2, 1, 0))
bx_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/bx_r_x0.dat').reshape(nt_2d, nlz, nly)); bx_r_x0 = np.transpose(bx_r_x0, axes=(2, 1, 0))
bx_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/bx_r_y0.dat').reshape(nt_2d, nlz, nlx)); bx_r_y0 = np.transpose(bx_r_y0, axes=(2, 1, 0))
                                                                                                                                                            
by_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/by_r_z0.dat').reshape(nt_2d, nly, nlx)); by_r_z0 = np.transpose(by_r_z0, axes=(2, 1, 0))
by_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/by_r_x0.dat').reshape(nt_2d, nlz, nly)); by_r_x0 = np.transpose(by_r_x0, axes=(2, 1, 0))
by_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/by_r_y0.dat').reshape(nt_2d, nlz, nlx)); by_r_y0 = np.transpose(by_r_y0, axes=(2, 1, 0))
                                                                                                                                                            
bz_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/bz_r_z0.dat').reshape(nt_2d, nly, nlx)); bz_r_z0 = np.transpose(bz_r_z0, axes=(2, 1, 0))
bz_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/bz_r_x0.dat').reshape(nt_2d, nlz, nly)); bz_r_x0 = np.transpose(bz_r_x0, axes=(2, 1, 0))
bz_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/bz_r_y0.dat').reshape(nt_2d, nlz, nlx)); bz_r_y0 = np.transpose(bz_r_y0, axes=(2, 1, 0))
                                                                  
jx_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/jx_r_z0.dat').reshape(nt_2d, nly, nlx)); jx_r_z0 = np.transpose(jx_r_z0, axes=(2, 1, 0))
jx_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/jx_r_x0.dat').reshape(nt_2d, nlz, nly)); jx_r_x0 = np.transpose(jx_r_x0, axes=(2, 1, 0))
jx_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/jx_r_y0.dat').reshape(nt_2d, nlz, nlx)); jx_r_y0 = np.transpose(jx_r_y0, axes=(2, 1, 0))
                                                                                                                                                            
jy_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/jy_r_z0.dat').reshape(nt_2d, nly, nlx)); jy_r_z0 = np.transpose(jy_r_z0, axes=(2, 1, 0))
jy_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/jy_r_x0.dat').reshape(nt_2d, nlz, nly)); jy_r_x0 = np.transpose(jy_r_x0, axes=(2, 1, 0))
jy_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/jy_r_y0.dat').reshape(nt_2d, nlz, nlx)); jy_r_y0 = np.transpose(jy_r_y0, axes=(2, 1, 0))
                                                                                                                                                            
jz_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/jz_r_z0.dat').reshape(nt_2d, nly, nlx)); jz_r_z0 = np.transpose(jz_r_z0, axes=(2, 1, 0))
jz_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/jz_r_x0.dat').reshape(nt_2d, nlz, nly)); jz_r_x0 = np.transpose(jz_r_x0, axes=(2, 1, 0))
jz_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/jz_r_y0.dat').reshape(nt_2d, nlz, nlx)); jz_r_y0 = np.transpose(jz_r_y0, axes=(2, 1, 0))

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
plot_3d(u_r_z0[final_2d_idx,:,:], u_r_y0[final_2d_idx,:,:], u_r_x0[final_2d_idx,:,:], xx, yy, zz, (ux_r_z0[final_2d_idx,:,:], uy_r_z0[final_2d_idx,:,:]), (ux_r_y0[final_2d_idx,:,:], uz_r_y0[final_2d_idx,:,:]), (uy_r_x0[final_2d_idx,:,:], uz_r_x0[final_2d_idx,:,:]), xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\bm{u}| (t = %.2E)$' % tt_2d[final_2d_idx], cmp='sequential', streamline_density=1, streamline_width=0.5, streamline_color='w', save=outdir+'u.pdf')
plot_3d(b_r_z0[final_2d_idx,:,:], b_r_y0[final_2d_idx,:,:], b_r_x0[final_2d_idx,:,:], xx, yy, zz, (bx_r_z0[final_2d_idx,:,:], by_r_z0[final_2d_idx,:,:]), (bx_r_y0[final_2d_idx,:,:], bz_r_y0[final_2d_idx,:,:]), (by_r_x0[final_2d_idx,:,:], bz_r_x0[final_2d_idx,:,:]), xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\bm{B}| (t = %.2E)$' % tt_2d[final_2d_idx], cmp='sequential', streamline_density=1, streamline_width=0.5, streamline_color='w', save=outdir+'b.pdf')
plot_3d(bx_r_z0[final_2d_idx,:,:], bx_r_y0[final_2d_idx,:,:], bx_r_x0[final_2d_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$B_x (t = %.2E)$' % tt_2d[final_2d_idx], cmp='diverging', streamline_density=1, streamline_width=0.5, streamline_color='w', save=outdir+'bx.pdf')
plot_3d(by_r_z0[final_2d_idx,:,:], by_r_y0[final_2d_idx,:,:], by_r_x0[final_2d_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$B_y (t = %.2E)$' % tt_2d[final_2d_idx], cmp='diverging', streamline_density=1, streamline_width=0.5, streamline_color='w', save=outdir+'by.pdf')
plot_3d(bz_r_z0[final_2d_idx,:,:], bz_r_y0[final_2d_idx,:,:], bz_r_x0[final_2d_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$B_z (t = %.2E)$' % tt_2d[final_2d_idx], cmp='diverging', streamline_density=1, streamline_width=0.5, streamline_color='w', save=outdir+'bz.pdf')
plot_3d(ux_r_z0[final_2d_idx,:,:], ux_r_y0[final_2d_idx,:,:], ux_r_x0[final_2d_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$u_x (t = %.2E)$' % tt_2d[final_2d_idx], cmp='diverging', streamline_density=1, streamline_width=0.5, streamline_color='w', save=outdir+'ux.pdf')
plot_3d(uy_r_z0[final_2d_idx,:,:], uy_r_y0[final_2d_idx,:,:], uy_r_x0[final_2d_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$u_y (t = %.2E)$' % tt_2d[final_2d_idx], cmp='diverging', streamline_density=1, streamline_width=0.5, streamline_color='w', save=outdir+'uy.pdf')
plot_3d(uz_r_z0[final_2d_idx,:,:], uz_r_y0[final_2d_idx,:,:], uz_r_x0[final_2d_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$u_z (t = %.2E)$' % tt_2d[final_2d_idx], cmp='diverging', streamline_density=1, streamline_width=0.5, streamline_color='w', save=outdir+'uz.pdf')
if is2D:
  plot_3d(wz_r_z0[final_2d_idx,:,:], wz_r_y0[final_2d_idx,:,:], wz_r_x0[final_2d_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$(\Curl\bm{u})_z (t = %.2E)$' % tt_2d[final_2d_idx], cmp='RdBu_r', save=outdir+'w.pdf')
  plot_3d(jz_r_z0[final_2d_idx,:,:], jz_r_y0[final_2d_idx,:,:], jz_r_x0[final_2d_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$(\Curl\bm{B})_z (t = %.2E)$' % tt_2d[final_2d_idx], cmp='RdBu_r', save=outdir+'j.pdf')
else:
  plot_3d( w_r_z0[final_2d_idx,:,:],  w_r_y0[final_2d_idx,:,:],  w_r_x0[final_2d_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\Curl\bm{u}| (t = %.2E)$' % tt_2d[final_2d_idx], cmp='sequential', save=outdir+'w.pdf')
  plot_3d( j_r_z0[final_2d_idx,:,:],  j_r_y0[final_2d_idx,:,:],  j_r_x0[final_2d_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\Curl\bm{B}| (t = %.2E)$' % tt_2d[final_2d_idx], cmp='sequential', save=outdir+'j.pdf')


#--------------------------------------------------------#
#                       plot movie                       #
#--------------------------------------------------------#
if ismovie:
  # evenly space time
  nframe = min(200, int(tt_2d[:final_2d_idx].size))
  idx = np.unique([np.argmin(abs(tt_2d - np.linspace(tt_2d[0], tt_2d[-1], nframe)[i])) for i in range(0, nframe)])
  tt_2d  = tt_2d [idx]

  ux_r_z0 = np.take(ux_r_z0, idx, axis=0)
  ux_r_y0 = np.take(ux_r_y0, idx, axis=0)
  ux_r_x0 = np.take(ux_r_x0, idx, axis=0)
  uy_r_z0 = np.take(uy_r_z0, idx, axis=0)
  uy_r_y0 = np.take(uy_r_y0, idx, axis=0)
  uy_r_x0 = np.take(uy_r_x0, idx, axis=0)
  uz_r_z0 = np.take(uz_r_z0, idx, axis=0)
  uz_r_y0 = np.take(uz_r_y0, idx, axis=0)
  uz_r_x0 = np.take(uz_r_x0, idx, axis=0)

  wx_r_z0 = np.take(wx_r_z0, idx, axis=0)
  wx_r_y0 = np.take(wx_r_y0, idx, axis=0)
  wx_r_x0 = np.take(wx_r_x0, idx, axis=0)
  wy_r_z0 = np.take(wy_r_z0, idx, axis=0)
  wy_r_y0 = np.take(wy_r_y0, idx, axis=0)
  wy_r_x0 = np.take(wy_r_x0, idx, axis=0)
  wz_r_z0 = np.take(wz_r_z0, idx, axis=0)
  wz_r_y0 = np.take(wz_r_y0, idx, axis=0)
  wz_r_x0 = np.take(wz_r_x0, idx, axis=0)

  bx_r_z0 = np.take(bx_r_z0, idx, axis=0)
  bx_r_y0 = np.take(bx_r_y0, idx, axis=0)
  bx_r_x0 = np.take(bx_r_x0, idx, axis=0)
  by_r_z0 = np.take(by_r_z0, idx, axis=0)
  by_r_y0 = np.take(by_r_y0, idx, axis=0)
  by_r_x0 = np.take(by_r_x0, idx, axis=0)
  bz_r_z0 = np.take(bz_r_z0, idx, axis=0)
  bz_r_y0 = np.take(bz_r_y0, idx, axis=0)
  bz_r_x0 = np.take(bz_r_x0, idx, axis=0)

  jx_r_z0 = np.take(jx_r_z0, idx, axis=0)
  jx_r_y0 = np.take(jx_r_y0, idx, axis=0)
  jx_r_x0 = np.take(jx_r_x0, idx, axis=0)
  jy_r_z0 = np.take(jy_r_z0, idx, axis=0)
  jy_r_y0 = np.take(jy_r_y0, idx, axis=0)
  jy_r_x0 = np.take(jy_r_x0, idx, axis=0)
  jz_r_z0 = np.take(jz_r_z0, idx, axis=0)
  jz_r_y0 = np.take(jz_r_y0, idx, axis=0)
  jz_r_x0 = np.take(jz_r_x0, idx, axis=0)

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
  movie_3d(tt_2d, u_r_z0  , u_r_y0, u_r_x0    , xx, yy, zz, (ux_r_z0, uy_r_z0), (ux_r_y0, uz_r_y0), (uy_r_x0, uz_r_x0), xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\bm{u}|$', cmp='sequential', save=outdir+'u_anim.gif')
  movie_3d(tt_2d, b_r_z0  , b_r_y0, b_r_x0    , xx, yy, zz, (bx_r_z0, by_r_z0), (bx_r_y0, bz_r_y0), (by_r_x0, bz_r_x0), xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\bm{B}|$', cmp='sequential', save=outdir+'b_anim.gif')
  if is2D:
    movie_3d(tt_2d, wz_r_z0, wz_r_y0, wz_r_x0, xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$(\Curl\bm{u})_z$', cmp='RdBu_r', save=outdir+'w_anim.gif')
    movie_3d(tt_2d, jz_r_z0, jz_r_y0, jz_r_x0, xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$(\Curl\bm{B})_z$', cmp='RdBu_r', save=outdir+'j_anim.gif')
  else:
    movie_3d(tt_2d,  w_r_z0,  w_r_y0,  w_r_x0, xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\Curl\bm{u}|$', cmp='sequential', save=outdir+'w_anim.gif')
    movie_3d(tt_2d,  j_r_z0,  j_r_y0,  j_r_x0, xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\Curl\bm{B}|$', cmp='sequential', save=outdir+'j_anim.gif')


#------------------#
#   output data    #
#------------------#
from scipy.io import savemat

savemat(outdir + 'grid' , {'xx':xx, 'yy':yy, 'zz':zz})
savemat(outdir + 'u_r' , {
                            'tt'      :tt[final_2d_idx],
                            'ux_r_z0' :ux_r_z0[final_2d_idx,:,:],
                            'ux_r_x0' :ux_r_x0[final_2d_idx,:,:],
                            'ux_r_y0' :ux_r_y0[final_2d_idx,:,:],
                            #
                            'uy_r_z0' :uy_r_z0[final_2d_idx,:,:],
                            'uy_r_x0' :uy_r_x0[final_2d_idx,:,:],
                            'uy_r_y0' :uy_r_y0[final_2d_idx,:,:],
                            #
                            'uz_r_z0' :uz_r_z0[final_2d_idx,:,:],
                            'uz_r_x0' :uz_r_x0[final_2d_idx,:,:],
                            'uz_r_y0' :uz_r_y0[final_2d_idx,:,:],
                          })
savemat(outdir + 'b_r' , {
                            'tt'      :tt[final_2d_idx],
                            'bx_r_z0' :bx_r_z0[final_2d_idx,:,:],
                            'bx_r_x0' :bx_r_x0[final_2d_idx,:,:],
                            'bx_r_y0' :bx_r_y0[final_2d_idx,:,:],
                            #
                            'by_r_z0' :by_r_z0[final_2d_idx,:,:],
                            'by_r_x0' :by_r_x0[final_2d_idx,:,:],
                            'by_r_y0' :by_r_y0[final_2d_idx,:,:],
                            #
                            'bz_r_z0' :bz_r_z0[final_2d_idx,:,:],
                            'bz_r_x0' :bz_r_x0[final_2d_idx,:,:],
                            'bz_r_y0' :bz_r_y0[final_2d_idx,:,:],
                          })
savemat(outdir + 'w_r' , {
                            'tt'      :tt[final_2d_idx],
                            'wx_r_z0' :wx_r_z0[final_2d_idx,:,:],
                            'wx_r_x0' :wx_r_x0[final_2d_idx,:,:],
                            'wx_r_y0' :wx_r_y0[final_2d_idx,:,:],
                            #
                            'wy_r_z0' :wy_r_z0[final_2d_idx,:,:],
                            'wy_r_x0' :wy_r_x0[final_2d_idx,:,:],
                            'wy_r_y0' :wy_r_y0[final_2d_idx,:,:],
                            #
                            'wz_r_z0' :wz_r_z0[final_2d_idx,:,:],
                            'wz_r_x0' :wz_r_x0[final_2d_idx,:,:],
                            'wz_r_y0' :wz_r_y0[final_2d_idx,:,:],
                          })
savemat(outdir + 'j_r' , {
                            'tt'      :tt[final_2d_idx],
                            'jx_r_z0' :jx_r_z0[final_2d_idx,:,:],
                            'jx_r_x0' :jx_r_x0[final_2d_idx,:,:],
                            'jx_r_y0' :jx_r_y0[final_2d_idx,:,:],
                            #
                            'jy_r_z0' :jy_r_z0[final_2d_idx,:,:],
                            'jy_r_x0' :jy_r_x0[final_2d_idx,:,:],
                            'jy_r_y0' :jy_r_y0[final_2d_idx,:,:],
                            #
                            'jz_r_z0' :jz_r_z0[final_2d_idx,:,:],
                            'jz_r_x0' :jz_r_x0[final_2d_idx,:,:],
                            'jz_r_y0' :jz_r_y0[final_2d_idx,:,:],
                          })
