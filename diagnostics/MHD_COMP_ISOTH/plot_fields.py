# -*- coding: utf-8 -*-
from load import *
from fft import *
from plots import *

print('\nplotting fields\n')
outdir = './fig_fields/'

tt_fld = np.transpose(np.loadtxt(input_dir+'out2d/time.dat'+restart_num))[0]
nt_fld = tt_fld.size
if nt_fld == 1 : tt_fld = [tt_fld]

rho_r_z0 = np.transpose(np.fromfile(input_dir+'out2d/rho_r_z0.dat'+restart_num).reshape(nt_fld, nly, nlx)); rho_r_z0 = np.transpose(rho_r_z0, axes=(2, 1, 0))
rho_r_x0 = np.transpose(np.fromfile(input_dir+'out2d/rho_r_x0.dat'+restart_num).reshape(nt_fld, nlz, nly)); rho_r_x0 = np.transpose(rho_r_x0, axes=(2, 1, 0))
rho_r_y0 = np.transpose(np.fromfile(input_dir+'out2d/rho_r_y0.dat'+restart_num).reshape(nt_fld, nlz, nlx)); rho_r_y0 = np.transpose(rho_r_y0, axes=(2, 1, 0))

mx_r_z0  = np.transpose(np.fromfile(input_dir+'out2d/mx_r_z0.dat' +restart_num).reshape(nt_fld, nly, nlx)); mx_r_z0  = np.transpose(mx_r_z0 , axes=(2, 1, 0))
mx_r_x0  = np.transpose(np.fromfile(input_dir+'out2d/mx_r_x0.dat' +restart_num).reshape(nt_fld, nlz, nly)); mx_r_x0  = np.transpose(mx_r_x0 , axes=(2, 1, 0))
mx_r_y0  = np.transpose(np.fromfile(input_dir+'out2d/mx_r_y0.dat' +restart_num).reshape(nt_fld, nlz, nlx)); mx_r_y0  = np.transpose(mx_r_y0 , axes=(2, 1, 0))
                                                                                                                                            
my_r_z0  = np.transpose(np.fromfile(input_dir+'out2d/my_r_z0.dat' +restart_num).reshape(nt_fld, nly, nlx)); my_r_z0  = np.transpose(my_r_z0 , axes=(2, 1, 0))
my_r_x0  = np.transpose(np.fromfile(input_dir+'out2d/my_r_x0.dat' +restart_num).reshape(nt_fld, nlz, nly)); my_r_x0  = np.transpose(my_r_x0 , axes=(2, 1, 0))
my_r_y0  = np.transpose(np.fromfile(input_dir+'out2d/my_r_y0.dat' +restart_num).reshape(nt_fld, nlz, nlx)); my_r_y0  = np.transpose(my_r_y0 , axes=(2, 1, 0))
                                                                                                                                            
mz_r_z0  = np.transpose(np.fromfile(input_dir+'out2d/mz_r_z0.dat' +restart_num).reshape(nt_fld, nly, nlx)); mz_r_z0  = np.transpose(mz_r_z0 , axes=(2, 1, 0))
mz_r_x0  = np.transpose(np.fromfile(input_dir+'out2d/mz_r_x0.dat' +restart_num).reshape(nt_fld, nlz, nly)); mz_r_x0  = np.transpose(mz_r_x0 , axes=(2, 1, 0))
mz_r_y0  = np.transpose(np.fromfile(input_dir+'out2d/mz_r_y0.dat' +restart_num).reshape(nt_fld, nlz, nlx)); mz_r_y0  = np.transpose(mz_r_y0 , axes=(2, 1, 0))
                                                                                                                                            
wx_r_z0  = np.transpose(np.fromfile(input_dir+'out2d/wx_r_z0.dat' +restart_num).reshape(nt_fld, nly, nlx)); wx_r_z0  = np.transpose(wx_r_z0 , axes=(2, 1, 0))
wx_r_x0  = np.transpose(np.fromfile(input_dir+'out2d/wx_r_x0.dat' +restart_num).reshape(nt_fld, nlz, nly)); wx_r_x0  = np.transpose(wx_r_x0 , axes=(2, 1, 0))
wx_r_y0  = np.transpose(np.fromfile(input_dir+'out2d/wx_r_y0.dat' +restart_num).reshape(nt_fld, nlz, nlx)); wx_r_y0  = np.transpose(wx_r_y0 , axes=(2, 1, 0))
                                                                                                                                                             
wy_r_z0  = np.transpose(np.fromfile(input_dir+'out2d/wy_r_z0.dat' +restart_num).reshape(nt_fld, nly, nlx)); wy_r_z0  = np.transpose(wy_r_z0 , axes=(2, 1, 0))
wy_r_x0  = np.transpose(np.fromfile(input_dir+'out2d/wy_r_x0.dat' +restart_num).reshape(nt_fld, nlz, nly)); wy_r_x0  = np.transpose(wy_r_x0 , axes=(2, 1, 0))
wy_r_y0  = np.transpose(np.fromfile(input_dir+'out2d/wy_r_y0.dat' +restart_num).reshape(nt_fld, nlz, nlx)); wy_r_y0  = np.transpose(wy_r_y0 , axes=(2, 1, 0))
                                                                                                                                                             
wz_r_z0  = np.transpose(np.fromfile(input_dir+'out2d/wz_r_z0.dat' +restart_num).reshape(nt_fld, nly, nlx)); wz_r_z0  = np.transpose(wz_r_z0 , axes=(2, 1, 0))
wz_r_x0  = np.transpose(np.fromfile(input_dir+'out2d/wz_r_x0.dat' +restart_num).reshape(nt_fld, nlz, nly)); wz_r_x0  = np.transpose(wz_r_x0 , axes=(2, 1, 0))
wz_r_y0  = np.transpose(np.fromfile(input_dir+'out2d/wz_r_y0.dat' +restart_num).reshape(nt_fld, nlz, nlx)); wz_r_y0  = np.transpose(wz_r_y0 , axes=(2, 1, 0))
                                                                                                                                            
bx_r_z0  = np.transpose(np.fromfile(input_dir+'out2d/bx_r_z0.dat' +restart_num).reshape(nt_fld, nly, nlx)); bx_r_z0  = np.transpose(bx_r_z0 , axes=(2, 1, 0))
bx_r_x0  = np.transpose(np.fromfile(input_dir+'out2d/bx_r_x0.dat' +restart_num).reshape(nt_fld, nlz, nly)); bx_r_x0  = np.transpose(bx_r_x0 , axes=(2, 1, 0))
bx_r_y0  = np.transpose(np.fromfile(input_dir+'out2d/bx_r_y0.dat' +restart_num).reshape(nt_fld, nlz, nlx)); bx_r_y0  = np.transpose(bx_r_y0 , axes=(2, 1, 0))
                                                                                                                                                             
by_r_z0  = np.transpose(np.fromfile(input_dir+'out2d/by_r_z0.dat' +restart_num).reshape(nt_fld, nly, nlx)); by_r_z0  = np.transpose(by_r_z0 , axes=(2, 1, 0))
by_r_x0  = np.transpose(np.fromfile(input_dir+'out2d/by_r_x0.dat' +restart_num).reshape(nt_fld, nlz, nly)); by_r_x0  = np.transpose(by_r_x0 , axes=(2, 1, 0))
by_r_y0  = np.transpose(np.fromfile(input_dir+'out2d/by_r_y0.dat' +restart_num).reshape(nt_fld, nlz, nlx)); by_r_y0  = np.transpose(by_r_y0 , axes=(2, 1, 0))
                                                                                                                                                             
bz_r_z0  = np.transpose(np.fromfile(input_dir+'out2d/bz_r_z0.dat' +restart_num).reshape(nt_fld, nly, nlx)); bz_r_z0  = np.transpose(bz_r_z0 , axes=(2, 1, 0))
bz_r_x0  = np.transpose(np.fromfile(input_dir+'out2d/bz_r_x0.dat' +restart_num).reshape(nt_fld, nlz, nly)); bz_r_x0  = np.transpose(bz_r_x0 , axes=(2, 1, 0))
bz_r_y0  = np.transpose(np.fromfile(input_dir+'out2d/bz_r_y0.dat' +restart_num).reshape(nt_fld, nlz, nlx)); bz_r_y0  = np.transpose(bz_r_y0 , axes=(2, 1, 0))
                                                                                                                                            
jx_r_z0  = np.transpose(np.fromfile(input_dir+'out2d/jx_r_z0.dat' +restart_num).reshape(nt_fld, nly, nlx)); jx_r_z0  = np.transpose(jx_r_z0 , axes=(2, 1, 0))
jx_r_x0  = np.transpose(np.fromfile(input_dir+'out2d/jx_r_x0.dat' +restart_num).reshape(nt_fld, nlz, nly)); jx_r_x0  = np.transpose(jx_r_x0 , axes=(2, 1, 0))
jx_r_y0  = np.transpose(np.fromfile(input_dir+'out2d/jx_r_y0.dat' +restart_num).reshape(nt_fld, nlz, nlx)); jx_r_y0  = np.transpose(jx_r_y0 , axes=(2, 1, 0))
                                                                                                                                                             
jy_r_z0  = np.transpose(np.fromfile(input_dir+'out2d/jy_r_z0.dat' +restart_num).reshape(nt_fld, nly, nlx)); jy_r_z0  = np.transpose(jy_r_z0 , axes=(2, 1, 0))
jy_r_x0  = np.transpose(np.fromfile(input_dir+'out2d/jy_r_x0.dat' +restart_num).reshape(nt_fld, nlz, nly)); jy_r_x0  = np.transpose(jy_r_x0 , axes=(2, 1, 0))
jy_r_y0  = np.transpose(np.fromfile(input_dir+'out2d/jy_r_y0.dat' +restart_num).reshape(nt_fld, nlz, nlx)); jy_r_y0  = np.transpose(jy_r_y0 , axes=(2, 1, 0))
                                                                                                                                                             
jz_r_z0  = np.transpose(np.fromfile(input_dir+'out2d/jz_r_z0.dat' +restart_num).reshape(nt_fld, nly, nlx)); jz_r_z0  = np.transpose(jz_r_z0 , axes=(2, 1, 0))
jz_r_x0  = np.transpose(np.fromfile(input_dir+'out2d/jz_r_x0.dat' +restart_num).reshape(nt_fld, nlz, nly)); jz_r_x0  = np.transpose(jz_r_x0 , axes=(2, 1, 0))
jz_r_y0  = np.transpose(np.fromfile(input_dir+'out2d/jz_r_y0.dat' +restart_num).reshape(nt_fld, nlz, nlx)); jz_r_y0  = np.transpose(jz_r_y0 , axes=(2, 1, 0))

ux_r_z0  = mx_r_z0/rho_r_z0 
ux_r_x0  = mx_r_x0/rho_r_x0 
ux_r_y0  = mx_r_y0/rho_r_y0 
          
uy_r_z0  = my_r_z0/rho_r_z0
uy_r_x0  = my_r_x0/rho_r_x0 
uy_r_y0  = my_r_y0/rho_r_y0 
           
uz_r_z0  = mz_r_z0/rho_r_z0
uz_r_x0  = mz_r_x0/rho_r_x0 
uz_r_y0  = mz_r_y0/rho_r_y0 

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
plot_3d(rho_r_z0[final_fld_idx,:,:], rho_r_y0[final_fld_idx,:,:], rho_r_x0[final_fld_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$\rho (t = %.2E)$'     % tt_fld[final_fld_idx], cmp=parula_map, save=outdir+'rho.pdf')
plot_3d(u_r_z0  [final_fld_idx,:,:], u_r_y0  [final_fld_idx,:,:], u_r_x0  [final_fld_idx,:,:], xx, yy, zz, (ux_r_z0[final_fld_idx,:,:], uy_r_z0[final_fld_idx,:,:]), (ux_r_y0[final_fld_idx,:,:], uz_r_y0[final_fld_idx,:,:]), (uy_r_x0[final_fld_idx,:,:], uz_r_x0[final_fld_idx,:,:]), xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\bm{u}| (t = %.2E)$' % tt_fld[final_fld_idx], cmp=parula_map, streamline_density=1, streamline_width=0.5, streamline_color='w', save=outdir+'u.pdf')
plot_3d(b_r_z0  [final_fld_idx,:,:], b_r_y0  [final_fld_idx,:,:], b_r_x0  [final_fld_idx,:,:], xx, yy, zz, (bx_r_z0[final_fld_idx,:,:], by_r_z0[final_fld_idx,:,:]), (bx_r_y0[final_fld_idx,:,:], bz_r_y0[final_fld_idx,:,:]), (by_r_x0[final_fld_idx,:,:], bz_r_x0[final_fld_idx,:,:]), xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\bm{B}| (t = %.2E)$' % tt_fld[final_fld_idx], cmp=parula_map, streamline_density=1, streamline_width=0.5, streamline_color='w', save=outdir+'b.pdf')
if is2D:
  plot_3d(wz_r_z0[final_fld_idx,:,:], wz_r_y0[final_fld_idx,:,:], wz_r_x0[final_fld_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$(\Curl\bm{u})_z (t = %.2E)$' % tt_fld[final_fld_idx], cmp='RdBu_r', save=outdir+'w.pdf')
  plot_3d(jz_r_z0[final_fld_idx,:,:], jz_r_y0[final_fld_idx,:,:], jz_r_x0[final_fld_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$(\Curl\bm{B})_z (t = %.2E)$' % tt_fld[final_fld_idx], cmp='RdBu_r', save=outdir+'j.pdf')
else:
  plot_3d( w_r_z0[final_fld_idx,:,:],  w_r_y0[final_fld_idx,:,:],  w_r_x0[final_fld_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\Curl\bm{u}| (t = %.2E)$' % tt_fld[final_fld_idx], cmp=parula_map, save=outdir+'w.pdf')
  plot_3d( j_r_z0[final_fld_idx,:,:],  j_r_y0[final_fld_idx,:,:],  j_r_x0[final_fld_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\Curl\bm{B}| (t = %.2E)$' % tt_fld[final_fld_idx], cmp=parula_map, save=outdir+'j.pdf')


#--------------------------------------------------------#
#                       plot movie                       #
#--------------------------------------------------------#
if ismovie:
  # evenly space time
  # nframe = min(200, int(tt[:final_fld_idx].size))
  nframe = tt_fld[:final_fld_idx].size
  idx = np.unique([np.argmin(abs(tt_fld - np.linspace(tt_fld[0], tt_fld[-1], nframe)[i])) for i in range(0, nframe)])
  tt_fld  = tt_fld [idx]

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
  movie_3d(tt_fld, rho_r_z0, rho_r_y0, rho_r_x0, xx_fld, yy_fld, zz_fld, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$\rho$'    , cmp=parula_map, save=outdir+'rho_anim.gif')
  movie_3d(tt_fld, u_r_z0  , u_r_y0, u_r_x0    , xx_fld, yy_fld, zz_fld, (ux_r_z0, uy_r_z0), (ux_r_y0, uz_r_y0), (uy_r_x0, uz_r_x0), xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\bm{u}|$', cmp=parula_map, save=outdir+'u_anim.gif')
  movie_3d(tt_fld, b_r_z0  , b_r_y0, b_r_x0    , xx_fld, yy_fld, zz_fld, (bx_r_z0, by_r_z0), (bx_r_y0, bz_r_y0), (by_r_x0, bz_r_x0), xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\bm{B}|$', cmp=parula_map, save=outdir+'b_anim.gif')
  if is2D:
    movie_3d(tt_fld, wz_r_z0, wz_r_y0, wz_r_x0, xx_fld, yy_fld, zz_fld, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$(\Curl\bm{u})_z$', cmp='RdBu_r', save=outdir+'w_anim.gif')
    movie_3d(tt_fld, jz_r_z0, jz_r_y0, jz_r_x0, xx_fld, yy_fld, zz_fld, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$(\Curl\bm{B})_z$', cmp='RdBu_r', save=outdir+'j_anim.gif')
  else:
    movie_3d(tt_fld,  w_r_z0,  w_r_y0,  w_r_x0, xx_fld, yy_fld, zz_fld, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\Curl\bm{u}|$', cmp=parula_map, save=outdir+'w_anim.gif')
    movie_3d(tt_fld,  j_r_z0,  j_r_y0,  j_r_x0, xx_fld, yy_fld, zz_fld, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\Curl\bm{B}|$', cmp=parula_map, save=outdir+'j_anim.gif')


#------------------#
#   output data    #
#------------------#
from scipy.io import savemat

savemat(outdir + 'grid' , {'xx':xx, 'yy':yy, 'zz':zz})
savemat(outdir + 'rho_r' , {
                            'tt'       :tt[final_fld_idx],
                            'rho_r_z0' :rho_r_z0[final_fld_idx,:,:],
                            'rho_r_x0' :rho_r_x0[final_fld_idx,:,:],
                            'rho_r_y0' :rho_r_y0[final_fld_idx,:,:],
                          })
savemat(outdir + 'u_r' , {
                            'tt'       :tt[final_fld_idx],
                            'ux_r_z0'  :ux_r_z0 [final_fld_idx,:,:],
                            'ux_r_x0'  :ux_r_x0 [final_fld_idx,:,:],
                            'ux_r_y0'  :ux_r_y0 [final_fld_idx,:,:],
                            #
                            'uy_r_z0'  :uy_r_z0 [final_fld_idx,:,:],
                            'uy_r_x0'  :uy_r_x0 [final_fld_idx,:,:],
                            'uy_r_y0'  :uy_r_y0 [final_fld_idx,:,:],
                            #
                            'uz_r_z0'  :uz_r_z0 [final_fld_idx,:,:],
                            'uz_r_x0'  :uz_r_x0 [final_fld_idx,:,:],
                            'uz_r_y0'  :uz_r_y0 [final_fld_idx,:,:],
                          })
savemat(outdir + 'b_r' , {
                            'tt'       :tt[final_fld_idx],
                            'bx_r_z0'  :bx_r_z0 [final_fld_idx,:,:],
                            'bx_r_x0'  :bx_r_x0 [final_fld_idx,:,:],
                            'bx_r_y0'  :bx_r_y0 [final_fld_idx,:,:],
                            #
                            'by_r_z0'  :by_r_z0 [final_fld_idx,:,:],
                            'by_r_x0'  :by_r_x0 [final_fld_idx,:,:],
                            'by_r_y0'  :by_r_y0 [final_fld_idx,:,:],
                            #
                            'bz_r_z0'  :bz_r_z0 [final_fld_idx,:,:],
                            'bz_r_x0'  :bz_r_x0 [final_fld_idx,:,:],
                            'bz_r_y0'  :bz_r_y0 [final_fld_idx,:,:],
                          })
savemat(outdir + 'w_r' , {
                            'tt'       :tt[final_fld_idx],
                            'wx_r_z0'  :wx_r_z0 [final_fld_idx,:,:],
                            'wx_r_x0'  :wx_r_x0 [final_fld_idx,:,:],
                            'wx_r_y0'  :wx_r_y0 [final_fld_idx,:,:],
                            #
                            'wy_r_z0'  :wy_r_z0 [final_fld_idx,:,:],
                            'wy_r_x0'  :wy_r_x0 [final_fld_idx,:,:],
                            'wy_r_y0'  :wy_r_y0 [final_fld_idx,:,:],
                            #
                            'wz_r_z0'  :wz_r_z0 [final_fld_idx,:,:],
                            'wz_r_x0'  :wz_r_x0 [final_fld_idx,:,:],
                            'wz_r_y0'  :wz_r_y0 [final_fld_idx,:,:],
                          })
savemat(outdir + 'j_r' , {
                            'tt'       :tt[final_fld_idx],
                            'jx_r_z0'  :jx_r_z0 [final_fld_idx,:,:],
                            'jx_r_x0'  :jx_r_x0 [final_fld_idx,:,:],
                            'jx_r_y0'  :jx_r_y0 [final_fld_idx,:,:],
                            #
                            'jy_r_z0'  :jy_r_z0 [final_fld_idx,:,:],
                            'jy_r_x0'  :jy_r_x0 [final_fld_idx,:,:],
                            'jy_r_y0'  :jy_r_y0 [final_fld_idx,:,:],
                            #
                            'jz_r_z0'  :jz_r_z0 [final_fld_idx,:,:],
                            'jz_r_x0'  :jz_r_x0 [final_fld_idx,:,:],
                            'jz_r_y0'  :jz_r_y0 [final_fld_idx,:,:],
                          })
