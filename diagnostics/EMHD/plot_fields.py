# -*- coding: utf-8 -*-
from load import *
from fft import *
import sys
sys.path.append('../')
from plots import *

print('\nplotting fields\n')
outdir = './fig_fields/'

tt_fld = np.loadtxt(input_dir+'out2d/time.dat'+restart_num)
nt_fld = tt_fld.size
if nt_fld == 1 : tt_fld = [tt_fld]

bx_r_z0 = np.transpose(np.fromfile(input_dir+'out2d/bx_r_z0.dat'+restart_num).reshape(nt_fld, nly, nlx)); bx_r_z0 = np.transpose(bx_r_z0, axes=(2, 1, 0))
bx_r_x0 = np.transpose(np.fromfile(input_dir+'out2d/bx_r_x0.dat'+restart_num).reshape(nt_fld, nlz, nly)); bx_r_x0 = np.transpose(bx_r_x0, axes=(2, 1, 0))
bx_r_y0 = np.transpose(np.fromfile(input_dir+'out2d/bx_r_y0.dat'+restart_num).reshape(nt_fld, nlz, nlx)); bx_r_y0 = np.transpose(bx_r_y0, axes=(2, 1, 0))
                                                                                                                                                         
by_r_z0 = np.transpose(np.fromfile(input_dir+'out2d/by_r_z0.dat'+restart_num).reshape(nt_fld, nly, nlx)); by_r_z0 = np.transpose(by_r_z0, axes=(2, 1, 0))
by_r_x0 = np.transpose(np.fromfile(input_dir+'out2d/by_r_x0.dat'+restart_num).reshape(nt_fld, nlz, nly)); by_r_x0 = np.transpose(by_r_x0, axes=(2, 1, 0))
by_r_y0 = np.transpose(np.fromfile(input_dir+'out2d/by_r_y0.dat'+restart_num).reshape(nt_fld, nlz, nlx)); by_r_y0 = np.transpose(by_r_y0, axes=(2, 1, 0))
                                                                                                                                                         
bz_r_z0 = np.transpose(np.fromfile(input_dir+'out2d/bz_r_z0.dat'+restart_num).reshape(nt_fld, nly, nlx)); bz_r_z0 = np.transpose(bz_r_z0, axes=(2, 1, 0))
bz_r_x0 = np.transpose(np.fromfile(input_dir+'out2d/bz_r_x0.dat'+restart_num).reshape(nt_fld, nlz, nly)); bz_r_x0 = np.transpose(bz_r_x0, axes=(2, 1, 0))
bz_r_y0 = np.transpose(np.fromfile(input_dir+'out2d/bz_r_y0.dat'+restart_num).reshape(nt_fld, nlz, nlx)); bz_r_y0 = np.transpose(bz_r_y0, axes=(2, 1, 0))

jx_r_z0 = np.transpose(np.fromfile(input_dir+'out2d/jx_r_z0.dat'+restart_num).reshape(nt_fld, nly, nlx)); jx_r_z0 = np.transpose(jx_r_z0, axes=(2, 1, 0))
jx_r_x0 = np.transpose(np.fromfile(input_dir+'out2d/jx_r_x0.dat'+restart_num).reshape(nt_fld, nlz, nly)); jx_r_x0 = np.transpose(jx_r_x0, axes=(2, 1, 0))
jx_r_y0 = np.transpose(np.fromfile(input_dir+'out2d/jx_r_y0.dat'+restart_num).reshape(nt_fld, nlz, nlx)); jx_r_y0 = np.transpose(jx_r_y0, axes=(2, 1, 0))
                                                                                                                                                         
jy_r_z0 = np.transpose(np.fromfile(input_dir+'out2d/jy_r_z0.dat'+restart_num).reshape(nt_fld, nly, nlx)); jy_r_z0 = np.transpose(jy_r_z0, axes=(2, 1, 0))
jy_r_x0 = np.transpose(np.fromfile(input_dir+'out2d/jy_r_x0.dat'+restart_num).reshape(nt_fld, nlz, nly)); jy_r_x0 = np.transpose(jy_r_x0, axes=(2, 1, 0))
jy_r_y0 = np.transpose(np.fromfile(input_dir+'out2d/jy_r_y0.dat'+restart_num).reshape(nt_fld, nlz, nlx)); jy_r_y0 = np.transpose(jy_r_y0, axes=(2, 1, 0))
                                                                                                                                                         
jz_r_z0 = np.transpose(np.fromfile(input_dir+'out2d/jz_r_z0.dat'+restart_num).reshape(nt_fld, nly, nlx)); jz_r_z0 = np.transpose(jz_r_z0, axes=(2, 1, 0))
jz_r_x0 = np.transpose(np.fromfile(input_dir+'out2d/jz_r_x0.dat'+restart_num).reshape(nt_fld, nlz, nly)); jz_r_x0 = np.transpose(jz_r_x0, axes=(2, 1, 0))
jz_r_y0 = np.transpose(np.fromfile(input_dir+'out2d/jz_r_y0.dat'+restart_num).reshape(nt_fld, nlz, nlx)); jz_r_y0 = np.transpose(jz_r_y0, axes=(2, 1, 0))

#--------------------------------------------------------#
#                   plot final snapshot                  #
#--------------------------------------------------------#
b_r_z0 = np.sqrt(bx_r_z0**2 + by_r_z0**2 + bz_r_z0**2)
b_r_y0 = np.sqrt(bx_r_y0**2 + by_r_y0**2 + bz_r_y0**2)
b_r_x0 = np.sqrt(bx_r_x0**2 + by_r_x0**2 + bz_r_x0**2)
j_r_z0 = np.sqrt(jx_r_z0**2 + jy_r_z0**2 + jz_r_z0**2)
j_r_y0 = np.sqrt(jx_r_y0**2 + jy_r_y0**2 + jz_r_y0**2)
j_r_x0 = np.sqrt(jx_r_x0**2 + jy_r_x0**2 + jz_r_x0**2)
plot_3d(b_r_z0[final_fld_idx,:,:], b_r_y0[final_fld_idx,:,:], b_r_x0[final_fld_idx,:,:], xx, yy, zz, (bx_r_z0[final_fld_idx,:,:], by_r_z0[final_fld_idx,:,:]), (bx_r_y0[final_fld_idx,:,:], bz_r_y0[final_fld_idx,:,:]), (by_r_x0[final_fld_idx,:,:], bz_r_x0[final_fld_idx,:,:]), xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\bm{B}| (t = %.2E)$' % tt_fld[final_fld_idx], cmp='sequential', streamline_density=1, streamline_width=0.5, streamline_color='w', save=outdir+'b.pdf')
if is2D:
  plot_3d(jz_r_z0[final_fld_idx,:,:], jz_r_y0[final_fld_idx,:,:], jz_r_x0[final_fld_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$(\Curl\bm{B})_z (t = %.2E)$' % tt_fld[final_fld_idx], cmp='diverging', save=outdir+'j.pdf')
else:
  plot_3d( j_r_z0[final_fld_idx,:,:],  j_r_y0[final_fld_idx,:,:],  j_r_x0[final_fld_idx,:,:], xx, yy, zz, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\Curl\bm{B}| (t = %.2E)$' % tt_fld[final_fld_idx], cmp='sequential', save=outdir+'j.pdf')


#--------------------------------------------------------#
#                       plot movie                       #
#--------------------------------------------------------#
if ismovie:
  # evenly space time
  nframe = min(200, int(tt_fld[:final_fld_idx].size))
  idx = np.unique([np.argmin(abs(tt_fld - np.linspace(tt_fld[0], tt_fld[-1], nframe)[i])) for i in range(0, nframe)])
  tt_fld  = tt_fld [idx]

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

  b_r_z0 = np.sqrt(bx_r_z0**2 + by_r_z0**2 + bz_r_z0**2)
  b_r_y0 = np.sqrt(bx_r_y0**2 + by_r_y0**2 + bz_r_y0**2)
  b_r_x0 = np.sqrt(bx_r_x0**2 + by_r_x0**2 + bz_r_x0**2)
  j_r_z0 = np.sqrt(jx_r_z0**2 + jy_r_z0**2 + jz_r_z0**2)
  j_r_y0 = np.sqrt(jx_r_y0**2 + jy_r_y0**2 + jz_r_y0**2)
  j_r_x0 = np.sqrt(jx_r_x0**2 + jy_r_x0**2 + jz_r_x0**2)
  movie_3d(tt_fld, b_r_z0  , b_r_y0, b_r_x0    , xx_fld, yy_fld, zz_fld, (bx_r_z0, by_r_z0), (bx_r_y0, bz_r_y0), (by_r_x0, bz_r_x0), xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\bm{B}|$', cmp='sequential', save=outdir+'b_anim.gif')
  if is2D:
    movie_3d(tt_fld, jz_r_z0, jz_r_y0, jz_r_x0, xx_fld, yy_fld, zz_fld, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$(\Curl\bm{B})_z$', cmp='diverging', save=outdir+'j_anim.gif')
  else:
    movie_3d(tt_fld,  j_r_z0,  j_r_y0,  j_r_x0, xx_fld, yy_fld, zz_fld, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\Curl\bm{B}|$', cmp='sequential', save=outdir+'j_anim.gif')


#------------------#
#   output data    #
#------------------#
from scipy.io import savemat

savemat(outdir + 'grid' , {'xx':xx, 'yy':yy, 'zz':zz})
savemat(outdir + 'b_r' , {
                            'tt'      :tt[final_fld_idx],
                            'bx_r_z0' :bx_r_z0[final_fld_idx,:,:],
                            'bx_r_x0' :bx_r_x0[final_fld_idx,:,:],
                            'bx_r_y0' :bx_r_y0[final_fld_idx,:,:],
                            #
                            'by_r_z0' :by_r_z0[final_fld_idx,:,:],
                            'by_r_x0' :by_r_x0[final_fld_idx,:,:],
                            'by_r_y0' :by_r_y0[final_fld_idx,:,:],
                            #
                            'bz_r_z0' :bz_r_z0[final_fld_idx,:,:],
                            'bz_r_x0' :bz_r_x0[final_fld_idx,:,:],
                            'bz_r_y0' :bz_r_y0[final_fld_idx,:,:],
                          })
savemat(outdir + 'j_r' , {
                            'tt'      :tt[final_fld_idx],
                            'jx_r_z0' :jx_r_z0[final_fld_idx,:,:],
                            'jx_r_x0' :jx_r_x0[final_fld_idx,:,:],
                            'jx_r_y0' :jx_r_y0[final_fld_idx,:,:],
                            #
                            'jy_r_z0' :jy_r_z0[final_fld_idx,:,:],
                            'jy_r_x0' :jy_r_x0[final_fld_idx,:,:],
                            'jy_r_y0' :jy_r_y0[final_fld_idx,:,:],
                            #
                            'jz_r_z0' :jz_r_z0[final_fld_idx,:,:],
                            'jz_r_x0' :jz_r_x0[final_fld_idx,:,:],
                            'jz_r_y0' :jz_r_y0[final_fld_idx,:,:],
                          })
