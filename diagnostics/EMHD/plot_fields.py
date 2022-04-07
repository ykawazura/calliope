# -*- coding: utf-8 -*-
from load import *
from fft import *
from plots import *

print('\nplotting fields\n')
outdir = './fig_fields/'

# Load 2D cut
ncfile = netcdf.netcdf_file(input_dir+runname+'.out.2D.nc'+restart_num, 'r')   

tt_fld  = np.copy(ncfile.variables['tt' ][:]); tt_fld  = np.delete(tt_fld , ignored_points_fld, axis = 0)
xx_fld  = np.copy(ncfile.variables['xx' ][:])
yy_fld  = np.copy(ncfile.variables['yy' ][:])
zz_fld  = np.copy(ncfile.variables['zz' ][:])

bx_r_z0 = np.copy(ncfile.variables['bx_r_z0'][:]);  bx_r_z0 = np.delete(bx_r_z0, ignored_points_fld, axis = 0)
bx_r_x0 = np.copy(ncfile.variables['bx_r_x0'][:]);  bx_r_x0 = np.delete(bx_r_x0, ignored_points_fld, axis = 0)
bx_r_y0 = np.copy(ncfile.variables['bx_r_y0'][:]);  bx_r_y0 = np.delete(bx_r_y0, ignored_points_fld, axis = 0)

by_r_z0 = np.copy(ncfile.variables['by_r_z0'][:]);  by_r_z0 = np.delete(by_r_z0, ignored_points_fld, axis = 0)
by_r_x0 = np.copy(ncfile.variables['by_r_x0'][:]);  by_r_x0 = np.delete(by_r_x0, ignored_points_fld, axis = 0)
by_r_y0 = np.copy(ncfile.variables['by_r_y0'][:]);  by_r_y0 = np.delete(by_r_y0, ignored_points_fld, axis = 0)

bz_r_z0 = np.copy(ncfile.variables['bz_r_z0'][:]);  bz_r_z0 = np.delete(bz_r_z0, ignored_points_fld, axis = 0)
bz_r_x0 = np.copy(ncfile.variables['bz_r_x0'][:]);  bz_r_x0 = np.delete(bz_r_x0, ignored_points_fld, axis = 0)
bz_r_y0 = np.copy(ncfile.variables['bz_r_y0'][:]);  bz_r_y0 = np.delete(bz_r_y0, ignored_points_fld, axis = 0)

jx_r_z0 = np.copy(ncfile.variables['jx_r_z0'][:]);  jx_r_z0 = np.delete(jx_r_z0, ignored_points_fld, axis = 0)
jx_r_x0 = np.copy(ncfile.variables['jx_r_x0'][:]);  jx_r_x0 = np.delete(jx_r_x0, ignored_points_fld, axis = 0)
jx_r_y0 = np.copy(ncfile.variables['jx_r_y0'][:]);  jx_r_y0 = np.delete(jx_r_y0, ignored_points_fld, axis = 0)

jy_r_z0 = np.copy(ncfile.variables['jy_r_z0'][:]);  jy_r_z0 = np.delete(jy_r_z0, ignored_points_fld, axis = 0)
jy_r_x0 = np.copy(ncfile.variables['jy_r_x0'][:]);  jy_r_x0 = np.delete(jy_r_x0, ignored_points_fld, axis = 0)
jy_r_y0 = np.copy(ncfile.variables['jy_r_y0'][:]);  jy_r_y0 = np.delete(jy_r_y0, ignored_points_fld, axis = 0)

jz_r_z0 = np.copy(ncfile.variables['jz_r_z0'][:]);  jz_r_z0 = np.delete(jz_r_z0, ignored_points_fld, axis = 0)
jz_r_x0 = np.copy(ncfile.variables['jz_r_x0'][:]);  jz_r_x0 = np.delete(jz_r_x0, ignored_points_fld, axis = 0)
jz_r_y0 = np.copy(ncfile.variables['jz_r_y0'][:]);  jz_r_y0 = np.delete(jz_r_y0, ignored_points_fld, axis = 0)

ncfile.close()

#--------------------------------------------------------#
#                   plot final snapshot                  #
#--------------------------------------------------------#
b_r_z0 = np.sqrt(bx_r_z0**2 + by_r_z0**2 + bz_r_z0**2)
b_r_y0 = np.sqrt(bx_r_y0**2 + by_r_y0**2 + bz_r_y0**2)
b_r_x0 = np.sqrt(bx_r_x0**2 + by_r_x0**2 + bz_r_x0**2)
j_r_z0 = np.sqrt(jx_r_z0**2 + jy_r_z0**2 + jz_r_z0**2)
j_r_y0 = np.sqrt(jx_r_y0**2 + jy_r_y0**2 + jz_r_y0**2)
j_r_x0 = np.sqrt(jx_r_x0**2 + jy_r_x0**2 + jz_r_x0**2)
plot_3d(b_r_z0[final_fld_idx,:,:], b_r_y0[final_fld_idx,:,:], b_r_x0[final_fld_idx,:,:], xx_fld, yy_fld, zz_fld, (bx_r_z0[final_fld_idx,:,:], by_r_z0[final_fld_idx,:,:]), (bx_r_y0[final_fld_idx,:,:], bz_r_y0[final_fld_idx,:,:]), (by_r_x0[final_fld_idx,:,:], bz_r_x0[final_fld_idx,:,:]), xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\bm{B}| (t = %.2E)$' % tt_fld[final_fld_idx], cmp=parula_map, streamline_density=1, streamline_width=0.5, streamline_color='w', save=outdir+'b.pdf')
if is2D:
  plot_3d(jz_r_z0[final_fld_idx,:,:], jz_r_y0[final_fld_idx,:,:], jz_r_x0[final_fld_idx,:,:], xx_fld, yy_fld, zz_fld, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$(\Curl\bm{B})_z (t = %.2E)$' % tt_fld[final_fld_idx], cmp='RdBu_r', save=outdir+'j.pdf')
else:
  plot_3d( j_r_z0[final_fld_idx,:,:],  j_r_y0[final_fld_idx,:,:],  j_r_x0[final_fld_idx,:,:], xx_fld, yy_fld, zz_fld, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\Curl\bm{B}| (t = %.2E)$' % tt_fld[final_fld_idx], cmp=parula_map, save=outdir+'j.pdf')


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
  movie_3d(tt_fld, b_r_z0  , b_r_y0, b_r_x0    , xx_fld, yy_fld, zz_fld, (bx_r_z0, by_r_z0), (bx_r_y0, bz_r_y0), (by_r_x0, bz_r_x0), xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\bm{B}|$', cmp=parula_map, save=outdir+'b_anim.gif')
  if is2D:
    movie_3d(tt_fld, jz_r_z0, jz_r_y0, jz_r_x0, xx_fld, yy_fld, zz_fld, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$(\Curl\bm{B})_z$', cmp='RdBu_r', save=outdir+'j_anim.gif')
  else:
    movie_3d(tt_fld,  j_r_z0,  j_r_y0,  j_r_x0, xx_fld, yy_fld, zz_fld, xlab=xlab, ylab=ylab, zlab=zlab, title=r'$|\Curl\bm{B}|$', cmp=parula_map, save=outdir+'j_anim.gif')


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
