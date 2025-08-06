# -*- coding: utf-8 -*-
from load import *
from fft import *
import sys
sys.path.append('../')
from plots import *

print('\nplotting fields\n')
outdir = './fig_fields/'

tt_fld = np.loadtxt(input_dir+'out2d'+restart_num+'/time.dat')
nt_fld = tt_fld.size
if nt_fld == 1 : tt_fld = [tt_fld]

phi_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/phi_r_z0.dat').reshape(nt_fld, nly, nlx)); phi_r_z0 = np.transpose(phi_r_z0, axes=(2, 1, 0))
phi_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/phi_r_x0.dat').reshape(nt_fld, nlz, nly)); phi_r_x0 = np.transpose(phi_r_x0, axes=(2, 1, 0))
phi_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/phi_r_y0.dat').reshape(nt_fld, nlz, nlx)); phi_r_y0 = np.transpose(phi_r_y0, axes=(2, 1, 0))
                                                                   
psi_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/psi_r_z0.dat').reshape(nt_fld, nly, nlx)); psi_r_z0 = np.transpose(psi_r_z0, axes=(2, 1, 0))
psi_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/psi_r_x0.dat').reshape(nt_fld, nlz, nly)); psi_r_x0 = np.transpose(psi_r_x0, axes=(2, 1, 0))
psi_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/psi_r_y0.dat').reshape(nt_fld, nlz, nlx)); psi_r_y0 = np.transpose(psi_r_y0, axes=(2, 1, 0))
                                                                   
omg_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/omg_r_z0.dat').reshape(nt_fld, nly, nlx)); omg_r_z0 = np.transpose(omg_r_z0, axes=(2, 1, 0))
omg_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/omg_r_x0.dat').reshape(nt_fld, nlz, nly)); omg_r_x0 = np.transpose(omg_r_x0, axes=(2, 1, 0))
omg_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/omg_r_y0.dat').reshape(nt_fld, nlz, nlx)); omg_r_y0 = np.transpose(omg_r_y0, axes=(2, 1, 0))
                                                                   
jpa_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/jpa_r_z0.dat').reshape(nt_fld, nly, nlx)); jpa_r_z0 = np.transpose(jpa_r_z0, axes=(2, 1, 0))
jpa_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/jpa_r_x0.dat').reshape(nt_fld, nlz, nly)); jpa_r_x0 = np.transpose(jpa_r_x0, axes=(2, 1, 0))
jpa_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/jpa_r_y0.dat').reshape(nt_fld, nlz, nlx)); jpa_r_y0 = np.transpose(jpa_r_y0, axes=(2, 1, 0))
                                                                   
upa_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/upa_r_z0.dat').reshape(nt_fld, nly, nlx)); upa_r_z0 = np.transpose(upa_r_z0, axes=(2, 1, 0))
upa_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/upa_r_x0.dat').reshape(nt_fld, nlz, nly)); upa_r_x0 = np.transpose(upa_r_x0, axes=(2, 1, 0))
upa_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/upa_r_y0.dat').reshape(nt_fld, nlz, nlx)); upa_r_y0 = np.transpose(upa_r_y0, axes=(2, 1, 0))
                                                                   
bpa_r_z0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/bpa_r_z0.dat').reshape(nt_fld, nly, nlx)); bpa_r_z0 = np.transpose(bpa_r_z0, axes=(2, 1, 0))
bpa_r_x0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/bpa_r_x0.dat').reshape(nt_fld, nlz, nly)); bpa_r_x0 = np.transpose(bpa_r_x0, axes=(2, 1, 0))
bpa_r_y0 = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/bpa_r_y0.dat').reshape(nt_fld, nlz, nlx)); bpa_r_y0 = np.transpose(bpa_r_y0, axes=(2, 1, 0))
                                                                   
ux_r_z0  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/ux_r_z0.dat' ).reshape(nt_fld, nly, nlx)); ux_r_z0  = np.transpose(ux_r_z0 , axes=(2, 1, 0))
ux_r_x0  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/ux_r_x0.dat' ).reshape(nt_fld, nlz, nly)); ux_r_x0  = np.transpose(ux_r_x0 , axes=(2, 1, 0))
ux_r_y0  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/ux_r_y0.dat' ).reshape(nt_fld, nlz, nlx)); ux_r_y0  = np.transpose(ux_r_y0 , axes=(2, 1, 0))
                                                                                                                                               
uy_r_z0  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/uy_r_z0.dat' ).reshape(nt_fld, nly, nlx)); uy_r_z0  = np.transpose(uy_r_z0 , axes=(2, 1, 0))
uy_r_x0  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/uy_r_x0.dat' ).reshape(nt_fld, nlz, nly)); uy_r_x0  = np.transpose(uy_r_x0 , axes=(2, 1, 0))
uy_r_y0  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/uy_r_y0.dat' ).reshape(nt_fld, nlz, nlx)); uy_r_y0  = np.transpose(uy_r_y0 , axes=(2, 1, 0))
                                                                                                                                               
bx_r_z0  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/bx_r_z0.dat' ).reshape(nt_fld, nly, nlx)); bx_r_z0  = np.transpose(bx_r_z0 , axes=(2, 1, 0))
bx_r_x0  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/bx_r_x0.dat' ).reshape(nt_fld, nlz, nly)); bx_r_x0  = np.transpose(bx_r_x0 , axes=(2, 1, 0))
bx_r_y0  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/bx_r_y0.dat' ).reshape(nt_fld, nlz, nlx)); bx_r_y0  = np.transpose(bx_r_y0 , axes=(2, 1, 0))
                                                                                                                                                                
by_r_z0  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/by_r_z0.dat' ).reshape(nt_fld, nly, nlx)); by_r_z0  = np.transpose(by_r_z0 , axes=(2, 1, 0))
by_r_x0  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/by_r_x0.dat' ).reshape(nt_fld, nlz, nly)); by_r_x0  = np.transpose(by_r_x0 , axes=(2, 1, 0))
by_r_y0  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/by_r_y0.dat' ).reshape(nt_fld, nlz, nlx)); by_r_y0  = np.transpose(by_r_y0 , axes=(2, 1, 0))

#--------------------------------------------------------#
#                   plot final snapshot                  #
#--------------------------------------------------------#
upe_r_z0 = np.sqrt(ux_r_z0**2 + uy_r_z0**2)
upe_r_y0 = np.sqrt(ux_r_y0**2 + uy_r_y0**2)
upe_r_x0 = np.sqrt(ux_r_x0**2 + uy_r_x0**2)
bpe_r_z0 = np.sqrt(bx_r_z0**2 + by_r_z0**2)
bpe_r_y0 = np.sqrt(bx_r_y0**2 + by_r_y0**2)
bpe_r_x0 = np.sqrt(bx_r_x0**2 + by_r_x0**2)
plot_3d(phi_r_z0[final_fld_idx,:,:], phi_r_y0[final_fld_idx,:,:], phi_r_x0[final_fld_idx,:,:], xx, yy, zz, xlab=r'$x/\rho_\mathrm{H}$', ylab=r'$y/\rho_\mathrm{H}$', zlab='$'+zlab+'$', title=r'$\Phi (t = %.2E)$' % tt_fld[final_fld_idx], cmp='diverging', save=outdir+'phi.pdf')
plot_3d(omg_r_z0[final_fld_idx,:,:], omg_r_y0[final_fld_idx,:,:], omg_r_x0[final_fld_idx,:,:], xx, yy, zz, xlab=r'$x/\rho_\mathrm{H}$', ylab=r'$y/\rho_\mathrm{H}$', zlab='$'+zlab+'$', title=r'$\nbl_\+^2\Phi (t = %.2E)$' % tt_fld[final_fld_idx], cmp='diverging', save=outdir+'omg.pdf')
plot_3d(psi_r_z0[final_fld_idx,:,:], psi_r_y0[final_fld_idx,:,:], psi_r_x0[final_fld_idx,:,:], xx, yy, zz, xlab=r'$x/\rho_\mathrm{H}$', ylab=r'$y/\rho_\mathrm{H}$', zlab='$'+zlab+'$', title=r'$\Psi (t = %.2E)$' % tt_fld[final_fld_idx], cmp='diverging', save=outdir+'psi.pdf')
plot_3d(jpa_r_z0[final_fld_idx,:,:], jpa_r_y0[final_fld_idx,:,:], jpa_r_x0[final_fld_idx,:,:], xx, yy, zz, xlab=r'$x/\rho_\mathrm{H}$', ylab=r'$y/\rho_\mathrm{H}$', zlab='$'+zlab+'$', title=r'$\nbl_\+^2\Psi (t = %.2E)$' % tt_fld[final_fld_idx], cmp='diverging', save=outdir+'jpa.pdf')
plot_3d(upa_r_z0[final_fld_idx,:,:], upa_r_y0[final_fld_idx,:,:], upa_r_x0[final_fld_idx,:,:], xx, yy, zz, xlab=r'$x/\rho_\mathrm{H}$', ylab=r'$y/\rho_\mathrm{H}$', zlab='$'+zlab+'$', title=r'$u_\| (t = %.2E)$' % tt_fld[final_fld_idx], cmp='diverging', save=outdir+'upa.pdf')
plot_3d(bpa_r_z0[final_fld_idx,:,:], bpa_r_y0[final_fld_idx,:,:], bpa_r_x0[final_fld_idx,:,:], xx, yy, zz, xlab=r'$x/\rho_\mathrm{H}$', ylab=r'$y/\rho_\mathrm{H}$', zlab='$'+zlab+'$', title=r'$\delta B_\| (t = %.2E)$' % tt_fld[final_fld_idx], cmp='diverging', save=outdir+'bpa.pdf')
plot_3d(upe_r_z0[final_fld_idx,:,:], upe_r_y0[final_fld_idx,:,:], upe_r_x0[final_fld_idx,:,:], xx, yy, zz, xlab=r'$x/\rho_\mathrm{H}$', ylab=r'$y/\rho_\mathrm{H}$', zlab='$'+zlab+'$', title=r'$|\bm{u}_\+| (t = %.2E)$' % tt_fld[final_fld_idx], cmp='sequential', save=outdir+'upe.pdf')
plot_3d(bpe_r_z0[final_fld_idx,:,:], bpe_r_y0[final_fld_idx,:,:], bpe_r_x0[final_fld_idx,:,:], xx, yy, zz, xlab=r'$x/\rho_\mathrm{H}$', ylab=r'$y/\rho_\mathrm{H}$', zlab='$'+zlab+'$', title=r'$|\delta\bm{B}_\+| (t = %.2E)$' % tt_fld[final_fld_idx], cmp='sequential', save=outdir+'bpe.pdf')


#--------------------------------------------------------#
#                       plot movie                       #
#--------------------------------------------------------#
if ismovie:
	# evenly space time
  nframe = min(100, int(tt_fld[:final_fld_idx].size))
  idx = np.unique([np.argmin(abs(tt_fld - np.linspace(tt_fld[0], tt_fld[-1], nframe)[i])) for i in range(0, nframe)])
  tt_fld  = tt_fld [idx]

  phi_r_z0 = np.take(phi_r_z0, idx, axis=0)
  phi_r_x0 = np.take(phi_r_x0, idx, axis=0)
  omg_r_z0 = np.take(omg_r_z0, idx, axis=0)
  omg_r_x0 = np.take(omg_r_x0, idx, axis=0)
  psi_r_z0 = np.take(psi_r_z0, idx, axis=0)
  psi_r_x0 = np.take(psi_r_x0, idx, axis=0)
  jpa_r_z0 = np.take(jpa_r_z0, idx, axis=0)
  jpa_r_x0 = np.take(jpa_r_x0, idx, axis=0)
  ux_r_z0  = np.take( ux_r_z0, idx, axis=0)
  ux_r_y0  = np.take( ux_r_y0, idx, axis=0)
  ux_r_x0  = np.take( ux_r_x0, idx, axis=0)
  uy_r_z0  = np.take( uy_r_z0, idx, axis=0)
  uy_r_y0  = np.take( uy_r_y0, idx, axis=0)
  uy_r_x0  = np.take( uy_r_x0, idx, axis=0)
  bx_r_z0  = np.take( bx_r_z0, idx, axis=0)
  bx_r_y0  = np.take( bx_r_y0, idx, axis=0)
  bx_r_x0  = np.take( bx_r_x0, idx, axis=0)
  by_r_z0  = np.take( by_r_z0, idx, axis=0)
  by_r_y0  = np.take( by_r_y0, idx, axis=0)
  by_r_x0  = np.take( by_r_x0, idx, axis=0)
  upe_r_z0 = np.sqrt(ux_r_z0**2 + uy_r_z0**2)
  upe_r_y0 = np.sqrt(ux_r_y0**2 + uy_r_y0**2)
  upe_r_x0 = np.sqrt(ux_r_x0**2 + uy_r_x0**2)
  bpe_r_z0 = np.sqrt(bx_r_z0**2 + by_r_z0**2)
  bpe_r_y0 = np.sqrt(bx_r_y0**2 + by_r_y0**2)
  bpe_r_x0 = np.sqrt(bx_r_x0**2 + by_r_x0**2)
  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
  movie_3d(tt_fld, phi_r_z0, phi_r_y0, phi_r_x0, xx_fld, yy_fld, zz_fld, xlab=r'$x/\rho_\mathrm{H}$', ylab=r'$y/\rho_\mathrm{H}$', zlab='$'+zlab+'$', title=r'$\Phi$', cmp='diverging', save=outdir+'phi_anim.gif')
  movie_3d(tt_fld, omg_r_z0, omg_r_y0, omg_r_x0, xx_fld, yy_fld, zz_fld, xlab=r'$x/\rho_\mathrm{H}$', ylab=r'$y/\rho_\mathrm{H}$', zlab='$'+zlab+'$', title=r'$\nbl_\+^2\Phi$', cmp='diverging', save=outdir+'omg_anim.gif')
  movie_3d(tt_fld, psi_r_z0, psi_r_y0, psi_r_x0, xx_fld, yy_fld, zz_fld, xlab=r'$x/\rho_\mathrm{H}$', ylab=r'$y/\rho_\mathrm{H}$', zlab='$'+zlab+'$', title=r'$\Psi$', cmp='diverging', save=outdir+'psi_anim.gif')
  movie_3d(tt_fld, jpa_r_z0, jpa_r_y0, jpa_r_x0, xx_fld, yy_fld, zz_fld, xlab=r'$x/\rho_\mathrm{H}$', ylab=r'$y/\rho_\mathrm{H}$', zlab='$'+zlab+'$', title=r'$\nbl_\+^2\Psi$', cmp='diverging', save=outdir+'jpa_anim.gif')
  movie_3d(tt_fld, upe_r_z0, upe_r_y0, upe_r_x0, xx_fld, yy_fld, zz_fld, xlab=r'$x/\rho_\mathrm{H}$', ylab=r'$y/\rho_\mathrm{H}$', zlab='$'+zlab+'$', title=r'$|\bm{u}_\+|$', cmp='sequential', save=outdir+'upe_anim.gif')
  movie_3d(tt_fld, bpe_r_z0, bpe_r_y0, bpe_r_x0, xx_fld, yy_fld, zz_fld, xlab=r'$x/\rho_\mathrm{H}$', ylab=r'$y/\rho_\mathrm{H}$', zlab='$'+zlab+'$', title=r'$|\delta\bm{B}_\+|$', cmp='sequential', save=outdir+'bpe_anim.gif')

#------------------#
#   output data    #
#------------------#
from scipy.io import savemat

savemat(outdir + 'grid' , {'xx':xx, 'yy':yy, 'zz':zz})
savemat(outdir + 'u_r' , {
                            'tt'       :tt[final_fld_idx],
                            'ux_r_z0'  :ux_r_z0 [final_fld_idx,:,:]/np.sqrt(upe2_sum[final_fld_idx]),
                            'ux_r_x0'  :ux_r_x0 [final_fld_idx,:,:]/np.sqrt(upe2_sum[final_fld_idx]),
                            'ux_r_y0'  :ux_r_y0 [final_fld_idx,:,:]/np.sqrt(upe2_sum[final_fld_idx]),
                            #
                            'uy_r_z0'  :uy_r_z0 [final_fld_idx,:,:]/np.sqrt(upe2_sum[final_fld_idx]),
                            'uy_r_x0'  :uy_r_x0 [final_fld_idx,:,:]/np.sqrt(upe2_sum[final_fld_idx]),
                            'uy_r_y0'  :uy_r_y0 [final_fld_idx,:,:]/np.sqrt(upe2_sum[final_fld_idx]),
                            #
                            'upa_r_z0' :upa_r_z0[final_fld_idx,:,:]/np.sqrt(upa2_sum[final_fld_idx]),
                            'upa_r_x0' :upa_r_x0[final_fld_idx,:,:]/np.sqrt(upa2_sum[final_fld_idx]),
                            'upa_r_y0' :upa_r_y0[final_fld_idx,:,:]/np.sqrt(upa2_sum[final_fld_idx]),
                          })
savemat(outdir + 'b_r' , {
                            'tt'       :tt[final_fld_idx],
                            'bx_r_z0'  :bx_r_z0 [final_fld_idx,:,:]/np.sqrt(bpe2_sum[final_fld_idx]),
                            'bx_r_x0'  :bx_r_x0 [final_fld_idx,:,:]/np.sqrt(bpe2_sum[final_fld_idx]),
                            'bx_r_y0'  :bx_r_y0 [final_fld_idx,:,:]/np.sqrt(bpe2_sum[final_fld_idx]),
                            #
                            'by_r_z0'  :by_r_z0 [final_fld_idx,:,:]/np.sqrt(bpe2_sum[final_fld_idx]),
                            'by_r_x0'  :by_r_x0 [final_fld_idx,:,:]/np.sqrt(bpe2_sum[final_fld_idx]),
                            'by_r_y0'  :by_r_y0 [final_fld_idx,:,:]/np.sqrt(bpe2_sum[final_fld_idx]),
                            #
                            'bpa_r_z0' :bpa_r_z0 [final_fld_idx,:,:]/np.sqrt(bpa2_sum[final_fld_idx]),
                            'bpa_r_x0' :bpa_r_x0 [final_fld_idx,:,:]/np.sqrt(bpa2_sum[final_fld_idx]),
                            'bpa_r_y0' :bpa_r_y0 [final_fld_idx,:,:]/np.sqrt(bpa2_sum[final_fld_idx]),
                          })
savemat(outdir + 'omg_r' , {
                            'tt'       :tt[final_fld_idx],
                            'omg_r_z0' :omg_r_z0 [final_fld_idx,:,:]/np.sqrt(np.average(omg_r_z0**2)),
                            'omg_r_x0' :omg_r_x0 [final_fld_idx,:,:]/np.sqrt(np.average(omg_r_z0**2)),
                            'omg_r_y0' :omg_r_y0 [final_fld_idx,:,:]/np.sqrt(np.average(omg_r_z0**2)),
                          })
savemat(outdir + 'jpa_r' , {
                            'tt'       :tt[final_fld_idx],
                            'jpa_r_z0' :jpa_r_z0 [final_fld_idx,:,:]/np.sqrt(np.average(jpa_r_z0**2)),
                            'jpa_r_x0' :jpa_r_x0 [final_fld_idx,:,:]/np.sqrt(np.average(jpa_r_z0**2)),
                            'jpa_r_y0' :jpa_r_y0 [final_fld_idx,:,:]/np.sqrt(np.average(jpa_r_z0**2)),
                          })
