# -*- coding: utf-8 -*-
# import numpy as np
# from mayavi import mlab
# from tvtk.pyface.light_manager import CameraLight
from parula import parula_map

from load import *
from fft import *
import sys
sys.path.append('../')
from plots import *

#-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-#
#-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-#
#                                                                                                         #
#  Note that the x and y axis in the plotted figure is swapped from those in the code. This is because,   #
#   we swapped x and y in the code to employ remapping in shearing coordinate [Umurhan & Regev 2004]      #
#                                                                                                         #
#-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-#
#-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-#

print('\nplotting fields\n')
outdir = './fig_fields/'

####### ignore these points #########
ignored_points_f = [0]
elongation = 1

####### load 2decomp file #########
# Load field in spectral coordinate
tt_f = np.transpose(np.loadtxt(input_dir+'out2d/time.dat'+restart_num))[1]
nt_f = tt_f.size

ux = np.transpose(np.fromfile(input_dir+'out2d/ux.dat'+restart_num).reshape(nt_f, nky, nkz, nkx, 2))
ux = ux[0, :, :, :, final_idx] + 1j*ux[1, :, :, :, final_idx]
uy = np.transpose(np.fromfile(input_dir+'out2d/uy.dat'+restart_num).reshape(nt_f, nky, nkz, nkx, 2))
uy = uy[0, :, :, :, final_idx] + 1j*uy[1, :, :, :, final_idx]
uz = np.transpose(np.fromfile(input_dir+'out2d/uz.dat'+restart_num).reshape(nt_f, nky, nkz, nkx, 2))
uz = uz[0, :, :, :, final_idx] + 1j*uz[1, :, :, :, final_idx]

bx = np.transpose(np.fromfile(input_dir+'out2d/bx.dat'+restart_num).reshape(nt_f, nky, nkz, nkx, 2))
bx = bx[0, :, :, :, final_idx] + 1j*bx[1, :, :, :, final_idx]
by = np.transpose(np.fromfile(input_dir+'out2d/by.dat'+restart_num).reshape(nt_f, nky, nkz, nkx, 2))
by = by[0, :, :, :, final_idx] + 1j*by[1, :, :, :, final_idx]
bz = np.transpose(np.fromfile(input_dir+'out2d/bz.dat'+restart_num).reshape(nt_f, nky, nkz, nkx, 2))
bz = bz[0, :, :, :, final_idx] + 1j*bz[1, :, :, :, final_idx]
print (bx[1, 1, 0])

####### Inverse FFT #########
print (nlx, nlz, nly)

f = np.zeros([nlx, nlz, int(nly/2 + 1)])
f[1:nkx, 1:nkz, 1:nky] = ux[1:nkx, 1:nkz, 1:nky]

print (np.fft.irfftn(f).shape)
# print (np.irfftn())


# print (np.transpose(ux, axes=(0, 2, 1)).shape)
# print (np.fft.irfftn(np.transpose(ux, axes=(0, 2, 1))).shape)

quit()

psi = np.transpose(np.fromfile(input_dir+'psi.dat'+restart_num).reshape(nt_f, nkz, nky, nkx, 2))
psi = psi[0, :, :, :, final_idx] + 1j*psi[1, :, :, :, final_idx]

omg = np.zeros(phi.shape, dtype=complex)
for i in np.arange(0, nkx):
	for j in np.arange(0, nky):
		omg[i, j, :] = -(ky[j]**2 + kx[i]**2)*phi[i, j, :]

jpa = np.zeros(psi.shape, dtype=complex)
for i in np.arange(0, nkx):
	for j in np.arange(0, nky):
		jpa[i, j, :] = -(ky[j]**2 + kx[i]**2)*psi[i, j, :]


####### Inverse FFT #########
phi = fft_fwd3d(phi[:, :, :])
omg = fft_fwd3d(omg[:, :, :])
psi = fft_fwd3d(psi[:, :, :])
jpa = fft_fwd3d(jpa[:, :, :])
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

#--------------------------------------------------------#
#                   plot final snapshot                  #
#--------------------------------------------------------#
def mayavi_3d(u, x, y, z, cmp='RdBu', save=False):
  mlab.clf()
  # fig = mlab.figure()
  fig = mlab.figure(size=(600, 500), bgcolor=(1, 1, 1), fgcolor=(0 ,0 ,0))
  umin = min(np.min(u[-1, :, :]), np.min(u[:, -1, :]), np.min(u[:, :, -1]))
  umax = max(np.max(u[-1, :, :]), np.max(u[:, -1, :]), np.max(u[:, :, -1]))
  z = elongation*z
  s=mlab.pipeline.surface(mlab.pipeline.scalar_field(u), extent=[x.min(), x.max(), y.min(), y.max(), z.min(), z.max()], vmin=umin, vmax=umax, colormap=cmp)
  s.module_manager.scalar_lut_manager.reverse_lut = True

  ax=mlab.axes(xlabel='y', ylabel='x', zlabel='z', line_width=1)
  ax.title_text_property.font_family = 'times'
  ax.title_text_property.font_size = 10
  ax.title_text_property.bold = 0
  ax.label_text_property.font_family = 'times'
  ax.label_text_property.font_size = 0
  ax.label_text_property.bold = 0

  mlab.outline(extent=[x.min(), x.max(), y.min(), y.max(), z.min(), z.max()], line_width=1)

  cb = mlab.colorbar(orientation='vertical')
  cb.title_text_property.font_family = 'times'
  cb.title_text_property.font_size = 5
  cb.title_text_property.bold = 1
  cb.label_text_property.font_family = 'times'
  cb.label_text_property.font_size = 10
  cb.label_text_property.bold = 1

  mlab.orientation_axes(xlabel='y', ylabel='x', zlabel='z')
  mlab.view(azimuth=88.0, elevation=0.0, distance=22.0, focalpoint=(-7, -7, 0))
  mlab.pitch(20)
  mlab.yaw(22)
  mlab.savefig(save)
  mlab.close()


mayavi_3d(phi, xx, yy, zz, save=outdir + 'phi_final_3d.png' )
mayavi_3d(phi, xx, yy, zz, save=outdir + 'phi_final_3d.png' )
mayavi_3d(omg, xx, yy, zz, save=outdir + 'omg_final_3d.png' )
mayavi_3d(psi, xx, yy, zz, save=outdir + 'psi_final_3d.png' )
mayavi_3d(jpa, xx, yy, zz, save=outdir + 'jpa_final_3d.png' )


#--------------------------------------------------------#
#                      plot 2D cuts                      #
#--------------------------------------------------------#
plot_3d_3cuts(np.transpose(phi, axes=(1,0,2)), yy, xx, zz, xslice=int(nly/2), yslice=int(nlx/2), zslice=int(nlz/2), xlab='$x$', ylab='$y$', zlab='$z$', title=r'$\Phi (t = %.2E)$' % tt_f[final_idx], cmp='RdBu_r', save=outdir+'phi.pdf')
plot_3d_3cuts(np.transpose(psi, axes=(1,0,2)), yy, xx, zz, xslice=int(nly/2), yslice=int(nlx/2), zslice=int(nlz/2), xlab='$x$', ylab='$y$', zlab='$z$', title=r'$\Psi (t = %.2E)$' % tt_f[final_idx], cmp='RdBu_r', save=outdir+'psi.pdf')
plot_3d_3cuts(np.transpose(omg, axes=(1,0,2)), yy, xx, zz, xslice=int(nly/2), yslice=int(nlx/2), zslice=int(nlz/2), xlab='$x$', ylab='$y$', zlab='$z$', title=r'$\nabla_\perp^2\Phi (t = %.2E)$' % tt_f[final_idx], cmp='RdBu_r', save=outdir+'omg.pdf')
plot_3d_3cuts(np.transpose(jpa, axes=(1,0,2)), yy, xx, zz, xslice=int(nly/2), yslice=int(nlx/2), zslice=int(nlz/2), xlab='$x$', ylab='$y$', zlab='$z$', title=r'$\nabla_\perp^2\Psi (t = %.2E)$' % tt_f[final_idx], cmp='RdBu_r', save=outdir+'jpa.pdf')
