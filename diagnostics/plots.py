# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use('Agg')

import numpy as np
import pylab
from matplotlib import animation
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import *
import matplotlib.ticker as ticker
from parula import parula_map
from cmcrameri import cm

interpolation        = 'none' # or 'nearest'
default_colormap_seq = cm.batlow 
default_colormap_div = cm.vik    
default_colormap_log = cm.lajolla

# setup some plot defaults
linewidth = 2
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', serif='Times')
plt.rc('font', size=30)
plt.rc('axes', linewidth=linewidth)
plt.rc('xtick', top=True)
plt.rc('xtick.major', width=linewidth)
plt.rc('xtick.major', size=15)
plt.rc('xtick.minor', width=linewidth)
plt.rc('xtick.minor', size=8)
plt.rc('ytick', right=True)
plt.rc('ytick.major', width=linewidth)
plt.rc('ytick.major', size=10)
plt.rc('ytick.minor', width=linewidth)
plt.rc('ytick.minor', size=5)
plt.rc('xtick', direction='in')
plt.rc('ytick', direction='in')
from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler(
        color=["#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#999999"]) # Okabe-Ito colorscheme
rcParams.update({'figure.autolayout': True})
from latex_preamble import preamble
rcParams['text.latex.preamble'] = preamble  


def plot_1d(x, y, xlab, title='', ylab='', ls='', ymin='', ymax='', term=True, save=False):
  fig = plt.figure(figsize=(9, 8))
  plt.plot(x, y, ls, lw=2)
  plt.xlabel(xlab)
  if len(ylab) > 0:
    plt.ylabel(ylab)
  if len(title) > 0:
    plt.title(title)

  # force exponential yticks
  ax = plt.gca()
  # ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

  ax.tick_params(which='both', direction='in')
  if ymin != '' and ymax != '':
    ax.set_ylim([ymin, ymax])
  elif ymin != '':
    ax.set_ylim(ymin=ymin)
  elif ymax != '':
    ax.set_ylim(ymax=ymax)

  if save:
    print ('--   ' + title + '   --')
    plt.savefig(save)

  if term:
    plt.show()
  return fig


def movie_1d(t, x, y, nframes, xlab='', ylab='', title='', save=False):

  fig = plt.figure(figsize=(9, 8))
  def animate(i):
    plt.cla()
    plt.plot(x, y[i], lw=2)
    plt.xlim(min(x), max(x))
    plt.ylim(min(y[i,:]), max(y[i,:]))
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title + ' $(t/\\tau_\mathrm{A} = $ %.2E' % t[i] + '$)$')
    ax = plt.gca()
    ax.tick_params(which='both', direction='in')

  ani = animation.FuncAnimation(fig, animate, np.arange(0, nframes), interval=1)
  if save:
    print ('--   ' + title + '   --')
    if save.split('.')[-1] == 'gif':
      ani.save(save, writer='imagemagick', fps=30);
    if save.split('.')[-1] == 'mp4':
      # ani.save(save, writer=animation.FFMpegWriter());
      ani.save(save, writer='ffmpeg', fps=30)
  # plt.show()


def plot_1d_many(xs, ys, xlab, legends, ls, legendloc='upper left', title='', ylab='', ymin='', ymax='', term=True, save=False):
  fig = plt.figure(figsize=(9, 8))

  for i, y in enumerate(ys):
    plt.plot(xs[i], ys[i], ls[i], lw=2, label=legends[i])

  legenaxd = plt.legend(fontsize=20, loc=legendloc)
  plt.xlabel(xlab)
  if len(ylab) > 0:
    plt.ylabel(ylab)
  if len(title) > 0:
    plt.title(title)

  if ymin != '' and ymax != '':
    plt.ylim([ymin, ymax])
  elif ymin != '':
    plt.ylim(ymin=ymin)
  elif ymax != '':
    plt.ylim(ymax=ymax)

  # force exponential yticks
  ax = plt.gca()
  # ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

  ax.tick_params(which='both', direction='in')

  if save:
    print ('--   ' + title + '   --')
    plt.savefig(save)

  if term:
    plt.show()
  return fig


def plot_1d_many_average(xs, ys, avg_start, avg_end, xlab, legends, ls, legendloc='upper left', title='', ylab='', term=True, save=False):
  fig = plt.figure(figsize=(9, 8))

  for i, y in enumerate(ys):
    plt.plot(xs[i], ys[i], ls[i], lw=2, label=legends[i])

  legenaxd = plt.legend(fontsize=20, loc=legendloc)
  plt.xlabel(xlab)
  if len(ylab) > 0:
    plt.ylabel(ylab)
  if len(title) > 0:
    plt.title(title)

  # force exponential yticks
  ax = plt.gca()
  # ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

  ax.tick_params(which='both', direction='in')

  pylab.axvspan(avg_start, avg_end, facecolor='r', alpha=0.2)

  if save:
    print ('--   ' + title + '   --')
    plt.savefig(save)

  if term:
    plt.show()
  return fig


def plot_log1d(x, y, xlab, title='', ylab='', ls='', term=True, save=False):
  fig = plt.figure(figsize=(9, 8))
  plt.loglog(x, y, ls, lw=2)
  plt.xlabel(xlab)
  if len(ylab) > 0:
    plt.ylabel(ylab)
  if len(title) > 0:
    plt.title(title)

  ax = plt.gca()
  ax.tick_params(which='both', direction='in')

  if save:
    print ('--   ' + title + '   --')
    plt.savefig(save)

  if term:
    plt.show()
  return fig


def plot_log1d_many(xs, ys, xlab, legends, ls, xmin='', xmax='', ymin='', ymax='', legendloc='upper left', title='', ylab='', term=True, save=False):
  fig = plt.figure(figsize=(9, 8))

  for i, y in enumerate(ys):
    plt.loglog(xs[i], ys[i], ls[i], lw=2, label=legends[i])

  legenaxd = plt.legend(fontsize=20, loc=legendloc)
  plt.xlabel(xlab)
  if xmin != '':
    plt.xlim(xmin=xmin)
  if xmax != '':
    plt.xlim(xmax=xmax)
  if ymin != '':
    plt.ylim(ymin=ymin)
  if ymax != '':
    plt.ylim(ymax=ymax)
  if len(ylab) > 0:
    plt.ylabel(ylab)
  if len(title) > 0:
    plt.title(title)

  # force exponential yticks
  ax = plt.gca()
  ax.tick_params(which='both', direction='in')


  if save:
    print ('--   ' + title + '   --')
    plt.savefig(save)

  if term:
    plt.show()
  return fig


def plot_semilogx1d_many(xs, ys, xlab, legends, ls, xmin='', xmax='', ymin='', ymax='', legendloc='upper left', title='', ylab='', term=True, save=False):
  fig = plt.figure(figsize=(9, 8))

  for i, y in enumerate(ys):
    plt.semilogx(xs[i], ys[i], ls[i], lw=2, label=legends[i])

  legenaxd = plt.legend(fontsize=20, loc=legendloc)
  plt.xlabel(xlab)
  if xmin != '':
    plt.xlim(xmin=xmin)
  if xmax != '':
    plt.xlim(xmax=xmax)
  if ymin != '':
    plt.ylim(ymin=ymin)
  if ymax != '':
    plt.ylim(ymax=ymax)
  if len(ylab) > 0:
    plt.ylabel(ylab)
  if len(title) > 0:
    plt.title(title)

  # force exponential yticks
  ax = plt.gca()
  ax.tick_params(which='both', direction='in')


  if save:
    print ('--   ' + title + '   --')
    plt.savefig(save)

  if term:
    plt.show()
  return fig


def movie_log1d_many(t, xs, ys, legends, ls, legendloc, xlab='', title='', save=False):

  fig = plt.figure(figsize=(9, 8))
  def animate(i):
    plt.cla()
    for j, x in enumerate(xs):
      plt.loglog(xs[j], ys[i, j], ls[j], lw=2, label=legends[j])

    ncol = 1 if len(legends) < 4 else 2
    legenaxd = plt.legend(fontsize=20, loc=legendloc, ncol=ncol)
    plt.xlabel(xlab)

    # force exponential yticks
    ax = plt.gca()
    ax.tick_params(which='both', direction='in')

  ani = animation.FuncAnimation(fig, animate, np.arange(0, t.size), interval=1)
  if save:
    print ('--   ' + title + '   --')
    if save.split('.')[-1] == 'gif':
      ani.save(save, writer='imagemagick', fps=30);
    if save.split('.')[-1] == 'mp4':
      # ani.save(save, writer=animation.FFMpegWriter());
      ani.save(save, writer='ffmpeg', fps=30)
  # plt.show()


def plot_2d(u, xin, yin, umin=None, umax=None, xlab='', ylab='', title='', cmp=default_colormap_seq, save=False):
  if umin is None:
    umin = u.min()
  if umax is None:
    umax = u.max()

  if cmp == 'sequential':
    cmp = default_colormap_seq
  if cmp == 'diverging':
    cmp = default_colormap_div

  fig = plt.figure(figsize=(9, 8))
  ax = plt.gca()
  ax.tick_params(which='both', direction='in')
  x, y = np.meshgrid(xin, yin)
  im = ax.imshow(u, cmap=cmp, vmin=umin, vmax=umax,
             extent=[x.min(), x.max(), y.min(), y.max()],
             interpolation=interpolation, origin='lower', aspect='equal')
  if contour:
    plt.contour(x, y, u, colors='k', vmin=umin, vmax=umax,
               extent=[x.min(), x.max(), y.min(), y.max()],
               interpolation=interpolation, origin='lower', aspect=aspect)
  # try:
    # plt.contour(x, y, u, linewidths=2)
  # except:
    # pass
  plt.axis([x.min(), x.max(), y.min(), y.max()])
  plt.xlabel(xlab)
  plt.ylabel(ylab)
  plt.locator_params(nbins=10)
  plt.title(title)
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cb = plt.colorbar(im, cax=cax)
  cb.formatter.set_powerlimits((0, 0))
  cb.ax.yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
  cb.ax.yaxis.set_offset_position('left')  
  cb.locator = ticker.MaxNLocator(nbins=9)
  cb.update_ticks()

  if save:
    print ('--   ' + title + '   --')
    plt.savefig(save)

  plt.show()
  return fig


def plot_SF2(SF2b, SF2u, lpar, lper, xlab='', ylab='', zlab='', title='', cmp=default_colormap_seq, save=False):

  fig, axes = plt.subplots(1, 2, figsize=(16, 9))

  ax = axes[0]
  ax.contour(lpar, lper, SF2b, levels=20, linewidths=3, origin='lower', aspect='equal')
  ax.set_title(r'$\mr{SF}^2_b$')
  ax.set_xlabel(xlab)
  ax.set_ylabel(ylab)
  ax.locator_params(axis="x", nbins=4)
  ax.locator_params(axis="y", nbins=4)

  ax = axes[1]
  ax.contour(lpar, lper, SF2u, levels=20, linewidths=3, origin='lower', aspect='equal')
  ax.tick_params(labelleft='off')
  plt.setp(ax.get_yticklabels(), visible=False)
  ax.set_title(r'$\mr{SF}^2_u$')
  ax.set_xlabel(xlab)
  ax.locator_params(axis="x", nbins=4)
  ax.locator_params(axis="y", nbins=4)

  plt.suptitle(title)
  plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)

  # plt.tight_layout()

  if save:
    print ('--   ' + title + '   --')
    plt.savefig(save)

  plt.show()
  return fig


def plot_log2d(u, xin, yin, umin=None, umax=None, xlab='', ylab='', title='', cmp=default_colormap_log, contour=False, save=False, aux=np.nan, aux_label=''):
  if umin is None:
    umin = np.log10(u).min()
  if umax is None:
    umax = np.log10(u).max()

  fig = plt.figure(figsize=(9, 8))
  ax = plt.gca()
  ax.tick_params(which='both', direction='in')
  x, y = np.meshgrid(xin, yin)
  aspect = (x.min() - x.max()) / (y.min() - y.max())

  im = ax.imshow(np.log10(u), cmap=cmp, vmin=umin, vmax=umax,
             extent=[x.min(), x.max(), y.min(), y.max()],
             interpolation=interpolation, origin='lower', aspect=aspect)
  if contour:
    plt.contour(x, y, np.log10(u), colors='k', vmin=umin, vmax=umax,
               extent=[x.min(), x.max(), y.min(), y.max()],
               interpolation=interpolation, origin='lower', aspect=aspect)
  # try:
    # plt.contour(x, y, u, linewidths=2)
  # except:
    # pass
  if np.isfinite(aux).any():
    ax.plot(xin, aux, 'k--', label=aux_label)
    legenaxd = plt.legend(fontsize=20, loc='lower right')

  plt.axis([x.min(), x.max(), y.min(), y.max()])
  plt.xlabel(xlab)
  plt.ylabel(ylab)
  plt.locator_params(nbins=10)
  plt.title(title)
  ax.set_xscale('log')
  ax.set_yscale('log')
  divider = make_axes_locatable(ax)
  plt.axis('auto')
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cb = plt.colorbar(im, cax=cax)
  cb.formatter.set_powerlimits((0, 0))
  cb.ax.yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
  cb.ax.yaxis.set_offset_position('left')  
  cb.locator = ticker.MaxNLocator(nbins=9)
  cb.update_ticks()

  if save:
    print ('--   ' + title + '   --')
    plt.savefig(save)

  plt.show()
  return fig


def plot_symlog2d(u, xin, yin, umin=None, umax=None, xlab='', ylab='', title='', cmp=default_colormap_log, contour=False, save=False, aux=np.nan, aux_label=''):
  if umin is None:
    umin = np.log10(u).min()
  if umax is None:
    umax = np.log10(u).max()

  aspect = (xin.min() - xin.max()) / (yin.min() - yin.max())
	
  fig = plt.figure(figsize=(9, 8))
  ax = plt.gca()
  ax.tick_params(which='both', direction='in')
  x, y = np.meshgrid(xin, yin)

  im = ax.imshow(np.log10(u), cmap=cmp, vmin=umin, vmax=umax,
             extent=[x.min(), x.max(), y.min(), y.max()],
             interpolation=interpolation, origin='lower')
  if contour:
    plt.contour(x, y, np.log10(u), colors='k', vmin=umin, vmax=umax,
               extent=[x.min(), x.max(), y.min(), y.max()],
               interpolation=interpolation, origin='lower')
  # try:
    # plt.contour(x, y, u, linewidths=2)
  # except:
    # pass
  if np.isfinite(aux).any():
    ax.plot(xin, aux, 'k--', label=aux_label)
    legenaxd = plt.legend(fontsize=20, loc='lower right')

  plt.axis([x.min(), x.max(), y.min(), y.max()])
  plt.xlabel(xlab)
  plt.ylabel(ylab)
  plt.locator_params(nbins=10)
  plt.title(title)
  ax.set_xscale('symlog')
  ax.set_yscale('symlog')
  plt.axis('auto')
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cb = plt.colorbar(im, cax=cax)
  cb.formatter.set_powerlimits((0, 0))
  cb.ax.yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
  cb.ax.yaxis.set_offset_position('left')  
  cb.locator = ticker.MaxNLocator(nbins=9)
  cb.update_ticks()

  if save:
    print ('--   ' + title + '   --')
    plt.savefig(save)

  plt.show()
  return fig


def plot_3d(u_z0, u_y0, u_x0, xin, yin, zin, xlab='', ylab='', zlab='', title='', cmp=default_colormap_seq, save=False):

  if cmp == 'sequential':
    cmp = default_colormap_seq
  if cmp == 'diverging':
    cmp = default_colormap_div

  fig = plt.figure(figsize=(20, 7))

  # z = 0 slice
  umin = u_z0.min()
  umax = u_z0.max()

  plt.subplot(131)
  ax = plt.gca()
  ax.tick_params(which='both', direction='in')
  x, y = np.meshgrid(xin, yin)
  im = ax.imshow(u_z0, cmap=cmp, vmin=umin, vmax=umax,
             extent=[x.min(), x.max(), y.min(), y.max()],
             interpolation=interpolation, origin='lower', aspect='equal')

  # streamline
  try:
    if u_z0_vec is not None:
      plt.streamplot(x, y, u_z0_vec[0], u_z0_vec[1], density=streamline_density, color=streamline_color, arrowstyle='-', linewidth=streamline_width)
  except:
    pass

  # try:
    # plt.contour(x, y, u[:, :, zslice], linewidths=2)
  # except:
    # pass
  plt.axis([x.min(), x.max(), y.min(), y.max()])
  plt.xlabel(xlab)
  plt.ylabel(ylab)
  plt.locator_params(nbins=10)
  plt.title(title + r'$(z = 0)$')
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cb = plt.colorbar(im, cax=cax)
  cb.formatter.set_powerlimits((0, 0))
  cb.ax.yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
  cb.ax.yaxis.set_offset_position('left')  
  cb.locator = ticker.MaxNLocator(nbins=9)
  cb.update_ticks()


  # y = 0 slice
  umin = u_y0.min()
  umax = u_y0.max()

  plt.subplot(132)
  ax = plt.gca()
  ax.tick_params(which='both', direction='in')
  x, z = np.meshgrid(xin, zin)
  im = ax.imshow(u_y0, cmap=cmp, vmin=umin, vmax=umax,
             extent=[x.min(), x.max(), z.min(), z.max()],
             interpolation=interpolation, origin='lower', aspect='equal')
  # try:
    # plt.contour(x, z, np.transpose(u[:, yslice, :]), linewidths=2)
  # except:
    # pass
  plt.axis([x.min(), x.max(), z.min(), z.max()])
  plt.xlabel(xlab)
  plt.ylabel(zlab)
  plt.locator_params(nbins=10)
  plt.title(title + r'$(y = 0)$')
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cb = plt.colorbar(im, cax=cax)
  cb.formatter.set_powerlimits((0, 0))
  cb.ax.yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
  cb.ax.yaxis.set_offset_position('left')  
  cb.locator = ticker.MaxNLocator(nbins=9)
  cb.update_ticks()


  # x = 0 slice
  umin = u_x0.min()
  umax = u_x0.max()

  plt.subplot(133)
  ax = plt.gca()
  ax.tick_params(which='both', direction='in')
  y, z = np.meshgrid(yin, zin)
  im = ax.imshow(u_x0, cmap=cmp, vmin=umin, vmax=umax,
             extent=[y.min(), y.max(), z.min(), z.max()],
             interpolation=interpolation, origin='lower', aspect='equal')

  # streamline
  try:
    if u_x0_vec is not None:
      plt.streamplot(y, z, u_x0_vec[0], u_x0_vec[1], density=streamline_density, color=streamline_color, arrowstyle='-', linewidth=streamline_width)
  except:
    pass
  # try:
    # plt.contour(y, z, np.transpose(u[xslice, :, :]), linewidths=2)
  # except:
    # pass
  plt.axis([y.min(), y.max(), z.min(), z.max()])
  plt.xlabel(ylab)
  plt.ylabel(zlab)
  plt.locator_params(nbins=10)
  plt.title(title + r'$(x = 0)$')
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cb = plt.colorbar(im, cax=cax)
  cb.formatter.set_powerlimits((0, 0))
  cb.ax.yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
  cb.ax.yaxis.set_offset_position('left')  
  cb.locator = ticker.MaxNLocator(nbins=9)
  cb.update_ticks()

  # plt.tight_layout()

  if save:
    print ('--   ' + title + '   --')
    plt.savefig(save)

  plt.show()
  return fig


def movie_3d(t, u_z0, u_y0, u_x0, xin, yin, zin, xlab='', ylab='', zlab='', title='', cmp=default_colormap_seq, save=False):

  if cmp == 'sequential':
    cmp = default_colormap_seq
  if cmp == 'diverging':
    cmp = default_colormap_div

  fig, [ax1, ax2, ax3] = plt.subplots(1, 3, figsize=(20, 7))

  ims = []
  def animate(i):

    # z = 0 slice
    ax1.cla()
    ax1.tick_params(which='both', direction='in')
    ax1.axis([xin.min(), xin.max(), yin.min(), yin.max()])
    ax1.set_xlim([xin.min(), xin.max()])
    ax1.set_ylim([yin.min(), yin.max()])
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab)
    ax1.locator_params(nbins=10)
    ax1.set_title(title + r'$(z = 0)$')

    umin = u_z0.min()
    umax = u_z0.max()
    x, y = np.meshgrid(xin, yin)
    im1 = ax1.imshow(u_z0[i, :, :], cmap=cmp,
               extent=[x.min(), x.max(), y.min(), y.max()],
               interpolation=interpolation, origin='lower', aspect='equal', animated=True)

    # streamline
    try:
      if u_z0_vec is not None:
        ax1.streamplot(x, y, u_z0_vec[0][i, :, :], u_z0_vec[1][i, :, :], density=streamline_density, color=streamline_color, arrowstyle='-', linewidth=streamline_width)
    except:
      pass
    # divider = make_axes_locatable(ax1)
    # cax = divider.append_axes("right", size="5%", pad=0.05)
    # cb = plt.colorbar(im1, cax=cax)
    # cb.formatter.set_powerlimits((0, 0))
    # cb.ax.yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
    # cb.ax.yaxis.set_offset_position('left')  
    # cb.locator = ticker.MaxNLocator(nbins=9)
    # cb.update_ticks()

    # y = 0 slice
    ax2.cla()
    ax2.tick_params(which='both', direction='in')
    ax2.axis([xin.min(), xin.max(), zin.min(), zin.max()])
    ax2.set_xlim([xin.min(), xin.max()])
    ax2.set_ylim([zin.min(), zin.max()])
    ax2.set_xlabel(xlab)
    ax2.set_ylabel(zlab)
    ax2.locator_params(nbins=10)
    ax2.set_title(title + r'$(y = 0)$')

    umin = u_y0.min()
    umax = u_y0.max()
    x, z = np.meshgrid(xin, zin)
    im2 = ax2.imshow(u_y0[i, :, :], cmap=cmp,
               extent=[x.min(), x.max(), z.min(), z.max()],
               interpolation=interpolation, origin='lower', aspect='equal')

    # streamline
    try:
      if u_y0_vec is not None:
        ax2.streamplot(x, z, u_y0_vec[0][i, :, :], u_y0_vec[1][i, :, :], density=streamline_density, color=streamline_color, arrowstyle='-', linewidth=streamline_width)
    except:
      pass
    # divider = make_axes_locatable(ax2)
    # cax = divider.append_axes("right", size="5%", pad=0.05)
    # cb = plt.colorbar(im2, cax=cax)
    # cb.formatter.set_powerlimits((0, 0))
    # cb.ax.yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
    # cb.ax.yaxis.set_offset_position('left')  
    # cb.locator = ticker.MaxNLocator(nbins=9)
    # cb.update_ticks()

    # x = 0 slice
    ax3.cla()
    ax3.tick_params(which='both', direction='in')
    ax3.axis([yin.min(), yin.max(), zin.min(), zin.max()])
    ax3.set_xlim([yin.min(), yin.max()])
    ax3.set_ylim([zin.min(), zin.max()])
    ax3.set_xlabel(ylab)
    ax3.set_ylabel(zlab)
    ax3.locator_params(nbins=10)
    ax3.set_title(title + r'$(x = 0)$')

    umin = u_x0.min()
    umax = u_x0.max()
    y, z = np.meshgrid(yin, zin)
    im3 = ax3.imshow(u_x0[i, :, :], cmap=cmp,
               extent=[y.min(), y.max(), z.min(), z.max()],
               interpolation=interpolation, origin='lower', aspect='equal')

    # streamline
    try:
      if u_x0_vec is not None:
        ax3.streamplot(y, z, u_x0_vec[0][i, :, :], u_x0_vec[1][i, :, :], density=streamline_density, color=streamline_color, arrowstyle='-', linewidth=streamline_width)
    except:
      pass
    # divider = make_axes_locatable(ax3)
    # cax = divider.append_axes("right", size="5%", pad=0.05)
    # cb = plt.colorbar(im3, cax=cax)
    # cb.formatter.set_powerlimits((0, 0))
    # cb.ax.yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
    # cb.ax.yaxis.set_offset_position('left')  
    # cb.locator = ticker.MaxNLocator(nbins=9)
    # cb.update_ticks()

    ims.append([im1, im2, im3])


  ani = animation.FuncAnimation(fig, animate, np.arange(0, t.size), interval=1)
  if save:
    print ('--   ' + title + '   --')
    if save.split('.')[-1] == 'gif':
      ani.save(save, writer='imagemagick', fps=30);
    if save.split('.')[-1] == 'mp4':
      # ani.save(save, writer=animation.FFMpegWriter());
      ani.save(save, writer='ffmpeg', fps=30)
  # plt.show()


def plot_3d_3cuts(u, xin, yin, zin, xslice, yslice, zslice, xlab='', ylab='', zlab='', title='', suptitle='', text='', cmp=default_colormap_log, save=False):

  if cmp == 'sequential':
    cmp = default_colormap_seq
  if cmp == 'diverging':
    cmp = default_colormap_div

  fig = plt.figure(figsize=(21, 7))

  # z = 0 slice
  umin = u[:, :, zslice].min()
  umax = u[:, :, zslice].max()

  plt.subplot(131)
  ax = plt.gca()
  ax.tick_params(which='both', direction='in')
  x, y = np.meshgrid(xin, yin)
  im = ax.imshow(np.transpose(u[:, :, zslice]), cmap=cmp, vmin=umin, vmax=umax,
             extent=[x.min(), x.max(), y.min(), y.max()],
             interpolation=interpolation, origin='lower', aspect='equal')
  # try:
    # plt.contour(x, y, u[:, :, zslice], linewidths=2)
  # except:
    # pass
  plt.suptitle(suptitle, fontsize=25)
  plt.axis([x.min(), x.max(), y.min(), y.max()])
  plt.xlabel(xlab)
  plt.ylabel(ylab)
  plt.locator_params(nbins=10)
  plt.title(title + r'$(z = %2.1f)$' % zin[zslice])
  plt.text(x.min() + (x.max()-x.min())*0.02, y.min() + (y.max()-y.min())*0.02, text, fontsize=15)
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cb = plt.colorbar(im, cax=cax)
  cb.formatter.set_powerlimits((0, 0))
  cb.ax.yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
  cb.ax.yaxis.set_offset_position('left')  
  cb.locator = ticker.MaxNLocator(nbins=9)
  cb.update_ticks()


  # y = 0 slice
  umin = u[:, yslice, :].min()
  umax = u[:, yslice, :].max()

  plt.subplot(132)
  ax = plt.gca()
  ax.tick_params(which='both', direction='in')
  x, z = np.meshgrid(xin, zin)
  im = ax.imshow(np.transpose(u[:, yslice, :]), cmap=cmp, vmin=umin, vmax=umax,
             extent=[x.min(), x.max(), z.min(), z.max()],
             interpolation=interpolation, origin='lower', aspect='equal')
  # try:
    # plt.contour(x, z, np.transpose(u[:, yslice, :]), linewidths=2)
  # except:
    # pass
  plt.axis([x.min(), x.max(), z.min(), z.max()])
  plt.xlabel(xlab)
  plt.ylabel(zlab)
  plt.locator_params(nbins=10)
  plt.title(title + r'$(y = %2.1f)$' % yin[yslice])
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cb = plt.colorbar(im, cax=cax)
  cb.formatter.set_powerlimits((0, 0))
  cb.ax.yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
  cb.ax.yaxis.set_offset_position('left')  
  cb.locator = ticker.MaxNLocator(nbins=9)
  cb.update_ticks()


  # x = 0 slice
  umin = u[xslice, :, :].min()
  umax = u[xslice, :, :].max()

  plt.subplot(133)
  ax = plt.gca()
  ax.tick_params(which='both', direction='in')
  y, z = np.meshgrid(yin, zin)
  im = ax.imshow(np.transpose(u[xslice, :, :]), cmap=cmp, vmin=umin, vmax=umax,
             extent=[y.min(), y.max(), z.min(), z.max()],
             interpolation=interpolation, origin='lower', aspect='equal')
  # try:
    # plt.contour(y, z, np.transpose(u[xslice, :, :]), linewidths=2)
  # except:
    # pass
  plt.axis([y.min(), y.max(), z.min(), z.max()])
  plt.xlabel(ylab)
  plt.ylabel(zlab)
  plt.locator_params(nbins=10)
  plt.title(title + r'$(x = %2.1f)$' % xin[xslice])
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cb = plt.colorbar(im, cax=cax)
  cb.formatter.set_powerlimits((0, 0))
  cb.ax.yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
  cb.ax.yaxis.set_offset_position('left')  
  cb.locator = ticker.MaxNLocator(nbins=9)
  cb.update_ticks()

  # plt.tight_layout()

  if save:
    print ('--   ' + title + '   --')
    plt.savefig(save)

  plt.show()
  return fig

