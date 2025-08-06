# -*- coding: utf-8 -*-
import warnings
warnings.filterwarnings('ignore')
#-------------------------------------------------------------#
#                      Matplotlib setting                     #
#-------------------------------------------------------------#
import pylab
import matplotlib as mpl
# mpl.use('Agg')
from matplotlib import animation
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import *
import matplotlib.ticker as ticker
from parula import parula_map
from scipy.io import loadmat

interpolation        = 'nearest' # or 'nearest'
default_colormap     = parula_map # or parula_map or 'viridis' for sequential colormaps and RdBu for diverging colormap
default_colormap_log = 'inferno' # or parula_map or 'viridis'

# setup some plot defaults
linewidth = 2
fontsize  = 30
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', serif='Times')
plt.rc('font', size=fontsize)
plt.rc('axes', linewidth=linewidth)
plt.rc('axes', labelsize=fontsize)
plt.rc('legend', fontsize=fontsize)
plt.rc('xtick', labelsize=fontsize)
plt.rc('xtick', top=True)
plt.rc('xtick.major', width=linewidth)
plt.rc('xtick.major', size=17)
plt.rc('xtick.minor', width=linewidth)
plt.rc('xtick.minor', visible=True)
plt.rc('xtick.minor', size=8)
plt.rc('ytick', labelsize=fontsize)
plt.rc('ytick', right=True)
plt.rc('ytick.major', width=linewidth)
plt.rc('ytick.major', size=17)
plt.rc('ytick.minor', width=linewidth)
plt.rc('ytick.minor', visible=True)
plt.rc('ytick.minor', size=8)
# plt.rc('xtick', labelbottom='off')
plt.rc('xtick', direction='in')
# plt.rc('ytick', labelleft='off')
plt.rc('ytick', direction='in')
rcParams.update({'figure.autolayout': True})
from latex_preamble import preamble
rcParams['text.latex.preamble'] = preamble  
tab10 = plt.get_cmap('tab10').colors
#-------------------------------------------------------------#

import glob
import numpy as np
from numpy.fft import fft, fftfreq

import re

def atoi(text):
  return int(text) if text.isdigit() else text

def natural_keys(text):
  return [ atoi(c) for c in re.split(r'(\d+)', text) ]

prefix = '../'

#---------------------------------------#
#              Load files
#---------------------------------------#
modes = np.loadtxt('../modes.txt')

filenames = sorted(glob.glob(prefix + '*.txt'), key=natural_keys)
nfld = len(filenames)

fig1, axes1 = plt.subplots(nfld//2, 2, figsize=(20, 17))
fig2, axes2 = plt.subplots(nfld//2, 2, figsize=(20, 17))

for ifld, filename in enumerate([filename for filename in filenames if not 'modes.txt' in filename]):
  fldname = filename.split(prefix)[1].split('.txt')[0]

  data = np.loadtxt(filename).T
  time = data[0]; time = time - time[0]

  print ('field = '+fldname+' @ t = [%.2e, %.2e]' % (time[0], time[-1]))

  nt = len(time)
  dt = np.unique(np.diff(time))[0]

  freq = 2.*np.pi*fftfreq(nt, dt)[:nt//2]

  # plot time history of a few modes
  ax1 = axes1[ifld%(axes1.shape[1]+1)][ifld//axes1.shape[0]]
  ax1.set_title(fldname)
  for imode in [0, 1, 2]:
    fldval  = data[2*imode+1] + 1j*data[2*imode+2]
    ax1.plot(time, np.real(fldval), lw=2, label=str(modes[imode]))
  ax1.set_xlabel(r'$t\Omega$')
  ax1.legend(frameon=False, loc='upper right', fontsize=25, handlelength=1.1, columnspacing=0.6)

  ax2 = axes2[ifld%(axes2.shape[1]+1)][ifld//axes2.shape[0]]
  ax2.set_title(fldname)
  for imode, mode in enumerate(modes):
    print ('mode = '+str(mode))

    fldval  = data[2*imode+1] + 1j*data[2*imode+2]
    fldvalk = fft(fldval*np.hanning(nt))

    if imode in [0, 1, 10]:
      ax2.loglog(freq, np.abs(np.real(fldvalk)[:nt//2]), label=str(modes[imode]))
  #  ax2.set_xlim(xmax=100.)
  ax2.set_xlabel(r'$\omega/\Omega$')
  ax2.legend(frameon=False, loc='upper right', fontsize=25, handlelength=1.1, columnspacing=0.6)

fig1.savefig('evolution.pdf')
fig2.savefig('spectra.pdf')
quit()


fig, axes = plt.subplots(3, 2, figsize=(20, 17))

# get Lz/lmd_MRI and Lz/H
lines = [line.rstrip('\n') for line in open(glob.glob('../../../calliope.out.std*')[0])]

line = [x for x in lines if 'Lz/lmd_MRI' in x]
lmd_MRI = float(line[0].split('Lz/lmd_MRI = ')[1])


#---------------------------------------#
#              First figure
#---------------------------------------#
ax = axes[0][0]
# load data
kp, \
u2,\
ux2,\
uy2,\
uz2,\
b2,\
bx2,\
by2,\
bz2,\
zp2,\
zm2 = np.loadtxt('../Ek_avg.txt').T

# Wtot = np.sum(upe2_perp) + np.sum(bpe2_perp)
# upe2_perp = upe2_perp/Wtot
# bpe2_perp = bpe2_perp/Wtot

ax.loglog([lmd_MRI, lmd_MRI], [1e-1, 1e+2], 'k-', lw=1)
ax.loglog(kp[1:], u2[1:]*kp[1:]**(+3./2.), lw=3, color=tab10[0], label=r'$E_{u}$')
ax.loglog(kp[1:], b2[1:]*kp[1:]**(+3./2.), lw=3, color=tab10[1], label=r'$E_{B}$')
ax.loglog(kp[1:], kp[1:]**(-5./3.+3./2.)/kp[1]**(-5./3.+3./2.)*b2[1:][0], 'k--', lw=1, label=r'-5/3')


#.............................#
#  Load & plot Sun & Bai data
#.............................#
kp1, \
u2sb = np.loadtxt('u.txt').T
kp2, \
b2sb = np.loadtxt('b.txt').T

kp1 = kp1/kp1[1]
kp2 = kp2/kp2[1]

u2sb = u2sb/kp1
b2sb = b2sb/kp2

norm = u2sb[1]/u2[1]
u2sb = u2sb/norm
b2sb = b2sb/norm

ax.loglog(kp1[1:], u2sb[1:]*kp1[1:]**(+3./2.), '--', lw=3, color=tab10[0], label=r'$E_u$ (SB21)')
ax.loglog(kp2[1:], b2sb[1:]*kp2[1:]**(+3./2.), '--', lw=3, color=tab10[1], label=r'$E_B$ (SB21)')
#  plt.loglog([31, 31], [1e-7, 1e-1], 'k:', lw=2)
			 

ax.set_xlim([1.0, 400.0])
ax.set_ylim([0.5, 30.0])
ax.set_ylabel(r'$k^{3/2}E(k)$')
ax.tick_params(labelbottom=False)
leg = ax.legend(frameon=False, loc='upper right', ncol=2, fontsize=28, handlelength=1.1, columnspacing=0.6)

axes[0][1].axis('off')

#---------------------------------------#
#              Second figure
#---------------------------------------#
ax = axes[1][0]
# load data
kp, \
kpar_b,\
kpar_u,\
b1_ovr_b0 = np.loadtxt('../kpar_avg.txt').T

ax.loglog([lmd_MRI, lmd_MRI], [0.5, 100], 'k-', lw=1)
ax.loglog(kp[1:], kpar_b[1:], lw=3, color=tab10[0], label=r'$k_\|^b$')
ax.loglog(kp[1:], kpar_u[1:], lw=3, color=tab10[1], label=r'$k_\|^u$')
ax.loglog(kp[1:], kp[1:]**(2./3.), 'k--', lw=1, label=r'$k^{2/3}$')
ax.loglog(kp[1:], kp[1:], 'k:', lw=1, label=r'$k^1$')

ax.set_xlim([1.0, 400.0])
ax.set_ylim([0.5, 100.0])
ax.set_ylabel(r'$k_\|$')
ax.tick_params(labelbottom=False)
leg = ax.legend(frameon=False, loc='upper left', ncol=2, fontsize=28, handlelength=1.1, columnspacing=0.6)

#---------------------------------------#
#              Third figure
#---------------------------------------#
ax = axes[1][1]

ax.loglog([lmd_MRI, lmd_MRI], [0.1, 20], 'k-', lw=1)
ax.loglog(kp[1:], b1_ovr_b0[1:], lw=3, color=tab10[2], label=r'$\delta B/B_0$')

ax.set_xlim([1.0, 400.0])
ax.set_ylim([0.1, 20.0])
ax.set_ylabel(r'$\delta B/B_0$')
ax.tick_params(labelbottom=False)
leg = ax.legend(frameon=False, loc='upper right', ncol=2, fontsize=28, handlelength=1.1, columnspacing=0.6)

#---------------------------------------#
#              Fourth figure
#---------------------------------------#
ax = axes[2][0]
# load data
kp, \
b1par2, \
b1prp2, \
u1par2, \
u1prp2 = np.loadtxt('../shearAW_vs_pseudoAW.txt').T

ax.loglog([lmd_MRI, lmd_MRI], [0.5, 30], 'k-', lw=1)
ax.loglog(kp[1:], u1prp2[1:]*kp[1:]**(-1.0 + 3./2.), lw=3, color=tab10[0], label=r'$E_{u_\+}$')
ax.loglog(kp[1:], b1prp2[1:]*kp[1:]**(-1.0 + 3./2.), lw=3, color=tab10[1], label=r'$E_{\delta B_\+}$')
ax.loglog(kp[1:], u1par2[1:]*kp[1:]**(-1.0 + 3./2.), lw=3, color=tab10[2], label=r'$E_{u_\|}$')
ax.loglog(kp[1:], b1par2[1:]*kp[1:]**(-1.0 + 3./2.), lw=3, color=tab10[3], label=r'$E_{\delta B_\|}$')

ax.set_xlim([1.0, 400.0])
ax.set_ylim([0.5, 30.0])
ax.set_ylabel(r'$k^{3/2}E(k)$')
ax.set_xlabel(r'$k$')
leg = ax.legend(frameon=False, loc='lower right', ncol=2, fontsize=28, handlelength=1.1, columnspacing=0.6)

#---------------------------------------#
#              Fifth figure
#---------------------------------------#
ax = axes[2][1]

ax.loglog([lmd_MRI, lmd_MRI], [0.1, 10.0], 'k-', lw=1)
ax.loglog(kp[1:], b1par2[1:]/b1prp2[1:], lw=3, color=tab10[0], label=r'$E_{\delta B_\|}/E_{\delta B_\+}$')
ax.loglog(kp[1:], u1par2[1:]/u1prp2[1:], lw=3, color=tab10[1], label=r'$E_{u_\|}/E_{u_\+}$')

ax.set_xlim([1.0, 400.0])
ax.set_ylim([0.1, 10.0])
ax.set_xlabel(r'$k$')
leg = ax.legend(frameon=False, loc='lower right', fontsize=28, handlelength=1.1, columnspacing=0.6)


plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
plt.savefig('plot.pdf')
