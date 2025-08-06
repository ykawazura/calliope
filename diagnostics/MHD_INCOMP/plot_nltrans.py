# -*- coding: utf-8 -*-
from load import *
import sys
sys.path.append('../')
from plots import *

print('\nplotting shell-to-shell transfer function\n')
outdir = './fig_nltrans/'

k0idx = 3

kpbin_log = kpbin_log_nltrans
kp        = kpbin_log
trans_uu  = trans_uu
trans_bb  = trans_bb
trans_ub  = trans_ub
trans_bu  = trans_bu

#--------------------------------------------------------#
#                     2D transfer map                    #
#--------------------------------------------------------#
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

transes = [trans_uu[final_nltrans_idx,k0idx:,k0idx:], trans_bb[final_nltrans_idx,k0idx:,k0idx:], trans_ub[final_nltrans_idx,k0idx:,k0idx:], trans_bu[final_nltrans_idx,k0idx:,k0idx:]]
names   = [r'$\calT_{uu}$', r'$\calT_{BB}$', r'$\calT_{uB}$', r'$\calT_{Bu}$']

for i, trans in enumerate(transes):
  col = i%2
  row = int(i/2)

  umax = np.max(abs(trans))

  ax = axes[row][col]
  X, Y = np.meshgrid(kp[k0idx:], kp[k0idx:])
  im = ax.pcolormesh(X, Y, trans.T, cmap='RdBu_r', vmin = -umax, vmax = umax)

  ax.set_title(names[i] + r' $ (t = %.2E)$' % tt_nltrans[final_nltrans_idx])

  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cb = plt.colorbar(im, cax=cax)
  cb.formatter.set_powerlimits((0, 0))
  cb.ax.yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
  cb.ax.yaxis.set_offset_position('left')  
  cb.locator = ticker.MaxNLocator(nbins=9)
  cb.update_ticks()

for ax in axes.flatten():
  ax.set_aspect('equal')
  ax.minorticks_on()

  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.set_xlim(kp[k0idx], kp[-1])
  ax.set_ylim(kp[k0idx], kp[-1])

axes[0][0].tick_params(labelbottom=False)
axes[0][1].tick_params(labelbottom=False)
axes[0][1].tick_params(labelleft  =False)
axes[1][1].tick_params(labelleft  =False)

axes[0][0].set_ylabel(r'$K$')
axes[1][0].set_ylabel(r'$K$')
axes[1][0].set_xlabel(r'$Q$')
axes[1][1].set_xlabel(r'$Q$')

plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
plt.savefig(outdir+'trans.pdf')

#------------------#
#   output data    #
#------------------#
from scipy.io import savemat

np.savetxt(outdir + 'tt_nltrans.txt'  , tt_nltrans, fmt='%E')
savemat(outdir + 'grid' , {'kpbin_log':kp})
savemat(outdir + 'trans' , {
                            'tt'         :tt_nltrans[final_nltrans_idx],
                            'nltrans_uu' :trans_uu [final_nltrans_idx,:,:],
                            'nltrans_bb' :trans_bb [final_nltrans_idx,:,:],
                            'nltrans_ub' :trans_ub [final_nltrans_idx,:,:],
                            'nltrans_bu' :trans_bu [final_nltrans_idx,:,:],
                          })


#--------------------------------------------------------#
#                         Fluxes                         #
#--------------------------------------------------------#

flx_uu = np.asarray([np.sum(trans_uu[final_nltrans_idx].T[i:, :i]) for i in np.arange(0, kp.size)])
flx_bb = np.asarray([np.sum(trans_bb[final_nltrans_idx].T[i:, :i]) for i in np.arange(0, kp.size)])
flx_ub = np.asarray([np.sum(trans_ub[final_nltrans_idx].T[i:, :i]) for i in np.arange(0, kp.size)])
flx_bu = np.asarray([np.sum(trans_bu[final_nltrans_idx].T[i:, :i]) for i in np.arange(0, kp.size)])

ys = [ 
      flx_uu[:],
      flx_bb[:],
      flx_ub[:],
      flx_bu[:],
      (flx_uu[:] + flx_bb[:] + flx_bu[:] + flx_ub[:]),
     ]
xs = [ 
      kp[:], 
      kp[:], 
      kp[:], 
      kp[:], 
      kp[:], 
     ]
ls = [ 
        '', 
        '', 
        '', 
        '', 
        '', 
     ]
legends = [ 
            r'$\Pi^{u^<}_{u^>}$', 
            r'$\Pi^{B^<}_{B^>}$', 
            r'$\Pi^{u^<}_{B^>}$', 
            r'$\Pi^{B^<}_{u^>}$', 
            r'$\Pi_\mr{tot}$', 
          ]

plot_semilogx1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt_nltrans[final_nltrans_idx], ylab='', term=True, save=outdir+'flux.pdf')
