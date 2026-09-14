# -*- coding: utf-8 -*-
from load import *
from fft import *
import sys
import sys
sys.path.append('../')
sys.path.append('../')
from plots import *

print('\nplotting kspectrum\n')
outdir = './fig_kspectrum/'

upe2_bin    = sum_negative_kz2d(upe2_bin)
bpe2_bin    = sum_negative_kz2d(bpe2_bin)
upa2_bin    = sum_negative_kz2d(upa2_bin)
bpa2_bin    = sum_negative_kz2d(bpa2_bin)
ux2_bin     = sum_negative_kz2d(ux2_bin)
uy2_bin     = sum_negative_kz2d(uy2_bin)
bx2_bin     = sum_negative_kz2d(bx2_bin)
by2_bin     = sum_negative_kz2d(by2_bin)
zppe2_bin   = sum_negative_kz2d(zppe2_bin)
zmpe2_bin   = sum_negative_kz2d(zmpe2_bin)
zppa2_bin   = sum_negative_kz2d(zppa2_bin)
zmpa2_bin   = sum_negative_kz2d(zmpa2_bin)
p_aw_bin    = sum_negative_kz2d(p_aw_bin)
p_compr_bin = sum_negative_kz2d(p_compr_bin)
dissip_aw_bin    = sum_negative_kz2d(dissip_aw_bin)
dissip_compr_bin = sum_negative_kz2d(dissip_compr_bin)
dissip_KAW_bin   = sum_negative_kz2d(dissip_KAW_bin)
dissip_ICW_bin   = sum_negative_kz2d(dissip_ICW_bin)
ntrans_upe_upe_l_bin = sum_negative_kz2d(ntrans_upe_upe_l_bin)
ntrans_bpe_upe_l_bin = sum_negative_kz2d(ntrans_bpe_upe_l_bin)
ntrans_bpe_bpe_l_bin = sum_negative_kz2d(ntrans_bpe_bpe_l_bin)
ntrans_upe_bpe_l_bin = sum_negative_kz2d(ntrans_upe_bpe_l_bin)
ntrans_upa_upa_l_bin = sum_negative_kz2d(ntrans_upa_upa_l_bin)
ntrans_bpa_upa_l_bin = sum_negative_kz2d(ntrans_bpa_upa_l_bin)
ntrans_bpa_bpa_l_bin = sum_negative_kz2d(ntrans_bpa_bpa_l_bin)
ntrans_upa_bpa_l_bin = sum_negative_kz2d(ntrans_upa_bpa_l_bin)
ntrans_upe_upe_g_bin = sum_negative_kz2d(ntrans_upe_upe_g_bin)
ntrans_bpe_upe_g_bin = sum_negative_kz2d(ntrans_bpe_upe_g_bin)
ntrans_bpe_bpe_g_bin = sum_negative_kz2d(ntrans_bpe_bpe_g_bin)
ntrans_upe_bpe_g_bin = sum_negative_kz2d(ntrans_upe_bpe_g_bin)
ntrans_upa_upa_g_bin = sum_negative_kz2d(ntrans_upa_upa_g_bin)
ntrans_bpa_upa_g_bin = sum_negative_kz2d(ntrans_bpa_upa_g_bin)
ntrans_bpa_bpa_g_bin = sum_negative_kz2d(ntrans_bpa_bpa_g_bin)
ntrans_upa_bpa_g_bin = sum_negative_kz2d(ntrans_upa_bpa_g_bin)

ntrans_aw_l_bin    = ntrans_upe_upe_l_bin + ntrans_bpe_upe_l_bin + ntrans_bpe_bpe_l_bin + ntrans_upe_bpe_l_bin
ntrans_aw_g_bin    = ntrans_upe_upe_g_bin + ntrans_bpe_upe_g_bin + ntrans_bpe_bpe_g_bin + ntrans_upe_bpe_g_bin
ntrans_compr_l_bin = ntrans_upa_upa_l_bin + ntrans_bpa_upa_l_bin + ntrans_bpa_bpa_l_bin + ntrans_upa_bpa_l_bin
ntrans_compr_g_bin = ntrans_upa_upa_g_bin + ntrans_bpa_upa_g_bin + ntrans_bpa_bpa_g_bin + ntrans_upa_bpa_g_bin

upe2_KAW_bin    = sum_negative_kz2d(upe2_KAW_bin)
bpe2_KAW_bin    = sum_negative_kz2d(bpe2_KAW_bin)
upa2_KAW_bin    = sum_negative_kz2d(upa2_KAW_bin)
bpa2_KAW_bin    = sum_negative_kz2d(bpa2_KAW_bin)
upe2_ICW_bin    = sum_negative_kz2d(upe2_ICW_bin)
bpe2_ICW_bin    = sum_negative_kz2d(bpe2_ICW_bin)
upa2_ICW_bin    = sum_negative_kz2d(upa2_ICW_bin)
bpa2_ICW_bin    = sum_negative_kz2d(bpa2_ICW_bin)

# Nonlinear fluxes of the free energy W and the generalized helicity X
# (cf. Meyrand et al. 2021, Fig. 6). Runs made before these were written to
# calliope.out.nc do not have them, hence the has_flux guard everywhere below.
#
# Pi(kprp) is the flux through a kprp shell, so it has to be summed over the
# FULL kz axis. sum_negative_kz2d folds +-kz with a 1/2 weight, which keeps the
# shape but halves the total, so the kz sum is taken here, before folding. The
# folded arrays are kept as well, for the combined figure below where every
# other curve is folded the same way.
if has_flux:
  flux_W_kprp = np.sum(flux_W_bin, axis = 1)                       # [t, kprp]
  flux_X_kprp = np.sum(flux_X_bin, axis = 1)                       # [t, kprp]

  # injection rates used to normalize Pi; p_*_sum are the true totals
  eps_W = np.mean(p_aw_sum + p_compr_sum)
  eps_X = np.mean(p_xhl_sum)

  flux_W_bin  = sum_negative_kz2d(flux_W_bin)
  flux_X_bin  = sum_negative_kz2d(flux_X_bin)

if nlz == nkz:
  kp_end = np.argmin(np.abs(kpbin - kpbin.max()*2./3.))
  if not is2D:
    kz_end = np.argmin(np.abs(kz[1:int(nkz/2)] - kz[1:int(nkz/2)].max()*2./3.))
else:
  kp_end = kpbin.size - 1
  kz_end = int(nkz/2)

#--------------------------------------------------------#
#                      plot 1D spectra                   #
#--------------------------------------------------------#
# kprp spectrum
ys = [ 
       np.sum(upe2_bin[final_idx, :, 1:kp_end], axis=0), 
       np.sum(bpe2_bin[final_idx, :, 1:kp_end], axis=0), 
       np.sum(upa2_bin[final_idx, :, 1:kp_end], axis=0), 
       np.sum(bpa2_bin[final_idx, :, 1:kp_end], axis=0),
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*np.sum(bpe2_bin[final_idx,:,1:kp_end], axis=0)[0]
	 ]
xs = [ 
        kpbin[1:kp_end], 
        kpbin[1:kp_end], 
        kpbin[1:kp_end], 
        kpbin[1:kp_end], 
        kpbin[1:kp_end], 
     ]
ls = [ 
        '', 
        '', 
        '', 
        '', 
        'k--', 
     ]
legends = [ 
            r'$E_{u_\+}$', 
            r'$E_{\delta B_\+}$',
            r'$E_{u_\|}$', 
            r'$E_{\delta B_\|}$',
            r'-5/3',
         ]
plot_log1d_many(xs, ys, xlab=r'$'+kplab+'$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', vline_at = 1./rho, term=True, save=outdir+'kprp_spectra.pdf')

# kprp spectrum by components
ys = [ 
       np.sum(upe2_bin[final_idx, :, 1:kp_end], axis=0), 
       np.sum(ux2_bin [final_idx, :, 1:kp_end], axis=0), 
       np.sum(uy2_bin [final_idx, :, 1:kp_end], axis=0), 
       np.sum(bpe2_bin[final_idx, :, 1:kp_end], axis=0), 
       np.sum(bx2_bin [final_idx, :, 1:kp_end], axis=0), 
       np.sum(by2_bin [final_idx, :, 1:kp_end], axis=0), 
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*np.sum(bpe2_bin[final_idx,:,1:kp_end], axis=0)[0]
     ]
xs = [ 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end]  
     ]
ls = [ 
        'r-', 
        'r--', 
        'r:', 
        'b-', 
        'b--', 
        'b:', 
        'k--', 
     ]
legends = [ 
            r'$E_{u_\+}$', 
            r'$E_{u_x}$', 
            r'$E_{u_y}$', 
            r'$E_{\delta B_\+}$',
            r'$E_{\delta B_x}$',
            r'$E_{\delta B_y}$',
            r'-5/3',
          ]
plot_log1d_many(xs, ys, xlab=r'$'+kplab+'$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', vline_at = 1./rho, term=True, save=outdir+'kprp_spectra_components.pdf')

# kprp spectrum by MRI injection rate and nonlinear transfer rate
ys = [ 
       np.sum(p_aw_bin          [final_idx,:,1:kp_end], axis=0),
       np.sum(p_compr_bin       [final_idx,:,1:kp_end], axis=0),
       np.sum(dissip_aw_bin     [final_idx,:,1:kp_end], axis=0),
       np.sum(dissip_compr_bin  [final_idx,:,1:kp_end], axis=0),
       np.sum(ntrans_aw_l_bin   [final_idx,:,1:kp_end], axis=0),
       np.sum(ntrans_aw_g_bin   [final_idx,:,1:kp_end], axis=0),
       np.sum(ntrans_compr_l_bin[final_idx,:,1:kp_end], axis=0),
       np.sum(ntrans_compr_g_bin[final_idx,:,1:kp_end], axis=0),
     ]
xs = [ 
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
      kpbin[1:kp_end]  
     ]
ls = [ 
        '', 
        '', 
        '', 
        '', 
        '', 
        '', 
        '', 
        '', 
     ]
legends = [ 
            r'$I_\mr{AW}$', 
            r'$I_\mr{compr}$', 
            r'$\calD_\mr{AW}$', 
            r'$\calD_\mr{compr}$', 
            r'$\calN_\mr{AW}^{<k_\+}$', 
            r'$\calN_\mr{AW}^{>k_\+}$', 
            r'$\calN_\mr{compr}^{<k_\+}$', 
            r'$\calN_\mr{compr}^{>k_\+}$', 
          ]
if has_flux:
  ys      += [
               np.sum(flux_W_bin[final_idx,:,1:kp_end], axis=0),
               np.sum(flux_X_bin[final_idx,:,1:kp_end], axis=0),
             ]
  xs      += [kpbin[1:kp_end], kpbin[1:kp_end]]
  # black, because the 8-colour prop_cycle has already wrapped around by here
  ls      += ['k-', 'k--']
  legends += [r'$\Pi_W$', r'$\Pi_X$']
plot_log1d_many(xs, ys, xlab=r'$'+kplab+'$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', vline_at = 1./rho, term=True, save=outdir+'kprp_spectra_flux.pdf')

# nonlinear flux of W through the kprp shells, one curve per output time
# (cf. Meyrand et al. 2021, Fig. 6)
if has_flux:
  plot_semilogx1d_time(kpbin[1:kp_end], flux_W_kprp[:, 1:kp_end]/eps_W, tt,
      xlab=r'$'+kplab+'$', clab=tlab, ylab=r'$\Pi_W/\varepsilon_W$',
      title=r'$t = %.2E - %.2E$' % (tt[0], tt[final_idx]),
      vline_at=1./rho, hline_at=(0., 1.), hline_ls=(':', '--'), mean_line=True,
      term=True, save=outdir+'kprp_flux_W.pdf')

  # same for the generalized helicity. Balanced runs (xhl_inj = 0) only have a
  # residual eps_X left over from the random forcing, so normalizing by it is
  # meaningless; the raw flux is plotted instead.
  if np.abs(eps_X) > 1.e-3*np.abs(eps_W):
    ys_X   = flux_X_kprp[:, 1:kp_end]/eps_X
    ylab_X = r'$\Pi_X/\varepsilon_X$'
    hl_X, hls_X = (0., 1.), (':', '--')
  else:
    ys_X   = flux_X_kprp[:, 1:kp_end]
    ylab_X = r'$\Pi_X$'
    hl_X, hls_X = 0., None          # the raw flux has no natural unit level
  plot_semilogx1d_time(kpbin[1:kp_end], ys_X, tt,
      xlab=r'$'+kplab+'$', clab=tlab, ylab=ylab_X,
      title=r'$t = %.2E - %.2E$' % (tt[0], tt[final_idx]),
      vline_at=1./rho, hline_at=hl_X, hline_ls=hls_X, mean_line=True,
      term=True, save=outdir+'kprp_flux_X.pdf')

# Elsasser fields
ys = [ 
       np.sum(zppe2_bin   [final_idx, :, 1:kp_end], axis=0), 
       np.sum(zmpe2_bin   [final_idx, :, 1:kp_end], axis=0), 
       np.sum(zppa2_bin   [final_idx, :, 1:kp_end], axis=0), 
       np.sum(zmpa2_bin   [final_idx, :, 1:kp_end], axis=0),
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*np.sum(zppe2_bin[final_idx,:,1:kp_end], axis=0)[0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*np.sum(zppe2_bin[final_idx,:,1:kp_end], axis=0)[0]
     ]
xs = [ 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end], 
      kpbin[1:kp_end]  
     ]
ls = [ 
        '', 
        '', 
        '', 
        '', 
        'k--', 
        'k--', 
     ]
legends = [ 
            r'$E_{Z^+_\+}$', 
            r'$E_{Z^-_\+}$', 
            r'$E_{Z^+_\|}$', 
            r'$E_{Z^-_\|}$', 
            r'-5/3',
            r'-3/2',
          ]
plot_log1d_many(xs, ys, xlab=r'$'+kplab+'$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', vline_at = 1./rho, term=True, save=outdir+'kprp_spectra_ELS.pdf')

# KAW and ICW decomposition
ys = [ 
       np.sum(upe2_KAW_bin[final_idx, :, 1:kp_end], axis=0), 
       np.sum(bpe2_KAW_bin[final_idx, :, 1:kp_end], axis=0), 
       np.sum(upa2_KAW_bin[final_idx, :, 1:kp_end], axis=0), 
       np.sum(bpa2_KAW_bin[final_idx, :, 1:kp_end], axis=0),
       np.sum(upe2_ICW_bin[final_idx, :, 1:kp_end], axis=0), 
       np.sum(bpe2_ICW_bin[final_idx, :, 1:kp_end], axis=0), 
       np.sum(upa2_ICW_bin[final_idx, :, 1:kp_end], axis=0), 
       np.sum(bpa2_ICW_bin[final_idx, :, 1:kp_end], axis=0),
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*np.sum(bpe2_KAW_bin[final_idx,:,1:kp_end], axis=0)[0]
	 ]
xs = [ 
        kpbin[1:kp_end], 
        kpbin[1:kp_end], 
        kpbin[1:kp_end], 
        kpbin[1:kp_end], 
        kpbin[1:kp_end], 
        kpbin[1:kp_end], 
        kpbin[1:kp_end], 
        kpbin[1:kp_end], 
        kpbin[1:kp_end], 
     ]
ls = [ 
        '', 
        '', 
        '', 
        '', 
        '', 
        '', 
        '', 
        '', 
        'k--', 
     ]
legends = [ 
            r'$E_{\widetilde{u}_\+}$', 
            r'$E_{\widetilde{\delta B}_\+}$',
            r'$E_{\widetilde{u}_\|}$', 
            r'$E_{\widetilde{\delta B}_\|}$',
            r'$E_{\overline{u}_\+}$', 
            r'$E_{\overline{\delta B}_\+}$',
            r'$E_{\overline{u}_\|}$', 
            r'$E_{\overline{\delta B}_\|}$',
            r'-5/3',
         ]
plot_log1d_many(xs, ys, xlab=r'$'+kplab+'$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', vline_at = 1./rho, term=True, save=outdir+'kprp_KAW_ICW_spectra.pdf')

# kz spectrum
if not is2D:
  ys = [ 
         np.sum(upe2_bin[final_idx, 1:kz_end, :kp_end], axis=1), 
         np.sum(bpe2_bin[final_idx, 1:kz_end, :kp_end], axis=1), 
         np.sum(upa2_bin[final_idx, 1:kz_end, :kp_end], axis=1), 
         np.sum(bpa2_bin[final_idx, 1:kz_end, :kp_end], axis=1), 
       ]
  xs = [ 
          kz[1:kz_end], 
          kz[1:kz_end], 
          kz[1:kz_end], 
          kz[1:kz_end], 
       ]
  ls = [ 
          '', 
          '', 
          '', 
          '', 
       ]
  legends = [ 
              r'$E_{u_\+}$', 
              r'$E_{\delta B_\+}$',
              r'$E_{u_\|}$', 
              r'$E_{\delta B_\|}$',
            ]
  plot_log1d_many(xs, ys, xlab='$'+kzlab+'$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kz_spectra.pdf')

  # Elsasser fields
  ys = [ 
         np.sum(zppe2_bin[final_idx, 1:kz_end, :kp_end], axis=1), 
         np.sum(zmpe2_bin[final_idx, 1:kz_end, :kp_end], axis=1), 
         np.sum(zppa2_bin[final_idx, 1:kz_end, :kp_end], axis=1), 
         np.sum(zmpa2_bin[final_idx, 1:kz_end, :kp_end], axis=1), 
       ]
  xs = [ 
          kz[1:kz_end], 
          kz[1:kz_end], 
          kz[1:kz_end], 
          kz[1:kz_end], 
       ]
  ls = [ 
          '', 
          '', 
          '', 
          '', 
       ]
  legends = [ 
              r'$E_{Z^+_\+}$', 
              r'$E_{Z^-_\+}$', 
              r'$E_{Z^+_\|}$', 
              r'$E_{Z^-_\|}$', 
            ]
  plot_log1d_many(xs, ys, xlab='$'+kzlab+'$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kz_spectra_ELS.pdf')

  # MRI injection rate and nonlinear transfer rate
  ys = [ 
          np.sum(p_aw_bin          [final_idx,1:kz_end,:kp_end], axis=1),
          np.sum(p_compr_bin       [final_idx,1:kz_end,:kp_end], axis=1),
          np.sum(dissip_aw_bin     [final_idx,1:kz_end,:kp_end], axis=1),
          np.sum(dissip_compr_bin  [final_idx,1:kz_end,:kp_end], axis=1),
          np.sum(ntrans_aw_l_bin   [final_idx,1:kz_end,:kp_end], axis=1),
          np.sum(ntrans_aw_g_bin   [final_idx,1:kz_end,:kp_end], axis=1),
          np.sum(ntrans_compr_l_bin[final_idx,1:kz_end,:kp_end], axis=1),
          np.sum(ntrans_compr_g_bin[final_idx,1:kz_end,:kp_end], axis=1),
       ]
  xs = [ 
          kz[1:kz_end], 
          kz[1:kz_end], 
          kz[1:kz_end], 
          kz[1:kz_end], 
          kz[1:kz_end], 
          kz[1:kz_end], 
          kz[1:kz_end], 
          kz[1:kz_end], 
       ]
  ls = [ 
          '', 
          '', 
          '', 
          '', 
          '', 
          '', 
          '', 
          '', 
       ]
  legends = [ 
              r'$I_\mr{AW}$', 
              r'$I_\mr{compr}$', 
              r'$\calD_\mr{AW}$', 
              r'$\calD_\mr{compr}$', 
              r'$\calN_\mr{AW}^{<k_\+}$', 
              r'$\calN_\mr{AW}^{>k_\+}$', 
              r'$\calN_\mr{compr}^{<k_\+}$', 
              r'$\calN_\mr{compr}^{>k_\+}$', 
            ]
  plot_log1d_many(xs, ys, xlab='$'+kzlab+'$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kz_spectra_flux.pdf')

#--------------------------------------------------------#
#                      plot 2D spectra                   #
#--------------------------------------------------------#
if not is2D:
  plot_log2d(upe2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$E_{u_{\+}}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'upe2.pdf')
  plot_log2d(bpe2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$E_{\delta B_\+}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'bpe2.pdf')
  plot_log2d(upa2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$E_{u_{\|}}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'upa2.pdf')
  plot_log2d(bpa2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$E_{\delta B_\|}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'bpa2.pdf')
  plot_log2d(zppe2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$E_{Z_{\+}^+}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'zppe2.pdf')
  plot_log2d(zmpe2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$E_{Z_{\+}^-}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'zmpe2.pdf')
  plot_log2d(zppa2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$E_{Z_{\|}^+}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'zppa2.pdf')
  plot_log2d(zmpa2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$E_{Z_{\|}^-}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'zmpa2.pdf')
  plot_log2d(dissip_KAW_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab=r'$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$\calD_\mr{KAW}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'dissip_KAW.pdf')
  plot_log2d(dissip_ICW_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab=r'$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$\calD_\mr{ICW}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'dissip_ICW.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ekprp.txt'  , np.column_stack((kpbin[:kp_end], 
                                                       np.sum(upe2_bin            [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(bpe2_bin            [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(upa2_bin            [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(bpa2_bin            [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ux2_bin             [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(uy2_bin             [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(bx2_bin             [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(by2_bin             [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(zppe2_bin           [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(zmpe2_bin           [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(zppa2_bin           [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(zmpa2_bin           [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(p_aw_bin            [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(p_compr_bin         [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(dissip_aw_bin       [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(dissip_compr_bin    [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(dissip_KAW_bin      [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(dissip_ICW_bin      [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_upe_upe_l_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_bpe_upe_l_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_bpe_bpe_l_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_upe_bpe_l_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_upa_upa_l_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_bpa_upa_l_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_bpa_bpa_l_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_upa_bpa_l_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_upe_upe_g_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_bpe_upe_g_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_bpe_bpe_g_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_upe_bpe_g_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_upa_upa_g_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_bpa_upa_g_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_bpa_bpa_g_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_upa_bpa_g_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(upe2_KAW_bin        [final_idx,:kz_end,:kp_end], axis=0), 
                                                       np.sum(bpe2_KAW_bin        [final_idx,:kz_end,:kp_end], axis=0), 
                                                       np.sum(upa2_KAW_bin        [final_idx,:kz_end,:kp_end], axis=0), 
                                                       np.sum(bpa2_KAW_bin        [final_idx,:kz_end,:kp_end], axis=0),
                                                       np.sum(upe2_ICW_bin        [final_idx,:kz_end,:kp_end], axis=0), 
                                                       np.sum(bpe2_ICW_bin        [final_idx,:kz_end,:kp_end], axis=0), 
                                                       np.sum(upa2_ICW_bin        [final_idx,:kz_end,:kp_end], axis=0), 
                                                       np.sum(bpa2_ICW_bin        [final_idx,:kz_end,:kp_end], axis=0),
                                                     )), fmt='%E')
if not is2D:
  np.savetxt(outdir + 'Ekz.txt'  , np.column_stack((kz[:kz_end], 
                                                         np.sum(upe2_bin            [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(bpe2_bin            [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(upa2_bin            [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(bpa2_bin            [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ux2_bin             [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(uy2_bin             [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(bx2_bin             [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(by2_bin             [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(zppe2_bin           [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(zmpe2_bin           [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(zppa2_bin           [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(zmpa2_bin           [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(p_aw_bin            [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(p_compr_bin         [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(dissip_aw_bin       [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(dissip_compr_bin    [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(dissip_KAW_bin      [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(dissip_ICW_bin      [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_upe_upe_l_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_bpe_upe_l_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_bpe_bpe_l_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_upe_bpe_l_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_upa_upa_l_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_bpa_upa_l_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_bpa_bpa_l_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_upa_bpa_l_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_upe_upe_g_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_bpe_upe_g_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_bpe_bpe_g_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_upe_bpe_g_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_upa_upa_g_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_bpa_upa_g_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_bpa_bpa_g_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_upa_bpa_g_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(upe2_KAW_bin        [final_idx,:kz_end,:kp_end], axis=1), 
                                                         np.sum(bpe2_KAW_bin        [final_idx,:kz_end,:kp_end], axis=1), 
                                                         np.sum(upa2_KAW_bin        [final_idx,:kz_end,:kp_end], axis=1), 
                                                         np.sum(bpa2_KAW_bin        [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(upe2_ICW_bin        [final_idx,:kz_end,:kp_end], axis=1), 
                                                         np.sum(bpe2_ICW_bin        [final_idx,:kz_end,:kp_end], axis=1), 
                                                         np.sum(upa2_ICW_bin        [final_idx,:kz_end,:kp_end], axis=1), 
                                                         np.sum(bpa2_ICW_bin        [final_idx,:kz_end,:kp_end], axis=1),
                                                       )), fmt='%E')

del upe2_bin
del bpe2_bin
del upa2_bin
del bpa2_bin
