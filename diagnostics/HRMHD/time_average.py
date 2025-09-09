# -*- coding: utf-8 -*-
from load import *
from fft import *
import sys
sys.path.append('../')
from plots import *
from numpy import trapezoid
from scipy import interpolate

avg_start = 500
avg_end   = -1

def time_average(x, y, axis=0): # x: 1D array, y: any-D array 
  return trapezoid(y, x, axis=axis)/(x[-1] - x[0])
  # return trapz(y, x, axis=axis)/(x[-1] - x[0])

def std_dev(x, y): # x, y: 1D array
  y_intp = interpolate.interp1d(x, y)
  return np.std(y_intp(np.linspace(x[0], x[-1], x.size)), ddof=0)

def cov(x, y1, y2): # x, y: 1D array
  y1_intp  = interpolate.interp1d(x, y1)
  y2_intp  = interpolate.interp1d(x, y2)
  y1_intp_ = y1_intp(np.linspace(x[0], x[-1], x.size))
  y2_intp_ = y2_intp(np.linspace(x[0], x[-1], x.size))
  
  return np.cov( np.stack((y1_intp_, y2_intp_), axis=0) )

if avg_end == -1:
	avg_end = np.argwhere(tt == tt[avg_end])[0][0]

##########################################################
#              average energy time evolution             #
##########################################################
print('\nplotting energy\n')
outdir = './fig_energy/'

Waw        = upe2_sum + bpe2_sum
Wcompr     = upa2_sum + bpa2_sum
Wzpe       = zppe2_sum + zmpe2_sum
Wzpa       = zppa2_sum + zmpa2_sum
Hzpe       = zppe2_sum - zmpe2_sum
Hzpa       = zppa2_sum - zmpa2_sum
Wtot       = Waw + Wcompr
Waw_dot    = upe2dot_sum + bpe2dot_sum
Wcompr_dot = upa2dot_sum + bpa2dot_sum
Wtot_dot   = Waw_dot + Wcompr_dot
Daw        = upe2dissip_sum + bpe2dissip_sum
Dcompr     = upa2dissip_sum + bpa2dissip_sum
Dtot       = Daw + Dcompr
Paw        = p_aw_sum
Pcompr     = p_compr_sum
Ptot       = Paw + Pcompr
D_KAW      = upe2_KAW_dissip_sum + bpe2_KAW_dissip_sum + upa2_KAW_dissip_sum + bpa2_KAW_dissip_sum
D_ICW      = upe2_ICW_dissip_sum + bpe2_ICW_dissip_sum + upa2_ICW_dissip_sum + bpa2_ICW_dissip_sum
D_KAW_prp  = KAW_dissip_prp_sum
D_KAW_par  = KAW_dissip_par_sum
D_ICW_prp  = ICW_dissip_prp_sum
D_ICW_par  = ICW_dissip_par_sum

Waw_avg        = time_average(tt[avg_start:avg_end], Waw       [avg_start:avg_end], axis=0)
Wcompr_avg     = time_average(tt[avg_start:avg_end], Wcompr    [avg_start:avg_end], axis=0)
Wzpe_avg       = time_average(tt[avg_start:avg_end], Wzpe      [avg_start:avg_end], axis=0)
Wzpa_avg       = time_average(tt[avg_start:avg_end], Wzpa      [avg_start:avg_end], axis=0)
Hzpe_avg       = time_average(tt[avg_start:avg_end], Hzpe      [avg_start:avg_end], axis=0)
Hzpa_avg       = time_average(tt[avg_start:avg_end], Hzpa      [avg_start:avg_end], axis=0)
Wtot_avg       = time_average(tt[avg_start:avg_end], Wtot      [avg_start:avg_end], axis=0)
Waw_dot_avg    = time_average(tt[avg_start:avg_end], Waw_dot   [avg_start:avg_end], axis=0)
Wcompr_dot_avg = time_average(tt[avg_start:avg_end], Wcompr_dot[avg_start:avg_end], axis=0)
Wtot_dot_avg   = time_average(tt[avg_start:avg_end], Wtot_dot  [avg_start:avg_end], axis=0)
Daw_avg        = time_average(tt[avg_start:avg_end], Daw       [avg_start:avg_end], axis=0)
Dcompr_avg     = time_average(tt[avg_start:avg_end], Dcompr    [avg_start:avg_end], axis=0)
Dtot_avg       = time_average(tt[avg_start:avg_end], Dtot      [avg_start:avg_end], axis=0)
Paw_avg        = time_average(tt[avg_start:avg_end], Paw       [avg_start:avg_end], axis=0)
Pcompr_avg     = time_average(tt[avg_start:avg_end], Pcompr    [avg_start:avg_end], axis=0)
Ptot_avg       = time_average(tt[avg_start:avg_end], Ptot      [avg_start:avg_end], axis=0)
D_KAW_avg      = time_average(tt[avg_start:avg_end], D_KAW     [avg_start:avg_end], axis=0)
D_ICW_avg      = time_average(tt[avg_start:avg_end], D_ICW     [avg_start:avg_end], axis=0)
D_KAW_prp_avg  = time_average(tt[avg_start:avg_end], D_KAW_prp [avg_start:avg_end], axis=0)
D_KAW_par_avg  = time_average(tt[avg_start:avg_end], D_KAW_par [avg_start:avg_end], axis=0)
D_ICW_prp_avg  = time_average(tt[avg_start:avg_end], D_ICW_prp [avg_start:avg_end], axis=0)
D_ICW_par_avg  = time_average(tt[avg_start:avg_end], D_ICW_par [avg_start:avg_end], axis=0)
upe2_sum_avg   = time_average(tt[avg_start:avg_end], upe2_sum  [avg_start:avg_end], axis=0)
bpe2_sum_avg   = time_average(tt[avg_start:avg_end], bpe2_sum  [avg_start:avg_end], axis=0)
upa2_sum_avg   = time_average(tt[avg_start:avg_end], upa2_sum  [avg_start:avg_end], axis=0)
bpa2_sum_avg   = time_average(tt[avg_start:avg_end], bpa2_sum  [avg_start:avg_end], axis=0)
zppe2_sum_avg  = time_average(tt[avg_start:avg_end], zppe2_sum [avg_start:avg_end], axis=0)
zmpe2_sum_avg  = time_average(tt[avg_start:avg_end], zmpe2_sum [avg_start:avg_end], axis=0)
zppa2_sum_avg  = time_average(tt[avg_start:avg_end], zppa2_sum [avg_start:avg_end], axis=0)
zmpa2_sum_avg  = time_average(tt[avg_start:avg_end], zmpa2_sum [avg_start:avg_end], axis=0)

Waw_err        = std_dev(tt[avg_start:avg_end], Waw       [avg_start:avg_end])
Wcompr_err     = std_dev(tt[avg_start:avg_end], Wcompr    [avg_start:avg_end])
Wzpe_err       = std_dev(tt[avg_start:avg_end], Wzpe      [avg_start:avg_end])
Wzpa_err       = std_dev(tt[avg_start:avg_end], Wzpa      [avg_start:avg_end])
Hzpe_err       = std_dev(tt[avg_start:avg_end], Hzpe      [avg_start:avg_end])
Hzpa_err       = std_dev(tt[avg_start:avg_end], Hzpa      [avg_start:avg_end])
Wtot_err       = std_dev(tt[avg_start:avg_end], Wtot      [avg_start:avg_end])
Waw_dot_err    = std_dev(tt[avg_start:avg_end], Waw_dot   [avg_start:avg_end])
Wcompr_dot_err = std_dev(tt[avg_start:avg_end], Wcompr_dot[avg_start:avg_end])
Wtot_dot_err   = std_dev(tt[avg_start:avg_end], Wtot_dot  [avg_start:avg_end])
Daw_err        = std_dev(tt[avg_start:avg_end], Daw       [avg_start:avg_end])
Dcompr_err     = std_dev(tt[avg_start:avg_end], Dcompr    [avg_start:avg_end])
Dtot_err       = std_dev(tt[avg_start:avg_end], Dtot      [avg_start:avg_end])
Paw_err        = std_dev(tt[avg_start:avg_end], Paw       [avg_start:avg_end])
Pcompr_err     = std_dev(tt[avg_start:avg_end], Pcompr    [avg_start:avg_end])
Ptot_err       = std_dev(tt[avg_start:avg_end], Ptot      [avg_start:avg_end])
D_KAW_err      = std_dev(tt[avg_start:avg_end], D_KAW     [avg_start:avg_end])
D_ICW_err      = std_dev(tt[avg_start:avg_end], D_ICW     [avg_start:avg_end])
upe2_sum_err   = std_dev(tt[avg_start:avg_end], upe2_sum  [avg_start:avg_end])
bpe2_sum_err   = std_dev(tt[avg_start:avg_end], bpe2_sum  [avg_start:avg_end])
upa2_sum_err   = std_dev(tt[avg_start:avg_end], upa2_sum  [avg_start:avg_end])
bpa2_sum_err   = std_dev(tt[avg_start:avg_end], bpa2_sum  [avg_start:avg_end])
zppe2_sum_err  = std_dev(tt[avg_start:avg_end], zppe2_sum [avg_start:avg_end])
zmpe2_sum_err  = std_dev(tt[avg_start:avg_end], zmpe2_sum [avg_start:avg_end])
zppa2_sum_err  = std_dev(tt[avg_start:avg_end], zppa2_sum [avg_start:avg_end])
zmpa2_sum_err  = std_dev(tt[avg_start:avg_end], zmpa2_sum [avg_start:avg_end])
D_KAW_prp_err  = std_dev(tt[avg_start:avg_end], D_KAW_prp [avg_start:avg_end])
D_KAW_par_err  = std_dev(tt[avg_start:avg_end], D_KAW_par [avg_start:avg_end])
D_ICW_prp_err  = std_dev(tt[avg_start:avg_end], D_ICW_prp [avg_start:avg_end])
D_ICW_par_err  = std_dev(tt[avg_start:avg_end], D_ICW_par [avg_start:avg_end])

s =     r'average over t     = [%.3E' % tt[avg_start] + ', %.3E' % tt[avg_end] + ']' + '\n'
s = s + r'average over index = [' + str(avg_start) + ', ' + str(avg_end) + ']' + '\n'
s = s + r' error of ratio is calculated by (a + da)/(b + db) ~ a/b*[1 + sqrt( (da/a)^2 + (db/b)^2 - 2cov/(a*b))]' + '\n\n' 
s = s + r'  Wtot             = %.3E \pm %.3E'  % (Wtot_avg      , Wtot_err      ) + '\n'
s = s + r'    (Waw           = %.3E \pm %.3E)' % (Waw_avg       , Waw_err       ) + '\n'
s = s + r'    (Wcompr        = %.3E \pm %.3E)' % (Wcompr_avg    , Wcompr_err    ) + '\n'
s = s + r'    (Wzpe          = %.3E \pm %.3E)' % (Wzpe_avg      , Wzpe_err      ) + '\n'
s = s + r'    (Wzpa          = %.3E \pm %.3E)' % (Wzpa_avg      , Wzpa_err      ) + '\n'
s = s + r'  Wtot_dot         = %.3E \pm %.3E'  % (Wtot_dot_avg  , Wtot_dot_err  ) + '\n'
s = s + r'    (Waw_dot       = %.3E \pm %.3E)' % (Waw_dot_avg   , Waw_dot_err   ) + '\n'
s = s + r'    (Wcompr_dot    = %.3E \pm %.3E)' % (Wcompr_dot_avg, Wcompr_dot_err) + '\n'
s = s + r'  Ptot             = %.3E \pm %.3E'  % (Ptot_avg      , Ptot_err      ) + '\n'
s = s + r'    Paw            = %.3E \pm %.3E'  % (Paw_avg       , Paw_err       ) + '\n'
s = s + r'    Paw            = %.3E \pm %.3E'  % (Paw_avg       , Paw_err       ) + '\n'
s = s + r'  Dtot             = %.3E \pm %.3E'  % (Dtot_avg      , Dtot_err      ) + '\n'
s = s + r'    Daw            = %.3E \pm %.3E'  % (Daw_avg       , Daw_err       ) + '\n'
s = s + r'    Dcompr         = %.3E \pm %.3E'  % (Dcompr_avg    , Dcompr_err    ) + '\n'
s = s + r'    D_KAW          = %.3E \pm %.3E'  % (D_KAW_avg     , D_KAW_err     ) + '\n'
s = s + r'      D_KAW_prp    = %.3E \pm %.3E'  % (D_KAW_prp_avg , D_KAW_prp_err ) + '\n'
s = s + r'      D_KAW_par    = %.3E \pm %.3E'  % (D_KAW_par_avg , D_KAW_par_err ) + '\n'
s = s + r'    D_ICW          = %.3E \pm %.3E'  % (D_ICW_avg     , D_ICW_err     ) + '\n'
s = s + r'      D_ICW_prp    = %.3E \pm %.3E'  % (D_ICW_prp_avg , D_ICW_prp_err ) + '\n'
s = s + r'      D_ICW_par    = %.3E \pm %.3E'  % (D_ICW_par_avg , D_ICW_par_err ) + '\n'
s = s + r'  Hzpe             = %.3E \pm %.3E'  % (Hzpe_avg      , Hzpe_err      ) + '\n'
s = s + r'  Hzpa             = %.3E \pm %.3E'  % (Hzpa_avg      , Hzpa_err      ) + '\n'
print (s)

f = open('time_average.txt','w') 
f.write(s) 
f.close() 


# plot energy balance
ys = [ 
       upe2dot_sum + bpe2dot_sum + upa2dot_sum + bpa2dot_sum, 
       upe2dissip_sum + bpe2dissip_sum, 
       upa2dissip_sum + bpa2dissip_sum, 
       -p_aw_sum, 
       -p_compr_sum, 
        upe2dot_sum + bpe2dot_sum + upa2dot_sum + bpa2dot_sum + upe2dissip_sum + bpe2dissip_sum + upa2dissip_sum + bpa2dissip_sum - p_aw_sum - p_compr_sum,
        np.full([tt[avg_start:avg_end].size], Wtot_dot_avg),
        np.full([tt[avg_start:avg_end].size], Daw_avg),
        np.full([tt[avg_start:avg_end].size], Dcompr_avg),
    ]
xs = [
			 tt,
			 tt,
			 tt,
			 tt,
			 tt,
			 tt,
			 tt[avg_start:avg_end],
			 tt[avg_start:avg_end],
			 tt[avg_start:avg_end],
		 ]
ls = [ 
				'', 
				'', 
				'', 
				'', 
				'', 
				'k--', 
				'', 
				'', 
				'', 
		 ]
legends = [ 
            r'$\rmd W/\rmd t$', 
            r'$D_\mr{AW}$', 
            r'$D_\mr{compr}$', 
            r'$-P_\mr{AW}$', 
            r'$-P_\mr{compr}$', 
			 r'balance', 
			 '', 
			 '', 
			 '', 
		 ]
plot_1d_many_average(xs, ys, tt[avg_start], tt[avg_end], xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'balance_all_avg.pdf')

# plot energy change
ys = [ 
       upe2_sum, 
       bpe2_sum, 
       upa2_sum, 
       bpa2_sum, 
			 np.full([tt[avg_start:avg_end].size], upe2_sum_avg),
			 np.full([tt[avg_start:avg_end].size], bpe2_sum_avg),
			 np.full([tt[avg_start:avg_end].size], upa2_sum_avg),
			 np.full([tt[avg_start:avg_end].size], bpa2_sum_avg),
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
xs = [
       tt,
       tt,
       tt,
       tt,
			 tt[avg_start:avg_end],
			 tt[avg_start:avg_end],
			 tt[avg_start:avg_end],
			 tt[avg_start:avg_end],
     ]
legends = [ 
       r'$W_{u_\+}$', 
       r'$W_{\delta B_\+}$', 
       r'$W_{u_\|}$', 
       r'$W_{\delta B_\|}$', 
			 '',
			 '',
			 '',
			 '',
     ]
plot_1d_many_average(xs, ys, tt[avg_start], tt[avg_end], xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'energy_all_avg.pdf')

# plot KAW & ICW dissipation
ys = [ 
       D_KAW, 
       D_ICW, 
			 np.full([tt[avg_start:avg_end].size], D_KAW_avg),
			 np.full([tt[avg_start:avg_end].size], D_ICW_avg),
     ]
xs = [
       tt,
       tt,
       tt[avg_start:avg_end],
       tt[avg_start:avg_end],
     ]
ls = [ 
        '', 
				'', 
				'', 
				'', 
     ]
legends = [ 
       r'$D_\mr{KAW}$', 
       r'$D_\mr{ICW}$', 
				'', 
				'', 
     ]
plot_1d_many_average(xs, ys, tt[avg_start], tt[avg_end], xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'KAW-ICW_dissipation_avg.pdf')

# plot helicity change
ys = [ 
       (zppe2_sum - zmpe2_sum)/(zppe2_sum + zmpe2_sum), 
       (zppa2_sum - zmpa2_sum)/(zppa2_sum + zmpa2_sum), 
        np.full([tt[avg_start:avg_end].size], (zppe2_sum_avg - zmpe2_sum_avg)/(zppe2_sum_avg + zmpe2_sum_avg)),
        np.full([tt[avg_start:avg_end].size], (zppa2_sum_avg - zmpa2_sum_avg)/(zppa2_sum_avg + zmpa2_sum_avg)),
     ]
ls = [ 
       '', '', 
        '', '',
     ]
xs = [
       tt, tt,
			 tt[avg_start:avg_end],
			 tt[avg_start:avg_end],
     ]
legends = [ 
       r'$H_\+ := \f{\int\rmd^3\bm{x}[(Z_\+^+)^2 - (Z_\+^-)^2]}{\int\rmd^3\bm{x}[(Z_\+^+)^2 + (Z_\+^-)^2]}$', 
       r'$H_\| := \f{\int\rmd^3\bm{x}[(Z_\|^+)^2 - (Z_\|^-)^2]}{\int\rmd^3\bm{x}[(Z_\|^+)^2 + (Z_\|^-)^2]}$', 
			 '',
			 '',
     ]
plot_1d_many_average(xs, ys, tt[avg_start], tt[avg_end], xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'helicity_avg.pdf')

# plot KAW, ICW dissipation perp and para
ys = [ 
       D_KAW_par/(D_KAW_prp+D_KAW_par), 
       D_ICW_par/(D_ICW_prp+D_ICW_par), 
        np.full([tt[avg_start:avg_end].size], D_KAW_par_avg/(D_KAW_prp_avg+D_KAW_par_avg)),
        np.full([tt[avg_start:avg_end].size], D_ICW_par_avg/(D_ICW_prp_avg+D_ICW_par_avg)),
     ]
ls = [ 
       '', 
       '', 
       '', 
       '', 
     ]
xs = [
       tt,
       tt,
       tt[avg_start:avg_end],
       tt[avg_start:avg_end],
     ]
legends = [ 
       'KAW', 
       'ICW', 
       '', 
       '', 
     ]
plot_1d_many_average(xs, ys, tt[avg_start], tt[avg_end], xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'KAW-ICW_dissipation-perp-para_avg.pdf')

##########################################################
#                    average kspectrum                   #
##########################################################
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

if nlz == nkz:
  kp_end = np.argmin(np.abs(kpbin - kpbin.max()*2./3.))
  if not is2D:
    kz_end = np.argmin(np.abs(kz[1:int(nkz/2)] - kz[1:int(nkz/2)].max()*2./3.))
else:
  kp_end = kpbin.size - 1
  kz_end = int(nkz/2)

upe2_bin_avg             = time_average(tt[avg_start:avg_end], upe2_bin            [avg_start:avg_end], axis=0)
bpe2_bin_avg             = time_average(tt[avg_start:avg_end], bpe2_bin            [avg_start:avg_end], axis=0)
upa2_bin_avg             = time_average(tt[avg_start:avg_end], upa2_bin            [avg_start:avg_end], axis=0)
bpa2_bin_avg             = time_average(tt[avg_start:avg_end], bpa2_bin            [avg_start:avg_end], axis=0)
ux2_bin_avg              = time_average(tt[avg_start:avg_end], ux2_bin             [avg_start:avg_end], axis=0)
uy2_bin_avg              = time_average(tt[avg_start:avg_end], uy2_bin             [avg_start:avg_end], axis=0)
bx2_bin_avg              = time_average(tt[avg_start:avg_end], bx2_bin             [avg_start:avg_end], axis=0)
by2_bin_avg              = time_average(tt[avg_start:avg_end], by2_bin             [avg_start:avg_end], axis=0)
zppe2_bin_avg            = time_average(tt[avg_start:avg_end], zppe2_bin           [avg_start:avg_end], axis=0)
zmpe2_bin_avg            = time_average(tt[avg_start:avg_end], zmpe2_bin           [avg_start:avg_end], axis=0)
zppa2_bin_avg            = time_average(tt[avg_start:avg_end], zppa2_bin           [avg_start:avg_end], axis=0)
zmpa2_bin_avg            = time_average(tt[avg_start:avg_end], zmpa2_bin           [avg_start:avg_end], axis=0)
p_aw_bin_avg             = time_average(tt[avg_start:avg_end], p_aw_bin            [avg_start:avg_end], axis=0)
p_compr_bin_avg          = time_average(tt[avg_start:avg_end], p_compr_bin         [avg_start:avg_end], axis=0)
dissip_aw_bin_avg        = time_average(tt[avg_start:avg_end], dissip_aw_bin       [avg_start:avg_end], axis=0)
dissip_compr_bin_avg     = time_average(tt[avg_start:avg_end], dissip_compr_bin    [avg_start:avg_end], axis=0)
dissip_KAW_bin_avg       = time_average(tt[avg_start:avg_end], dissip_KAW_bin      [avg_start:avg_end], axis=0)
dissip_ICW_bin_avg       = time_average(tt[avg_start:avg_end], dissip_ICW_bin      [avg_start:avg_end], axis=0)
ntrans_upe_upe_l_bin_avg = time_average(tt[avg_start:avg_end], ntrans_upe_upe_l_bin[avg_start:avg_end], axis=0)
ntrans_bpe_upe_l_bin_avg = time_average(tt[avg_start:avg_end], ntrans_bpe_upe_l_bin[avg_start:avg_end], axis=0)
ntrans_bpe_bpe_l_bin_avg = time_average(tt[avg_start:avg_end], ntrans_bpe_bpe_l_bin[avg_start:avg_end], axis=0)
ntrans_upe_bpe_l_bin_avg = time_average(tt[avg_start:avg_end], ntrans_upe_bpe_l_bin[avg_start:avg_end], axis=0)
ntrans_upa_upa_l_bin_avg = time_average(tt[avg_start:avg_end], ntrans_upa_upa_l_bin[avg_start:avg_end], axis=0)
ntrans_bpa_upa_l_bin_avg = time_average(tt[avg_start:avg_end], ntrans_bpa_upa_l_bin[avg_start:avg_end], axis=0)
ntrans_bpa_bpa_l_bin_avg = time_average(tt[avg_start:avg_end], ntrans_bpa_bpa_l_bin[avg_start:avg_end], axis=0)
ntrans_upa_bpa_l_bin_avg = time_average(tt[avg_start:avg_end], ntrans_upa_bpa_l_bin[avg_start:avg_end], axis=0)
ntrans_upe_upe_g_bin_avg = time_average(tt[avg_start:avg_end], ntrans_upe_upe_g_bin[avg_start:avg_end], axis=0)
ntrans_bpe_upe_g_bin_avg = time_average(tt[avg_start:avg_end], ntrans_bpe_upe_g_bin[avg_start:avg_end], axis=0)
ntrans_bpe_bpe_g_bin_avg = time_average(tt[avg_start:avg_end], ntrans_bpe_bpe_g_bin[avg_start:avg_end], axis=0)
ntrans_upe_bpe_g_bin_avg = time_average(tt[avg_start:avg_end], ntrans_upe_bpe_g_bin[avg_start:avg_end], axis=0)
ntrans_upa_upa_g_bin_avg = time_average(tt[avg_start:avg_end], ntrans_upa_upa_g_bin[avg_start:avg_end], axis=0)
ntrans_bpa_upa_g_bin_avg = time_average(tt[avg_start:avg_end], ntrans_bpa_upa_g_bin[avg_start:avg_end], axis=0)
ntrans_bpa_bpa_g_bin_avg = time_average(tt[avg_start:avg_end], ntrans_bpa_bpa_g_bin[avg_start:avg_end], axis=0)
ntrans_upa_bpa_g_bin_avg = time_average(tt[avg_start:avg_end], ntrans_upa_bpa_g_bin[avg_start:avg_end], axis=0)

ntrans_aw_l_bin_avg    = ntrans_upe_upe_l_bin_avg + ntrans_bpe_upe_l_bin_avg + ntrans_bpe_bpe_l_bin_avg + ntrans_upe_bpe_l_bin_avg
ntrans_aw_g_bin_avg    = ntrans_upe_upe_g_bin_avg + ntrans_bpe_upe_g_bin_avg + ntrans_bpe_bpe_g_bin_avg + ntrans_upe_bpe_g_bin_avg
ntrans_compr_l_bin_avg = ntrans_upa_upa_l_bin_avg + ntrans_bpa_upa_l_bin_avg + ntrans_bpa_bpa_l_bin_avg + ntrans_upa_bpa_l_bin_avg
ntrans_compr_g_bin_avg = ntrans_upa_upa_g_bin_avg + ntrans_bpa_upa_g_bin_avg + ntrans_bpa_bpa_g_bin_avg + ntrans_upa_bpa_g_bin_avg

upe2_KAW_bin_avg         = time_average(tt[avg_start:avg_end], upe2_KAW_bin        [avg_start:avg_end], axis=0)
bpe2_KAW_bin_avg         = time_average(tt[avg_start:avg_end], bpe2_KAW_bin        [avg_start:avg_end], axis=0)
upa2_KAW_bin_avg         = time_average(tt[avg_start:avg_end], upa2_KAW_bin        [avg_start:avg_end], axis=0)
bpa2_KAW_bin_avg         = time_average(tt[avg_start:avg_end], bpa2_KAW_bin        [avg_start:avg_end], axis=0)
upe2_ICW_bin_avg         = time_average(tt[avg_start:avg_end], upe2_ICW_bin        [avg_start:avg_end], axis=0)
bpe2_ICW_bin_avg         = time_average(tt[avg_start:avg_end], bpe2_ICW_bin        [avg_start:avg_end], axis=0)
upa2_ICW_bin_avg         = time_average(tt[avg_start:avg_end], upa2_ICW_bin        [avg_start:avg_end], axis=0)
bpa2_ICW_bin_avg         = time_average(tt[avg_start:avg_end], bpa2_ICW_bin        [avg_start:avg_end], axis=0)

#--------------------------------------------------------#
#                      plot 1D spectra                   #
#--------------------------------------------------------#
# kprp spectrum
ys = [ 
			 np.sum(upe2_bin_avg[:, 1:kp_end], axis=0), 
			 np.sum(bpe2_bin_avg[:, 1:kp_end], axis=0), 
       np.sum(upa2_bin_avg[:, 1:kp_end], axis=0), 
       np.sum(bpa2_bin_avg[:, 1:kp_end], axis=0),
        kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*np.sum(bpe2_bin_avg[:,1:kp_end], axis=0)[0]
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
plot_log1d_many(xs, ys, xlab='$'+kplab+'$', legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'kprp_spectra_avg.pdf')

# kprp spectrum by components
ys = [ 
       np.sum(upe2_bin_avg[:, 1:kp_end], axis=0), 
       np.sum(ux2_bin_avg [:, 1:kp_end], axis=0), 
       np.sum(uy2_bin_avg [:, 1:kp_end], axis=0), 
       np.sum(bpe2_bin_avg[:, 1:kp_end], axis=0), 
       np.sum(bx2_bin_avg [:, 1:kp_end], axis=0), 
       np.sum(by2_bin_avg [:, 1:kp_end], axis=0), 
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*np.sum(bpe2_bin_avg[:,1:kp_end], axis=0)[0]
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
plot_log1d_many(xs, ys, xlab='$'+kplab+'$', legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'kprp_spectra_components_avg.pdf')

# kprp spectrum by MRI injection rate and nonlinear transfer rate
ys = [ 
       np.sum(p_aw_bin_avg      [:,1:kp_end], axis=0),
       np.sum(p_compr_bin_avg   [:,1:kp_end], axis=0),
      -np.sum(dissip_aw_bin_avg     [:,1:kp_end], axis=0),
      -np.sum(dissip_compr_bin_avg  [:,1:kp_end], axis=0),
       np.sum(ntrans_aw_l_bin_avg   [:,1:kp_end], axis=0),
       np.sum(ntrans_aw_g_bin_avg   [:,1:kp_end], axis=0),
       np.sum(ntrans_compr_l_bin_avg[:,1:kp_end], axis=0),
       np.sum(ntrans_compr_g_bin_avg[:,1:kp_end], axis=0),
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
            r'$-\calD_\mr{AW}$', 
            r'$-\calD_\mr{compr}$', 
            r'$\calN_\mr{AW}^{<k_\+}$', 
            r'$\calN_\mr{AW}^{>k_\+}$', 
            r'$\calN_\mr{compr}^{<k_\+}$', 
            r'$\calN_\mr{compr}^{>k_\+}$', 
          ]
plot_log1d_many(xs, ys, xlab='$'+kplab+'$', legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'kprp_spectra_flux_avg.pdf')

# Elsasser fields
ys = [ 
       np.sum(zppe2_bin_avg   [:, 1:kp_end], axis=0), 
       np.sum(zmpe2_bin_avg   [:, 1:kp_end], axis=0), 
       np.sum(zppa2_bin_avg   [:, 1:kp_end], axis=0), 
       np.sum(zmpa2_bin_avg   [:, 1:kp_end], axis=0),
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*np.sum(zppe2_bin_avg[:,1:kp_end], axis=0)[0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*np.sum(zppe2_bin_avg[:,1:kp_end], axis=0)[0]
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
plot_log1d_many(xs, ys, xlab='$'+kplab+'$', legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'kprp_spectra_ELS_avg.pdf')

# KAW and ICW decomposition
ys = [ 
       np.sum(upe2_KAW_bin_avg[:, 1:kp_end], axis=0), 
       np.sum(bpe2_KAW_bin_avg[:, 1:kp_end], axis=0), 
       np.sum(upa2_KAW_bin_avg[:, 1:kp_end], axis=0), 
       np.sum(bpa2_KAW_bin_avg[:, 1:kp_end], axis=0),
       np.sum(upe2_ICW_bin_avg[:, 1:kp_end], axis=0), 
       np.sum(bpe2_ICW_bin_avg[:, 1:kp_end], axis=0), 
       np.sum(upa2_ICW_bin_avg[:, 1:kp_end], axis=0), 
       np.sum(bpa2_ICW_bin_avg[:, 1:kp_end], axis=0),
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*np.sum(bpe2_KAW_bin_avg[:,1:kp_end], axis=0)[0]
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
plot_log1d_many(xs, ys, xlab=r'$'+kplab+'$', legends=legends, ls=ls, legendloc='lower left', ylab='', vline_at = 1./rho, term=True, save=outdir+'kprp_KAW_ICW_spectra_avg.pdf')

# kz spectrum
if not is2D:
  ys = [ 
         np.sum(upe2_bin_avg[1:kz_end, :kp_end], axis=1), 
         np.sum(bpe2_bin_avg[1:kz_end, :kp_end], axis=1), 
         np.sum(upa2_bin_avg[1:kz_end, :kp_end], axis=1), 
         np.sum(bpa2_bin_avg[1:kz_end, :kp_end], axis=1), 
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
  plot_log1d_many(xs, ys, xlab='$'+kzlab+'$', legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'kz_spectra_avg.pdf')

  # Elsasser fields
  ys = [ 
         np.sum(zppe2_bin_avg[1:kz_end, :kp_end], axis=1), 
         np.sum(zmpe2_bin_avg[1:kz_end, :kp_end], axis=1), 
         np.sum(zppa2_bin_avg[1:kz_end, :kp_end], axis=1), 
         np.sum(zmpa2_bin_avg[1:kz_end, :kp_end], axis=1), 
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
  plot_log1d_many(xs, ys, xlab='$'+kzlab+'$', legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'kz_spectra_ELS_avg.pdf')

  # MRI injection rate and nonlinear transfer rate
  ys = [ 
          np.sum(p_aw_bin_avg      [1:kz_end,:kp_end], axis=1),
          np.sum(p_compr_bin_avg   [1:kz_end,:kp_end], axis=1),
         -np.sum(dissip_aw_bin_avg     [1:kz_end,:kp_end], axis=1),
         -np.sum(dissip_compr_bin_avg  [1:kz_end,:kp_end], axis=1),
          np.sum(ntrans_aw_l_bin_avg   [1:kz_end,:kp_end], axis=1),
          np.sum(ntrans_aw_g_bin_avg   [1:kz_end,:kp_end], axis=1),
          np.sum(ntrans_compr_l_bin_avg[1:kz_end,:kp_end], axis=1),
          np.sum(ntrans_compr_g_bin_avg[1:kz_end,:kp_end], axis=1),
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
              r'$-\calD_\mr{AW}$', 
              r'$-\calD_\mr{compr}$', 
              r'$\calN_\mr{AW}^{<k_\+}$', 
              r'$\calN_\mr{AW}^{>k_\+}$', 
              r'$\calN_\mr{compr}^{<k_\+}$', 
              r'$\calN_\mr{compr}^{>k_\+}$', 
            ]
  plot_log1d_many(xs, ys, xlab='$'+kzlab+'$', legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'kz_spectra_flux_avg.pdf')

#--------------------------------------------------------#
#                      plot 2D spectra                   #
#--------------------------------------------------------#
if not is2D:
  plot_log2d(upe2_bin_avg[1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$E_{u_{\+}}$', save=outdir + 'upe2_avg.pdf')
  plot_log2d(bpe2_bin_avg[1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$E_{\delta B_\+}$', save=outdir + 'bpe2_avg.pdf')
  plot_log2d(upa2_bin_avg[1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$E_{u_{\|}}$', save=outdir + 'upa2_avg.pdf')
  plot_log2d(bpa2_bin_avg[1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$E_{\delta B_\|}$', save=outdir + 'bpa2_avg.pdf')
  plot_log2d(zppe2_bin_avg[1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$E_{Z_{\+}^+}$', save=outdir + 'zppe2_avg.pdf')
  plot_log2d(zmpe2_bin_avg[1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$E_{Z_{\+}^-}$', save=outdir + 'zmpe2_avg.pdf')
  plot_log2d(zppa2_bin_avg[1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$E_{Z_{\|}^+}$', save=outdir + 'zppa2_avg.pdf')
  plot_log2d(zmpa2_bin_avg[1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$E_{Z_{\|}^-}$', save=outdir + 'zmpa2_avg.pdf')
  plot_log2d(dissip_KAW_bin_avg[1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab=r'$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$\calD_\mr{KAW}$', save=outdir + 'dissip_KAW_avg.pdf')
  plot_log2d(dissip_ICW_bin_avg[1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab=r'$'+kplab+'$', ylab='$'+kzlab+'$', 
      title=r'$\calD_\mr{ICW}$', save=outdir + 'dissip_ICW_avg.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ekprp_avg.txt'  , np.column_stack((kpbin[:kp_end], 
                                                       np.sum(upe2_bin_avg            [:kz_end,:kp_end], axis=0),
                                                       np.sum(bpe2_bin_avg            [:kz_end,:kp_end], axis=0),
                                                       np.sum(upa2_bin_avg            [:kz_end,:kp_end], axis=0),
                                                       np.sum(bpa2_bin_avg            [:kz_end,:kp_end], axis=0),
                                                       np.sum(ux2_bin_avg             [:kz_end,:kp_end], axis=0),
                                                       np.sum(uy2_bin_avg             [:kz_end,:kp_end], axis=0),
                                                       np.sum(bx2_bin_avg             [:kz_end,:kp_end], axis=0),
                                                       np.sum(by2_bin_avg             [:kz_end,:kp_end], axis=0),
                                                       np.sum(zppe2_bin_avg           [:kz_end,:kp_end], axis=0),
                                                       np.sum(zmpe2_bin_avg           [:kz_end,:kp_end], axis=0),
                                                       np.sum(zppa2_bin_avg           [:kz_end,:kp_end], axis=0),
                                                       np.sum(zmpa2_bin_avg           [:kz_end,:kp_end], axis=0),
                                                       np.sum(p_aw_bin_avg            [:kz_end,:kp_end], axis=0),
                                                       np.sum(p_compr_bin_avg         [:kz_end,:kp_end], axis=0),
                                                       np.sum(dissip_aw_bin_avg       [:kz_end,:kp_end], axis=0),
                                                       np.sum(dissip_compr_bin_avg    [:kz_end,:kp_end], axis=0),
                                                       np.sum(dissip_KAW_bin_avg      [:kz_end,:kp_end], axis=0),
                                                       np.sum(dissip_ICW_bin_avg      [:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_upe_upe_l_bin_avg[:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_bpe_upe_l_bin_avg[:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_bpe_bpe_l_bin_avg[:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_upe_bpe_l_bin_avg[:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_upa_upa_l_bin_avg[:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_bpa_upa_l_bin_avg[:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_bpa_bpa_l_bin_avg[:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_upa_bpa_l_bin_avg[:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_upe_upe_g_bin_avg[:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_bpe_upe_g_bin_avg[:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_bpe_bpe_g_bin_avg[:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_upe_bpe_g_bin_avg[:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_upa_upa_g_bin_avg[:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_bpa_upa_g_bin_avg[:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_bpa_bpa_g_bin_avg[:kz_end,:kp_end], axis=0),
                                                       np.sum(ntrans_upa_bpa_g_bin_avg[:kz_end,:kp_end], axis=0),
                                                       np.sum(upe2_KAW_bin_avg        [:kz_end,:kp_end], axis=0), 
                                                       np.sum(bpe2_KAW_bin_avg        [:kz_end,:kp_end], axis=0), 
                                                       np.sum(upa2_KAW_bin_avg        [:kz_end,:kp_end], axis=0), 
                                                       np.sum(bpa2_KAW_bin_avg        [:kz_end,:kp_end], axis=0),
                                                       np.sum(upe2_ICW_bin_avg        [:kz_end,:kp_end], axis=0), 
                                                       np.sum(bpe2_ICW_bin_avg        [:kz_end,:kp_end], axis=0), 
                                                       np.sum(upa2_ICW_bin_avg        [:kz_end,:kp_end], axis=0), 
                                                       np.sum(bpa2_ICW_bin_avg        [:kz_end,:kp_end], axis=0),
                                                     )), fmt='%E')
if not is2D:
  np.savetxt(outdir + 'Ekz_avg.txt'  , np.column_stack((kz[:kz_end], 
                                                         np.sum(upe2_bin_avg            [:kz_end,:kp_end], axis=1),
                                                         np.sum(bpe2_bin_avg            [:kz_end,:kp_end], axis=1),
                                                         np.sum(upa2_bin_avg            [:kz_end,:kp_end], axis=1),
                                                         np.sum(bpa2_bin_avg            [:kz_end,:kp_end], axis=1),
                                                         np.sum(ux2_bin_avg             [:kz_end,:kp_end], axis=1),
                                                         np.sum(uy2_bin_avg             [:kz_end,:kp_end], axis=1),
                                                         np.sum(bx2_bin_avg             [:kz_end,:kp_end], axis=1),
                                                         np.sum(by2_bin_avg             [:kz_end,:kp_end], axis=1),
                                                         np.sum(zppe2_bin_avg           [:kz_end,:kp_end], axis=1),
                                                         np.sum(zmpe2_bin_avg           [:kz_end,:kp_end], axis=1),
                                                         np.sum(zppa2_bin_avg           [:kz_end,:kp_end], axis=1),
                                                         np.sum(zmpa2_bin_avg           [:kz_end,:kp_end], axis=1),
                                                         np.sum(p_aw_bin_avg            [:kz_end,:kp_end], axis=1),
                                                         np.sum(p_compr_bin_avg         [:kz_end,:kp_end], axis=1),
                                                         np.sum(dissip_aw_bin_avg       [:kz_end,:kp_end], axis=1),
                                                         np.sum(dissip_compr_bin_avg    [:kz_end,:kp_end], axis=1),
                                                         np.sum(dissip_KAW_bin_avg      [:kz_end,:kp_end], axis=1),
                                                         np.sum(dissip_ICW_bin_avg      [:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_upe_upe_l_bin_avg[:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_bpe_upe_l_bin_avg[:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_bpe_bpe_l_bin_avg[:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_upe_bpe_l_bin_avg[:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_upa_upa_l_bin_avg[:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_bpa_upa_l_bin_avg[:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_bpa_bpa_l_bin_avg[:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_upa_bpa_l_bin_avg[:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_upe_upe_g_bin_avg[:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_bpe_upe_g_bin_avg[:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_bpe_bpe_g_bin_avg[:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_upe_bpe_g_bin_avg[:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_upa_upa_g_bin_avg[:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_bpa_upa_g_bin_avg[:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_bpa_bpa_g_bin_avg[:kz_end,:kp_end], axis=1),
                                                         np.sum(ntrans_upa_bpa_g_bin_avg[:kz_end,:kp_end], axis=1),
                                                         np.sum(upe2_KAW_bin_avg        [:kz_end,:kp_end], axis=1), 
                                                         np.sum(bpe2_KAW_bin_avg        [:kz_end,:kp_end], axis=1), 
                                                         np.sum(upa2_KAW_bin_avg        [:kz_end,:kp_end], axis=1), 
                                                         np.sum(bpa2_KAW_bin_avg        [:kz_end,:kp_end], axis=1),
                                                         np.sum(upe2_ICW_bin_avg        [:kz_end,:kp_end], axis=1), 
                                                         np.sum(bpe2_ICW_bin_avg        [:kz_end,:kp_end], axis=1), 
                                                         np.sum(upa2_ICW_bin_avg        [:kz_end,:kp_end], axis=1), 
                                                         np.sum(bpa2_ICW_bin_avg        [:kz_end,:kp_end], axis=1),
                                                       )), fmt='%E')


del upe2_bin
del bpe2_bin
