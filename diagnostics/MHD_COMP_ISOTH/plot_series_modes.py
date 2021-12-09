# -*- coding: utf-8 -*-
from load import *
from plots import *

print('\nplotting series modes\n')
outdir = './fig_series_modes/'

# plot energy balance

for i, fn in enumerate(fldnames_unique):
  ys = []; xs = []; ls = []; legends = []
  for j, mode in enumerate(modes_unique):
    ys.append(np.real(f_[fn].T[j]))
    xs.append(tt_unique)
    ls.append('')
    legends.append('Re '+ label_series[fn] + r' $(k_x, k_y, k_z) = (%.2f, %.2f, %.2f)$' % (mode[0], mode[1], mode[2]))
  plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + fn + '.pdf')
