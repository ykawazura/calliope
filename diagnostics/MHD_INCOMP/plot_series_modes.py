# -*- coding: utf-8 -*-
from load import *
from plots import *

print('\nplotting series modes\n')
outdir = './fig_series_modes/'

if shear_flg == 1:
  b00 = np.max(b0[0])
  np.savetxt(outdir + 'modes.txt' , modes_unique*b00, fmt='%E')
else:
  np.savetxt(outdir + 'modes.txt' , modes_unique, fmt='%E')

# plot energy balance
for i, fn in enumerate(fldnames_unique):
  ys = []; xs = []; ls = []; legends = []
  for j, mode in enumerate(modes_unique):
    ys.append(np.abs(f_[fn].T[j]))
    xs.append(tt_unique)
    ls.append('')
    legends.append('Re '+ label_series[fn] + r' $(k_x, k_y, k_z) = (%.1f, %.1f, %.1f)$' % (mode[0], mode[1], mode[2]))
  plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + fn + '.pdf')
  np.savetxt(outdir + fn + '.txt' , np.vstack((tt_unique, np.asarray(ys))).T, fmt='%E')
