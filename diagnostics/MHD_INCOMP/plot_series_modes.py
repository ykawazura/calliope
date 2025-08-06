# -*- coding: utf-8 -*-
from load import *
import sys
sys.path.append('../')
from plots import *

print('\nplotting series modes\n')
outdir = './fig_series_modes/'

tab10 = plt.get_cmap('tab10').colors

if shear_flg == 1:
  b00 = np.max(b0[0])
  np.savetxt(outdir + 'modes.txt' , modes_unique*b00, fmt='%E')
else:
  np.savetxt(outdir + 'modes.txt' , modes_unique, fmt='%E')

# omg = np.fft.fftfreq(tt_unique.size, d=tt_unique[1] - tt_unique[0])
# window = 0.54 - 0.46*np.cos(2.0*np.pi*((tt_unique - tt_unique[0])/(tt_unique[-1] - tt_unique[0])))

# plot energy balance
for i, fn in enumerate(fldnames_unique):
  print (fn)
  ys = []; xs = []; ls = []; legends = []
  for j, mode in enumerate(modes_unique):
    print (mode)
    ys.append(np.real(f_[fn].T[j]))
    ys.append(np.imag(f_[fn].T[j]))
    xs.append(tt_unique)
    ls.append('')
    legends.append('Re '+ label_series[fn] + r' $(k_x, k_y, k_z) = (%.1f, %.1f, %.1f)$' % (mode[0], mode[1], mode[2]))

  plot_1d_many(xs, ys[:, 0], xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + fn + '.pdf')

  # Output will be (time, $(real of mode 1), $(imag of mode 1_, $(real of mode 2), $(imag of mode 2), ...)
  np.savetxt(outdir + fn + '.txt' , np.vstack((tt_unique, np.asarray(ys))).T, fmt='%E')
