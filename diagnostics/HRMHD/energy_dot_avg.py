# -*- coding: utf-8 -*-

#####################################################
## main program for making plots from AstroGK data ##
#####################################################
import numpy as np
from scipy.integrate import simps, trapz
from scipy import interpolate

from numba import jit

@jit
def time_average(x, y, axis=0): # x: 1D array, y: any-D array 
  return trapz(y, x, axis=axis)/(x[-1] - x[0])
  # return trapz(y, x, axis=axis)/(x[-1] - x[0])

avg_start = 77
avg_end   = -1

##########################################################
#              average energy time evolution             #
##########################################################
print('\nplotting energy\n')
outdir = './fig_energy/'

# load data
time        = np.transpose(np.loadtxt(outdir + 'energies.txt'    ))[0]
upe2dot_sum = np.transpose(np.loadtxt(outdir + 'energies.txt'    ))[5]
bpe2dot_sum = np.transpose(np.loadtxt(outdir + 'energies.txt'    ))[6]

energy_dot = upe2dot_sum + bpe2dot_sum

energy_dot_avg = time_average(time[avg_start:avg_end], energy_dot[avg_start:avg_end], axis=0)
print ('energy_dot     = %.3E' % energy_dot_avg + ' ! This must be zero when stationary')
