# -*- coding: utf-8 -*-

#####################################################
## main program for making plots from AstroGK data ##
#####################################################
import numpy as np
from scipy.integrate import simps, trapz
from scipy import interpolate

from numba import jit

def time_average(x, y, axis=0): # x: 1D array, y: any-D array 
  return trapz(y, x, axis=axis)/(x[-1] - x[0])
  # return trapz(y, x, axis=axis)/(x[-1] - x[0])

avg_start = 101
avg_end   = -1

##########################################################
#              average energy time evolution             #
##########################################################
print('\nplotting energy\n')
outdir = './fig_energy/'

# load data
time      = np.transpose(np.loadtxt(outdir + 'energies.txt'    ))[0]
u2dot_sum = np.transpose(np.loadtxt(outdir + 'energies.txt'    ))[3]
b2dot_sum = np.transpose(np.loadtxt(outdir + 'energies.txt'    ))[4]

energy_dot = u2dot_sum + b2dot_sum

energy_dot_avg = time_average(time[avg_start:avg_end], energy_dot[avg_start:avg_end], axis=0)
print ('energy_dot     = %.3E' % energy_dot_avg + ' ! This must be zero when stationary')
