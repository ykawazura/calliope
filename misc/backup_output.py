# -*- coding: utf-8 -*-
'''
This script copies output files. You need to give two arguments, MODEL and restart_num. 
Otherwise, you will be asked to input them via stdin.
For example, running 
  [$] python misc/backup_output.py MHD_INCOMP 1 
at the root directory will copy 
  calliope.out.std -> calliope.out.std1
  calliope.in      -> calliope.in1
  calliope.out.nc  -> calliope.out.nc1
  out2d            -> out2d1
  restart          -> restart1
and so on.
When restart_num = 'all', NetCDF files are merged via nccat, and binary files in directories,
e.g., out2d, are marged.

'''
import warnings
warnings.filterwarnings('ignore')

import sys
import os.path
import os
import shutil
import glob
import numpy as np

argv = sys.argv

if len(argv) != 3:
  print ('MODEL = ')
  MODEL = input()

  print ('Restart number = ')
  restart_num = input()
else:
  MODEL, restart_num = argv[1:3]


if MODEL == 'MHD_INCOMP':
  filenames = [
                'calliope.out.std',
                'calliope.in',
                'cfl.dat',
                'rms.dat',
                'calliope.out.nc',
                'calliope.modes.out',
                'calliope.out.SF2.nc',
                'calliope.out.kpar.nc',
                'calliope.out.nltrans.nc',
              ]
  dirnames =  ['out2d', 'out3d', 'restart']


if MODEL == 'RRMHD' or MODEL == 'HRMHD' or MODEL == 'WILL' or MODEL == 'HIGH_BETA_RMHD':
  filenames = [
                'calliope.out.std',
                'calliope.in',
                'cfl.dat',
                'calliope.out.nc',
                #'calliope.modes.out',
              ]
  dirnames =  ['out2d', 'out3d', 'restart']



if restart_num != 'all':
  #==============================
  # Backup output for restart_num
  #==============================
  if any([os.path.exists(filename+restart_num) for filename in filenames]):
    print ('There is already a backup with restart_num = ' + restart_num)
  else:
    for filename in filenames:
      print ('copying ' + filename)
      shutil.copy(filename, filename+restart_num)

    for dirname in dirnames:
      print ('copying ' + dirname)
      os.system("cp -rf " + dirname + " " + dirname+restart_num)
else:
  #==============================
  # Merge all the backup output
  #==============================
  import subprocess
  import re

  def atoi(text):
    return int(text) if text.isdigit() else text

  def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]


  # Clear all merged files 
  command = "rm *tmp; rm -rf *_all"
  result = subprocess.run([command], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  if result.returncode == 0:
    print (command)
  else:
    print(result.stdout)

  # Merge netCDF files
  for filename in filenames:
    if '.nc' in filename:
      command = "ncrcat "

      for filename_expand in sorted(glob.glob(filename+'?*'), key=natural_keys):
        command += filename_expand + ' '

      command += filename+"_all"

      result = subprocess.run([command], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      if result.returncode == 0:
        print (command)
      else:
        print(result.stdout)

    if '.dat' in filename:
      command = "cat " + filename + "* > " + filename + "_all"
      result = subprocess.run([command], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

      if result.returncode == 0:
        print (command)
      else:
        print(result.stdout)


  # Merge time series file
  command = 'cat '

  flg = False
  for filename_expand in sorted(glob.glob('calliope.modes.out?*'), key=natural_keys):
    flg = True
    command += filename_expand + ' '

  command += '> calliope.modes.out_all'

  if flg:
    result = subprocess.run([command], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode == 0:
      print (command)
    else:
      print(result.stdout)


  # Merge binary files in subdirectories
  for dirname in dirnames:
    if dirname != 'restart':

      result = subprocess.run(["mkdir " + dirname + '_all'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      if result.returncode != 0:
        print(result.stdout)

      
      for filename in set([f.split('/')[1] for f in glob.glob(dirname+'?*/*')]):
        command = 'cat '

        for fn in [x for x in sorted(glob.glob(dirname+'?*/*'), key=natural_keys) if filename in x]:
          command += fn + ' '

        command += '> ' + dirname + '_all/' + filename

        result = subprocess.run([command], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode == 0:
          print (command)
        else:
          print(result.stdout)
