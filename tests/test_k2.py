#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test_k2.py
----------

Test the K2 de-trending core and the user interface.

'''

import everest
from everest.config import EVEREST_DAT
from k2plr.config import KPLR_ROOT
import os
import numpy as np
import shutil

def test_c1():
  '''
  Testing the de-trending for K2 campaign 1, a regular campaign with 2 breakpoints
  
  '''
  
  # Copy the TPF file to the correct directory
  dest = os.path.join(KPLR_ROOT, 'data', 'k2', 'target_pixel_files', '201367065')
  if not os.path.exists(dest):
    os.makedirs(dest)
  orig = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ktwo201367065-c01_lpd-targ.fits.gz')
  shutil.copy(orig, dest)

  # Copy the K2SFF file to the correct directory
  dest = os.path.join(KPLR_ROOT, "data", "k2sff", '201367065')
  if not os.path.exists(dest):
    os.makedirs(dest)
  orig = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hlsp_k2sff_k2_lightcurve_201367065-c01_kepler_v1_llc.fits')
  shutil.copy(orig, dest)
  
  # Create a CBV file and place it in the correct directory
  #
  # Note that 'Xc01.txt' was created as follows:
  #
  # data = np.load('~/.everest2/k2/cbv/c01/X.npz')
  # np.savetxt('Xc01.txt', np.hstack([data['time'].reshape(-1,1), 
  #            data['X']]), header = ",".join([str(b) for b in data['breakpoints']]))
  #
  data = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Xc01.txt'))
  with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Xc01.txt'), 'r') as f:
    header = f.readline()
    breakpoints = [int(b) for b in header[2:-1].split(",")]
  time = data[:,0]
  X = data[:,1:]
  dest = os.path.join(EVEREST_DAT, 'k2', 'cbv', 'c01')
  if not os.path.exists(dest):
    os.makedirs(dest)
  np.savez(os.path.join(dest, 'X.npz'), time = time, X = X, breakpoints = breakpoints)
  
  # Run the de-trending
  star1 = everest.rPLD(201367065, clobber = True, mission = 'k2', debug = False,
                       giter = 1, gmaxf = 3, lambda_arr = [1e0, 1e5, 1e10], oiter = 3,
                       pld_order = 2, get_hires = False, get_nearby = False)
  
  # Check!
  print("De-trended CDPP: %.3f ppm" % star1.cdpp)
  assert (star1.cdpp > 15.) and (star1.cdpp < 19.), "De-trended CDPP is different from benchmark value (17.302 ppm)."

  # Publish
  star1.publish()

  # Load the FITS file
  star2 = everest.Everest(201367065)
  
  # Compute the model
  star2.compute()
  
  # Check!
  print("De-trended CDPP (Stage 2): %.3f ppm" % star2.cdpp)
  assert np.abs(star2.cdpp - star1.cdpp) / star1.cdpp < 0.001, "De-trended CDPP is different from Stage 1 value (%.3f ppm)." % star1.cdpp
  
def test_c9():
  '''
  Testing the de-trending for K2 campaign 9, a split campaign (91, 92) with no other breakpoints
  
  '''
  
  # Copy the TPF files to the correct directory
  dest = os.path.join(KPLR_ROOT, 'data', 'k2', 'target_pixel_files', '221312395')
  if not os.path.exists(dest):
    os.makedirs(dest)
  orig = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ktwo221312395-c91_lpd-targ.fits.gz')
  shutil.copy(orig, dest)
  orig = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ktwo221312395-c92_lpd-targ.fits.gz')
  shutil.copy(orig, dest)

  # Copy the K2SFF files to the correct directory
  dest = os.path.join(KPLR_ROOT, "data", "k2sff", '221312395')
  if not os.path.exists(dest):
    os.makedirs(dest)
  orig = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hlsp_k2sff_k2_lightcurve_221312395-c91_kepler_v1_llc.fits')
  shutil.copy(orig, dest)
  orig = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hlsp_k2sff_k2_lightcurve_221312395-c92_kepler_v1_llc.fits')
  shutil.copy(orig, dest)

  # Run the de-trending
  star1 = everest.rPLD(221312395, clobber = True, mission = 'k2', debug = False,
                       giter = 1, gmaxf = 3, lambda_arr = [1e0, 1e5, 1e10], oiter = 3,
                       pld_order = 2, get_hires = False, get_nearby = False, aperture = 'k2sff_13')
  
  # Check!
  print("De-trended CDPP: %.3f ppm" % star1.cdpp)
  assert (star1.cdpp > 200.) and (star1.cdpp < 400.), "De-trended CDPP is different from benchmark value (352.3 ppm)."

  # Publish
  star1.publish()

  # Load the FITS file
  star2 = everest.Everest(221312395)
  
  # Compute the model
  star2.compute()
  
  # Check!
  print("De-trended CDPP (Stage 2): %.3f ppm" % star2.cdpp)
  assert np.abs(star2.cdpp - star1.cdpp) / star1.cdpp < 0.001, "De-trended CDPP is different from Stage 1 value (%.3f ppm)." % star1.cdpp
