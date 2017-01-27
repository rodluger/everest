#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test_detrend.py
---------------

Test the de-trending core.

'''

import everest
from k2plr.config import KPLR_ROOT
import os
import shutil

def test_detrend():
  '''
  
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

  # Run the de-trending
  star = everest.rPLD(201367065, clobber = True, debug = True, mission = 'k2',
                      giter = 1, gmaxf = 3, lambda_arr = [1e0, 1e5, 1e10], oiter = 3,
                      pld_order = 2, get_hires = False, get_nearby = False)
  
  # Check!
  print("De-trended CDPP: %.3f ppm" % star.cdpp)
  assert (star.cdpp > 15.) and (star.cdpp < 19.), "De-trended CDPP is different from benchmark value (17.302 ppm)."