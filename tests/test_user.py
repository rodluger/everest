#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test_user.py
------------

Test the user interface.

'''

import everest
import os
import shutil

def test_user():
  '''
  
  '''
  
  # Copy the sample K2 FITS file to the correct directory
  path = everest.missions.k2.TargetDirectory(201367065, 1)
  if not os.path.exists(path):
    os.makedirs(path)
  dest = os.path.join(path, everest.missions.k2.FITSFile(201367065, 1))
  orig = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hlsp_everest_k2_llc_201367065-c01_kepler_v2.0_lc.fits')
  shutil.copy(orig, dest)

  # Load the FITS file
  star = everest.Everest(201367065)
  
  # Compute the model
  star.compute()