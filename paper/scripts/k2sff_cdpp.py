#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
k2sff_cdpp.py
-------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
EVEREST_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, EVEREST_ROOT)
import everest
from everest.utils import RMS
import kplr
from kplr.config import KPLR_ROOT
import random
import numpy as np
import shutil
import subprocess
import warnings
from urllib.error import HTTPError
from scipy.signal import savgol_filter

for campaign in range(0,7):
  
  print("\nRunning campaign %02d..." % campaign)
  
  # Create file if it doesn't exist
  if not os.path.exists(os.path.join('CDPP', 'k2sff_C%02d.tsv' % campaign)):
    open(os.path.join('CDPP', 'k2sff_C%02d.tsv' % campaign), 'a').close()
  
  # Get all EPIC stars
  stars = list(np.loadtxt(os.path.join(EVEREST_ROOT, 'tables', 'C%02d.csv' % campaign), dtype = int))  

  # Now remove candidates and EBs
  ebs = set([int(eb.epic) for eb in everest.GetK2EBs()])
  planets = set([int(planet.epic_name[5:]) for planet in everest.GetK2Planets()])
  stars = list(set(stars) - (planets | ebs))
  nstars = len(stars)

  # Remove ones we've done
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    done = np.loadtxt(os.path.join('CDPP', 'k2sff_C%02d.tsv' % campaign), dtype = float)
  if len(done):
    done = [int(s) for s in done[:,0]]
  stars = list(set(stars) - set(done))
  n = len(done) + 1

  # Open the output file
  with open(os.path.join('CDPP', 'k2sff_C%02d.tsv' % campaign), 'a') as outfile:

    # Loop over all to get the CDPP
    for star in stars:

      # Progress
      sys.stdout.write('\rRunning target %d/%d...' % (n, nstars))
      sys.stdout.flush()
      n += 1
      
      # Get the cdpp
      try:
        s = kplr.K2SFF(star)
      except HTTPError:
        print("{:>09d} {:>15.3f} {:>15.3f}".format(star, 0, 0), file = outfile)
        continue
      rms = RMS(s.fcor / np.nanmedian(s.fcor), remove_outliers = True)
      flux_sv2 = s.fcor - savgol_filter(s.fcor, 49, 2) + np.nanmedian(s.fcor)
      rms_sv2 = RMS(flux_sv2 / np.nanmedian(flux_sv2), remove_outliers = True)
      print("{:>09d} {:>15.3f} {:>15.3f}".format(star, rms, rms_sv2), file = outfile)
      # Delete the lightcurve on disk
      os.remove(s._file)
      try:
        os.rmdir(os.path.dirname(s._file))
      except:
        # I'm getting an "OSError: [Errno 39] Directory not empty:" error
        # when trying to do this sometimes... But the directory *is* empty.
        # No idea why!
        pass