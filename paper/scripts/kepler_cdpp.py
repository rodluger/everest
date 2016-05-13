#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
kepler_cdpp.py
--------------

Computes the raw 6-hr CDPP for all original `Kepler` targets.

..warning:: This was copied over from an older version, and may need to be tweaked slightly.

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

# Start up kplr
client = kplr.API()

# Get all stars
stars = list(np.loadtxt(os.path.join(EVEREST_ROOT, 'tables', 'keplerstars.csv'), dtype = int)) 
nstars = len(stars)

# Remove ones we've done
with warnings.catch_warnings():
  warnings.simplefilter("ignore")
  done = np.loadtxt(os.path.join('CDPP', 'kepler.tsv'), dtype = float)
if len(done):
  done = [int(s) for s in done[:,0]]
stars = list(set(stars) - set(done))
n = len(done) + 1

# Open the output file
with open(os.path.join('CDPP', 'kepler.tsv'), 'a') as outfile:

  # Loop over all to get the CDPP
  for star in stars:

    # Progress
    sys.stdout.write('\rRunning target %d/%d...' % (n, nstars))
    sys.stdout.flush()
    n += 1
    
    # Get the cdpp
    try:
      s = kplr.K2SFF(star)
    except (HTTPError, TypeError, ValueError):
      continue
    
    # Get a random quarter
    lc = random.choice(s.get_light_curves(short_cadence = False))
    
    # Extract the timeseries
    with lc.open() as infile:
      time = infile[1].data.field('TIME')
      flux = infile[1].data.field('SAP_FLUX')
      bad = np.where(np.isnan(time) | np.isnan(flux))
      time = np.delete(time, bad)
      flux = np.delete(flux, bad)
    
    rms = RMS(flux / np.median(flux), remove_outliers = True)
    flux_sv2 = flux - savgol_filter(flux, 49, 2) + np.median(flux)
    rms_sv2 = RMS(flux_sv2 / np.nanmedian(flux_sv2), remove_outliers = True)  
    print("{:>09d} {:>15.3f} {:>15.3f}".format(star, rms, rms_sv2), file = outfile)
    
    # Delete the lightcurve on disk
    shutil.rmtree(os.path.join(kplr.config.KPLR_ROOT, 'data', 'lightcurves', '%09d' % star))