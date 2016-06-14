#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
everest_cdpp.py
---------------

Computes the 6-hr CDPP for all the `everest` de-trended light curves.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
import everest
from everest.config import EVEREST_SRC, EVEREST_DAT
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

run_name = 'default'

for campaign in range(8):
  
  print("\nRunning campaign %02d..." % campaign)
  
  # Create file if it doesn't exist
  if not os.path.exists(os.path.join('CDPP', 'everest_C%02d.tsv' % campaign)):
    open(os.path.join('CDPP', 'everest_C%02d.tsv' % campaign), 'a').close()
  
  # Get all EPIC stars
  stars = list(np.loadtxt(os.path.join(EVEREST_SRC, 'tables', 'C%02d.csv' % campaign), dtype = int))  
  nstars = len(stars)

  # Remove ones we've done
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    try:
      done, _, _ = np.loadtxt(os.path.join('CDPP', 'everest_C%02d.tsv' % campaign), dtype = float, unpack = True)
    except ValueError:
      done = []
  if len(done):
    done = [int(d) for d in done]
  stars = list(set(stars) - set(done))
  n = len(done) + 1

  # Open the output file (1)
  with open(os.path.join('CDPP', 'everest_C%02d.tsv' % campaign), 'a') as feverest:

    # Open the output file (2)
    with open(os.path.join('CDPP', 'k2raw_C%02d.tsv' % campaign), 'a') as fraw:

      # Loop over all to get the CDPP
      for star in stars:

        # Progress
        sys.stdout.write('\rRunning target %d/%d...' % (n, nstars))
        sys.stdout.flush()
        n += 1
      
        # Get the cdpp
        try:
          data = np.load(os.path.join(EVEREST_DAT, 'output', 'C%02d' % campaign, '%d' % star, run_name, 'data.npz'))
        except KeyboardInterrupt:
          sys.exit()
        #except:
        #  continue
        rms_raw, rms_raw_savgol, rms_evr, rms_evr_savgol, rms_pht = data['rms']
        satsev = data['satsev'][()]
        crwdsev = data['crwdsev'][()]
        kepmag = data['kepmag'][()]
      
        print("{:>09d} {:>15.3f} {:>15.3f}".format(star, rms_evr, rms_evr_savgol), file = feverest)
        print("{:>09d} {:>15.3f} {:>15.3f} {:>15.3f} {:>15.3f} {:>01d} {:>01d}".format(star, kepmag, rms_raw, rms_raw_savgol, rms_pht, satsev, crwdsev), file = fraw)