#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`pipelines.py` - Pipeline comparison
--------------------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from ...math import CDPP6
from ...config import EVEREST_SRC
import os, sys, shutil
import k2plr
from k2plr.config import KPLR_ROOT
from urllib.error import HTTPError
import matplotlib.pyplot as pl
import numpy as np
import warnings
import logging
log = logging.getLogger(__name__)

Pipelines = ['everest1', 'k2sff', 'k2sc']

def plot(ID, pipeline = 'everest1', show = True):
  '''
  
  '''
  
  # Get the data
  log.info('Downloading %s light curve for %d...' % (pipeline, ID))
  if pipeline.lower() == 'everest1':
    s = k2plr.EVEREST(ID, version = 1)
    time = s.time
    flux = s.flux
  elif pipeline.lower() == 'k2sff':
    s = k2plr.K2SFF(ID)
    time = s.time
    flux = s.fcor
  elif pipeline.lower() == 'k2sc':
    s = k2plr.K2SC(ID)
    time = s.time
    flux = s.pdcflux
  else:
    raise ValueError('Invalid pipeline: `%s`.' % pipeline)
  
  # Remove nans
  mask = np.where(np.isnan(flux))[0]
  time = np.delete(time, mask)
  flux = np.delete(flux, mask)
  
  # Plot it
  fig, ax = pl.subplots(1, figsize = (10, 4))
  fig.subplots_adjust(bottom = 0.15)
  ax.plot(time, flux, "k.", markersize = 3, alpha = 0.5)
  
  # Axis limits
  N = int(0.995 * len(flux))
  hi, lo = flux[np.argsort(flux)][[N,-N]]
  fsort = flux[np.argsort(flux)]
  pad = (hi - lo) * 0.1
  ylim = (lo - pad, hi + pad)
  ax.set_ylim(ylim)
  
  # Show the CDPP
  ax.annotate('%.2f ppm' % CDPP6(flux), 
              xy = (0.98, 0.975), xycoords = 'axes fraction', 
              ha = 'right', va = 'top', fontsize = 12, color = 'r', zorder = 99)
  
  # Appearance
  ax.margins(0, None)
  ax.set_xlabel("Time (BJD - 2454833)", fontsize = 16)
  ax.set_ylabel("%s Flux" % pipeline.upper(), fontsize = 16)
  fig.canvas.set_window_title("%s: EPIC %d" % (pipeline.upper(), ID))
  
  if show:
    pl.show()
    pl.close()
  else:
    return fig, ax

def get_cdpp(campaign, pipeline = 'everest1'):
  '''
  
  '''
  
  # Check pipeline
  assert pipeline.lower() in Pipelines, 'Invalid pipeline: `%s`.' % pipeline
  
  # Create file if it doesn't exist
  file = os.path.join(EVEREST_SRC, 'missions', 'k2', 'tables', 'c%02d_%s.tsv' % (int(campaign), pipeline))
  if not os.path.exists(file):
    open(file, 'a').close()

  # Get all EPIC stars
  from .aux import GetK2Campaign
  stars = GetK2Campaign(campaign, epics_only = True) 
  nstars = len(stars)

  # Remove ones we've done
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    done = np.loadtxt(file, dtype = float)
  if len(done):
    done = [int(s) for s in done[:,0]]
  stars = list(set(stars) - set(done))
  n = len(done) + 1

  # Open the output file
  with open(file, 'a', 1) as outfile:

    # Loop over all to get the CDPP
    for EPIC in stars:

      # Progress
      sys.stdout.write('\rRunning target %d/%d...' % (n, nstars))
      sys.stdout.flush()
      n += 1
      
      # Get the CDPP
      try:
        if pipeline.lower() == 'everest1':
          flux = k2plr.EVEREST(EPIC, version = 1).flux
        elif pipeline.lower() == 'k2sff':
          flux = k2plr.K2SFF(EPIC).fcor
        elif pipeline.lower() == 'k2sc':
          flux = k2plr.K2SC(EPIC).pdcflux
        mask = np.where(np.isnan(flux))[0]
        flux = np.delete(flux, mask)
        cdpp = CDPP6(flux)
      except (HTTPError, TypeError, ValueError):
        print("{:>09d} {:>15.3f}".format(EPIC, 0), file = outfile)
        continue

      # Log to file
      print("{:>09d} {:>15.3f}".format(EPIC, cdpp), file = outfile)