#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
run.py
------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .plot import Plot
from .data import GetK2Stars, GetK2Data
from .compute import Compute
from .utils import ExceptionHook, ExceptionHookPDB
import sys
import logging
log = logging.getLogger(__name__)

def Run(EPIC, debug = True, **kwargs):
  '''
  Compute and plot data for a given target.
  
  '''
  
  # Set up our custom exception handlers
  if debug:
    sys.excepthook = ExceptionHookPDB
  else:
    sys.excepthook = ExceptionHook
  
  # Compute
  data = Compute(EPIC, **kwargs)
  
  # Plot the output
  Plot(data)

def RunCampaign(campaign, **kwargs):
  '''
  Compute and plot data for all targets in a given campaign.
  
  '''
  
  # Get all star IDs for this campaign
  stars = GetK2Stars()[campaign]
  nstars = len(stars)
  
  # Download the TPF data for each one
  for i, EPIC in enumerate(stars):
    print("Downloading data for EPIC %d (%d/%d)..." % (EPIC, i + 1, nstars))
    GetK2Data(EPIC)
  
  # TODO: Now submit a cluster job