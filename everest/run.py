#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
run.py
------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .plot import Plot
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

def RunCampaign(campaign, debug = False, **kwargs):
  '''
  Compute and plot data for all targets in a given campaign.
  
  '''
  
  # Set up our custom exception handlers
  if debug:
    sys.excepthook = ExceptionHookPDB
  else:
    sys.excepthook = ExceptionHook
  
  raise Exception("Not yet implemented!")
