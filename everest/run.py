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
from .utils import ExceptionHook, ExceptionHookPDB, FunctionWrapper
from .pool import Pool
import os
EVEREST_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
from kplr.config import KPLR_ROOT
import sys
import subprocess
import imp

def Run(EPIC, kwargs_file = os.path.join(EVEREST_ROOT, 'scripts', 'kwargs.py')):
  '''
  Compute and plot data for a given target.
  
  '''
  
  # Get the kwargs
  kwargs = imp.load_source("kwargs", kwargs_file).kwargs
  
  # Set up our custom exception handlers
  if kwargs['debug']:
    sys.excepthook = ExceptionHookPDB
  else:
    sys.excepthook = ExceptionHook
  
  # Compute
  data = Compute(EPIC, **kwargs)
  
  # Plot the output
  Plot(data)

def RunCampaign(campaign, nodes = 5, ppn = 12, walltime = 100, 
                email = 'rodluger@gmail.com', 
                kwargs_file = os.path.join(EVEREST_ROOT, 'scripts', 'kwargs.py')):
  '''
  Submits a cluster job to compute and plot data for all targets in a given campaign.
  
  '''
  
  # Get all star IDs for this campaign
  stars = GetK2Stars()[campaign]
  nstars = len(stars)
  
  # Download the TPF data for each one
  for i, EPIC in enumerate(stars):
    print("Downloading data for EPIC %d (%d/%d)..." % (EPIC, i + 1, nstars))
    if not os.path.exists(os.path.join(KPLR_ROOT, 'data', 'everest', 
                          str(EPIC), str(EPIC) + '.npz')):
      try:
        GetK2Data(EPIC)
      except:
        # Some targets could be corrupted
        continue
        
  # Submit the cluster job      
  pbsfile = os.path.join(EVEREST_ROOT, 'everest', 'run.pbs')
  str_n = 'nodes=%d:ppn=%d,feature=%dcore' % (nodes, ppn, ppn)
  str_w = 'walltime=%d:00:00' % walltime
  str_v = 'EVEREST_ROOT=%s,NODES=%d,KWARGS_FILE=%s,CAMPAIGN=%d' % (EVEREST_ROOT, 
          nodes, os.path.abspath(kwargs_file), campaign)
  str_out = os.path.join(EVEREST_ROOT, 'C%02d.log' % campaign)
  qsub_args = ['qsub', pbsfile, 
               '-v', str_v, 
               '-o', str_out,
               '-j', 'oe', 
               '-N', 'C%02d' % campaign, 
               '-l', str_n,
               '-l', str_w,
               '-M', email,
               '-m', 'a']
            
  # Now we submit the job
  print("Submitting the job...")
  subprocess.call(qsub_args)

def _ComputeAndPlot(EPIC, **kwargs):
  '''
  Wrapper around ``Compute()`` and ``Plot()``.
  
  '''
  
  data = Compute(EPIC, **kwargs)
  Plot(data)
  return True

def _RunCampaign(campaign, kwargs_file):
  '''
  The actual function that runs a given campaign; this must
  be called from ``run.pbs``.
  
  '''
  
  # Set up our custom exception handler
  sys.excepthook = ExceptionHook
  
  # Initialize our multiprocessing pool
  with Pool() as pool:

    # Get the kwargs
    kwargs = imp.load_source("kwargs", kwargs_file).kwargs
  
    # Get all the stars
    stars = GetK2Stars()[campaign]
    nstars = len(stars)
  
    # Compute and plot
    C = FunctionWrapper(_ComputeAndPlot, **kwargs)
    pool.map(C, stars)