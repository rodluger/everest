#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
run.py
------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .plot import Plot
from .data import GetK2Stars, GetK2Data, GetK2Planets, GetK2InjectionTestStars
from .compute import Compute
from .utils import ExceptionHook, ExceptionHookPDB, FunctionWrapper
from .pool import Pool
import os
EVEREST_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
from kplr.config import KPLR_ROOT
import sys
import subprocess
import imp
import time

def DownloadCampaign(campaign, queue = 'build', email = 'rodluger@gmail.com', walltime = 8):
  '''
  Submits a cluster job to the build queue to download all TPFs for a given
  campaign.
  
  '''
          
  # Submit the cluster job      
  pbsfile = os.path.join(EVEREST_ROOT, 'everest', 'download.pbs')
  str_w = 'walltime=%d:00:00' % walltime
  str_v = 'EVEREST_ROOT=%s,CAMPAIGN=%d' % (EVEREST_ROOT, campaign)
  str_out = os.path.join(EVEREST_ROOT, 'DOWNLOAD_C%02d.log' % campaign)
  qsub_args = ['qsub', pbsfile, 
               '-q', queue,
               '-v', str_v, 
               '-o', str_out,
               '-j', 'oe', 
               '-N', 'DOWNLOAD_C%02d' % campaign,
               '-l', str_w,
               '-M', email,
               '-m', 'ae']
            
  # Now we submit the job
  print("Submitting the job...")
  subprocess.call(qsub_args)

def Run(EPIC, **kwargs):
  '''
  Wrapper around ``Compute()`` and ``Plot()``.
  
  '''
  
  data = Compute(EPIC, **kwargs)
  if data is not None:
    Plot(data)
    return True
  else:
    return False

def RunSingle(EPIC, debug = False, kwargs_file = os.path.join(EVEREST_ROOT, 'scripts', 'kwargs.py')):
  '''
  Compute and plot data for a given target.
  
  '''
  
  # Get the kwargs
  kwargs = imp.load_source("kwargs", kwargs_file).kwargs
  
  # Set up our custom exception handlers
  if debug:
    sys.excepthook = ExceptionHookPDB
  else:
    sys.excepthook = ExceptionHook
  
  # Run
  Run(EPIC, **kwargs)

def RunInjections(depth = 0.01, mask = False, queue = 'vsm',
                  nodes = 5, ppn = 12, walltime = 100, 
                  email = 'rodluger@gmail.com', 
                  kwargs_file = os.path.join(EVEREST_ROOT, 'scripts', 'kwargs.py')):
  '''
  Submits a cluster job to compute and plot data for a sample
  of targets chosen for transit injection tests.
  
  '''
  
  # Submit the cluster job   
  name = 'inject_%.4f%s' % (depth, ('m' if mask else 'u'))   
  pbsfile = os.path.join(EVEREST_ROOT, 'everest', 'runinjections.pbs')
  str_n = 'nodes=%d:ppn=%d,feature=%dcore' % (nodes, ppn, ppn)
  str_w = 'walltime=%d:00:00' % walltime
  str_v = 'EVEREST_ROOT=%s,NODES=%d,MASK=%d,DEPTH=%0.4f,KWARGS_FILE=%s' % (EVEREST_ROOT, 
          nodes, int(mask), depth, os.path.abspath(kwargs_file))
  str_out = os.path.join(EVEREST_ROOT, '%s.log' % name)
  qsub_args = ['qsub', pbsfile, 
               '-q', queue,
               '-v', str_v, 
               '-o', str_out,
               '-j', 'oe', 
               '-N', name, 
               '-l', str_n,
               '-l', str_w,
               '-M', email,
               '-m', 'ae']
            
  # Now we submit the job
  print("Submitting the job...")
  subprocess.call(qsub_args)

def RunCandidates(nodes = 5, ppn = 12, walltime = 100, queue = 'vsm',
                  email = 'rodluger@gmail.com', 
                  kwargs_file = os.path.join(EVEREST_ROOT, 'scripts', 'kwargs.py')):
  '''
  Submits a cluster job to compute and plot data for all targets hosting
  confirmed planets or planet candidates.
  
  '''
          
  # Submit the cluster job      
  pbsfile = os.path.join(EVEREST_ROOT, 'everest', 'runcandidates.pbs')
  str_n = 'nodes=%d:ppn=%d,feature=%dcore' % (nodes, ppn, ppn)
  str_w = 'walltime=%d:00:00' % walltime
  str_v = 'EVEREST_ROOT=%s,NODES=%d,KWARGS_FILE=%s' % (EVEREST_ROOT, 
          nodes, os.path.abspath(kwargs_file))
  str_out = os.path.join(EVEREST_ROOT, 'candidates.log')
  qsub_args = ['qsub', pbsfile, 
               '-q', queue, 
               '-v', str_v, 
               '-o', str_out,
               '-j', 'oe', 
               '-N', 'candidates', 
               '-l', str_n,
               '-l', str_w,
               '-M', email,
               '-m', 'ae']
            
  # Now we submit the job
  print("Submitting the job...")
  subprocess.call(qsub_args)

def RunCampaign(campaign, nodes = 5, ppn = 12, walltime = 100, 
                email = 'rodluger@gmail.com', queue = 'vsm',
                kwargs_file = os.path.join(EVEREST_ROOT, 'scripts', 'kwargs.py')):
  '''
  Submits a cluster job to compute and plot data for all targets in a given campaign.
  
  '''
          
  # Submit the cluster job      
  pbsfile = os.path.join(EVEREST_ROOT, 'everest', 'runcampaign.pbs')
  str_n = 'nodes=%d:ppn=%d,feature=%dcore' % (nodes, ppn, ppn)
  str_w = 'walltime=%d:00:00' % walltime
  str_v = 'EVEREST_ROOT=%s,NODES=%d,KWARGS_FILE=%s,CAMPAIGN=%d' % (EVEREST_ROOT, 
          nodes, os.path.abspath(kwargs_file), campaign)
  str_out = os.path.join(EVEREST_ROOT, 'C%02d.log' % campaign)
  qsub_args = ['qsub', pbsfile, 
               '-q', queue, 
               '-v', str_v, 
               '-o', str_out,
               '-j', 'oe', 
               '-N', 'C%02d' % campaign, 
               '-l', str_n,
               '-l', str_w,
               '-M', email,
               '-m', 'ae']
            
  # Now we submit the job
  print("Submitting the job...")
  subprocess.call(qsub_args)

# ---- PBS routines ---- #

def _RunCandidates(kwargs_file):
  '''
  The actual function that runs all candidates; this must
  be called from ``runcandidates.pbs``.
  
  '''
  
  # Set up our custom exception handler
  sys.excepthook = ExceptionHook
  
  # Initialize our multiprocessing pool
  with Pool() as pool:

    # Get the kwargs
    kwargs = imp.load_source("kwargs", kwargs_file).kwargs
  
    # Override ``mask_candidates``
    kwargs['mask_candidates'] = True
  
    # Get all the stars
    stars = [int(p.epic_name[5:]) for p in GetK2Planets()]
    new = [int(f[:9]) for f in os.listdir(os.path.join(EVEREST_ROOT, 'new')) if f.endswith('.npz')]
    stars = list(set(stars + new))
  
    # Compute and plot
    C = FunctionWrapper(Run, **kwargs)
    pool.map(C, stars)

def _RunInjections(kwargs_file, depth, mask):
  '''
  The actual function that runs injection tests; this must
  be called from ``runinjections.pbs``.
  
  '''
  
  # Set up our custom exception handler
  sys.excepthook = ExceptionHook
  
  # Initialize our multiprocessing pool
  with Pool() as pool:

    # Get the kwargs
    kwargs = imp.load_source("kwargs", kwargs_file).kwargs
    
    # Get all the stars
    stars = [int(s) for s in GetK2InjectionTestStars()]
    
    # Override ``inject`` and ``run_name``
    kwargs['inject'] = dict(t0 = 0., 
                            per = 3.56789, 
                            depth = depth, 
                            dur = 0.1,
                            mask = mask)
    kwargs['run_name'] = 'inject_%.4f%s' % (depth, ('m' if mask else 'u'))

    # Compute and plot
    C = FunctionWrapper(Run, **kwargs)
    pool.map(C, stars)

def _RunCampaign(campaign, kwargs_file):
  '''
  The actual function that runs a given campaign; this must
  be called from ``runcampaign.pbs``.
  
  '''
  
  # Set up our custom exception handler
  sys.excepthook = ExceptionHook
  
  # Initialize our multiprocessing pool
  with Pool() as pool:

    # Get the kwargs
    kwargs = imp.load_source("kwargs", kwargs_file).kwargs
        
    # Get all the stars
    stars = GetK2Stars()[campaign]
  
    # Compute and plot
    C = FunctionWrapper(Run, **kwargs)
    pool.map(C, stars)

def _DownloadCampaign(campaign):
  '''
  Download all stars from a given campaign. This is
  called from ``download.pbs``
  
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