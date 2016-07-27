#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys

# Version number
__version__ = "1.0"
__develop__ = False

# Was everest imported from setup.py?
try:
  __EVEREST_SETUP__
except NameError:
  __EVEREST_SETUP__ = False

# This is a regular everest run
if not __EVEREST_SETUP__:
  
  # Matplotlib backend hack (dev only)
  if __develop__:
    import platform
    if platform.system() == "Linux":
      import matplotlib as mpl
      mpl.use("Agg", warn=False)
    elif platform.system() == "Darwin":
      import matplotlib as mpl
      try:
        mpl.use("Qt4Agg", warn=False)
      except:
        pass

  # Make sure we can import our submodules
  import pysyzygy
  import k2plr

  # Import modules
  from . import config, compute, crowding, data, detrend, fits, gp, kernels, pool, \
                sources, transit, usertools, utils
  from .data import GetK2Data, GetK2Planets, GetK2EBs, GetK2Stars, Progress, Campaign
  from .pool import Pool
  from .compute import Compute
  from .run import DownloadCampaign, DownloadInjections, RunSingle, RunCampaign, \
                   RunCandidates, RunInjections, RunFITS
  from .fits import MakeFITS
  from .usertools import Everest
  from matplotlib.pyplot import show
  
  # Create the data directories if they don't exist
  if not os.path.exists(os.path.join(config.EVEREST_DAT, 'output')):
    os.makedirs(os.path.join(config.EVEREST_DAT, 'output'))
  if not os.path.exists(os.path.join(config.EVEREST_DAT, 'fits')):
    os.makedirs(os.path.join(config.EVEREST_DAT, 'fits'))
  if not os.path.exists(os.path.join(config.EVEREST_DAT, 'kwargs.py')):
    with open(os.path.join(config.EVEREST_DAT, 'kwargs.py'), 'w') as f:
      f.write(config.KWARGS_PY)