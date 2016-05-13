#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, unicode_literals

# MPL backend: force Agg for all Everest modules if running on a Linux machine
# In order for this to work, ``everest`` must be imported first!
# If on MacOS, try the Qt4Agg backend before the MacOSX backend, which is
# quite bugged (at least on my Mac!). In particular, it complains when trying
# to save JPEGs.
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

# Add our submodules to the PATH
import os, sys
for submodule in ['kplr', 'pysyzygy']:
  sys.path.insert(1, os.path.join(os.path.dirname(os.path.dirname(
                                  os.path.abspath(__file__))), submodule))

try:
  import pysyzygy
except Exception as e:
  if str(e).startswith("Can't find ``transitlib.so``"):
    raise Exception("Please compile ``pysyzygy`` by running ``make`` in '/everest/pysyzygy'.")

# Import modules
from . import compute, data, detrend, gp, kernels, pool, sources, transit, utils
from .data import GetK2Data, GetK2Planets, GetK2EBs, GetK2Stars
from .pool import Pool
from .compute import Compute
from .run import DownloadCampaign, DownloadInjections, RunSingle, RunCampaign, RunCandidates, RunInjections