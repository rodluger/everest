#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, unicode_literals

# Version number
__version__ = "2.0"

# Was everest imported from setup.py?
try:
  __EVEREST_SETUP__
except NameError:
  __EVEREST_SETUP__ = False

if __EVEREST_SETUP__:
  # Set up the individual missions
  from .missions import Missions
  for mission in Missions:
    Missions[mission].Setup()
else:
  # This is a regular everest run
  from . import config, data, dvs, gp, math, missions, model, pbs, pool, transit, utils
  from .model import *
  from .missions import *
  from .data import Statistics
  from .pbs import Download, Run, Status
  from .transit import Transit