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
  from . import missions
  for mission in missions.Missions:
    getattr(missions, mission).Setup()
else:
  # This is a regular everest run
  from . import basecamp, config, detrender, dvs, fits, gp, inject, math, \
                missions, pld, pool, transit, user, utils
  from .detrender import *
  from .inject import *
  from .missions import *
  from .transit import Transit
  from .user import Everest