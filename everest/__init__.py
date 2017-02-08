#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, unicode_literals

# Version number
__version__ = "2.0"
__subversion__ = "2.0.3"

# Was everest imported from setup.py?
try:
  __EVEREST_SETUP__
except NameError:
  __EVEREST_SETUP__ = False

if not __EVEREST_SETUP__:
  # This is a regular everest run

  # Import all modules
  from . import config
  from . import utils
  from . import math
  from . import transit
  from . import pool
  from . import fits
  from . import dvs
  from . import gp
  from . import missions
  from . import basecamp
  from . import detrender
  from . import inject
  from . import user
  
  # Import the good stuff
  from .detrender import *
  from .inject import *
  from .missions import *
  from .transit import Transit
  from .user import Everest, DVS
