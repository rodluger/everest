#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`config.py` - $PATH settings
------------------------------------

Stores information on the location of the :py:mod:`everest` source code,
the data files, and the MAST url.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os

# Dev version hacks
EVEREST_DEV = os.environ.get('EVEREST2_DEV', 0)
if EVEREST_DEV:
  import platform
  if platform.system() == "Linux":
    import matplotlib as mpl
    mpl.use("Agg", warn=False)

# Directories
EVEREST_DAT = os.path.expanduser(os.environ.get("EVEREST2_DATA_DIR", os.path.join("~", ".everest2")))                               
EVEREST_SRC = os.path.dirname(os.path.abspath(__file__))

# MAST url
EVEREST_FITS = os.environ.get('EVEREST2_FITS', None)
MAST_ROOT = 'https://archive.stsci.edu/missions/hlsp/everest/'