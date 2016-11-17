#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`config.py` - Misc settings
-----------------------------------

Stores information on the location of the :py:mod:`everest` source code,
the data files, and the MAST url.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os

#: Is this the development (pre-publication) version?
EVEREST_DEV = os.environ.get('EVEREST2_DEV', 0)
if EVEREST_DEV:
  # Dev version hack: enforce a non-ui backend
  import platform
  if platform.system() == "Linux":
    import matplotlib as mpl
    mpl.use("Agg", warn=False)

#: The :py:mod:`everest` data directory
EVEREST_DAT = os.path.expanduser(os.environ.get("EVEREST2_DATA_DIR", os.path.join("~", ".everest2")))                               
#: The :py:mod:`everest` source code directory
EVEREST_SRC = os.path.dirname(os.path.abspath(__file__))
#: The MAST url where the light curves are published
MAST_ROOT = 'https://archive.stsci.edu/missions/hlsp/everest/'

#: Everest quality bit: masked because a Kepler flag was raised
QUALITY_BAD = 23
#: Everest quality bit: masked because data was NaN
QUALITY_NAN = 24
#: Everest quality bit: masked because data was an outlier
QUALITY_OUT = 25