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
import urllib

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
#: The ``user@server:/path`` scp argument for accessing the FITS files (pre-publication only)
EVEREST_FITS = os.environ.get('EVEREST2_FITS', None)
#: The MAST url where the light curves are published
MAST_ROOT = 'https://archive.stsci.edu/missions/hlsp/everest/'

def MAST_VERSION():
  '''
  Returns the current :py:mod:`everest` version on MAST.
  
  '''
  
  url = MAST_ROOT + 'version.txt'
  r = urllib.request.Request(url)
  handler = urllib.request.urlopen(r)
  code = handler.getcode()
  if int(code) != 200:
    raise Exception("Error code {0} for URL '{1}'".format(code, url))
  data = handler.read().decode('utf-8')
  if data.endswith('\n'):
    data = data[:-1]

  return data

#: Everest quality bit: masked because a Kepler flag was raised
QUALITY_BAD = 23
#: Everest quality bit: masked because data was NaN
QUALITY_NAN = 24
#: Everest quality bit: masked because data was an outlier
QUALITY_OUT = 25