#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`config.py` - Misc settings
-----------------------------------

Stores information on the location of the :py:mod:`everest` source code,
the data files, and the MAST url.

'''

from __future__ import division, print_function, absolute_import, \
     unicode_literals
from . import __version__ as EVEREST_VERSION
import os
import six
from six.moves import urllib
import logging
log = logging.getLogger(__name__)

#: major.minor version number
EVEREST_MAJOR_MINOR = ".".join(EVEREST_VERSION.split(".")[:-1])

#: Is this the development (pre-publication) version?
EVEREST_DEV = int(os.environ.get('EVEREST2_DEV', 0))
if EVEREST_DEV:

    # Dev version hack: enforce a non-ui backend
    import platform
    import matplotlib as mpl
    if platform.system() == "Linux":
        mpl.use("Agg", warn=False)
    else:
        # Dev version hack: custom font
        mpl.rc('font', family='serif')
        mpl.rc('font', serif='Palatino Linotype')

#: The :py:mod:`everest` data directory
EVEREST_DAT = os.path.expanduser(os.environ.get("EVEREST2_DATA_DIR",
                                                os.path.join("~",
                                                             ".everest2")))
#: The :py:mod:`everest` source code directory
EVEREST_SRC = os.path.dirname(os.path.abspath(__file__))
#: The ``user@server:/path`` scp argument for accessing the FITS files
#: (pre-publication only)
EVEREST_FITS = os.environ.get('EVEREST2_FITS', None)
#: The directory containing the Kepler PRF files
KEPPRF_DIR = os.path.expanduser(os.environ.get(
    "KEPPRF_DIR", os.path.join("~", "src", "KeplerPRF")))

#: The MAST url where the light curves are published
MAST_ROOT = 'https://archive.stsci.edu/hlsps/everest/v2/'

#: Everest quality bit: masked because a Kepler flag was raised
QUALITY_BAD = 23
#: Everest quality bit: masked because data was NaN
QUALITY_NAN = 24
#: Everest quality bit: masked because data was an outlier
QUALITY_OUT = 25
#: Everest quality bit: masked in the original model (for recursive PLD only)
QUALITY_REC = 26
#: Everest quality bit: masked transit cadence
QUALITY_TRN = 27
