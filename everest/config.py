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

__all__ = ["EVEREST_DAT", "EVEREST_SRC", "MAST_ROOT", "KWARGS_PY"]

# Directories
EVEREST_DAT = os.path.expanduser(os.environ.get("EVEREST_DATA_DIR", os.path.join("~", ".everest")))                               
EVEREST_SRC = os.path.dirname(os.path.abspath(__file__))

# MAST url
EVEREST_FITS = os.environ.get('EVEREST_FITS', None)
MAST_ROOT = 'https://archive.stsci.edu/missions/hlsp/everest/v1/'

# Default kwargs file
KWARGS_PY = \
"""
#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
kwargs.py
---------

This file contains the kwargs passed to :py:func:`everest.compute.Compute()` when
running a cluster job with :py:func:`everest.run.RunCampaign()`.
See the :py:func:`everest.compute.Compute()` documentation for info on each
of these values.
   
'''

from __future__ import division, print_function, absolute_import, unicode_literals

kwargs = dict()
"""