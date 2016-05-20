#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`config.py` - PATH settings
-----------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os

__all__ = ["EVEREST_DAT", "EVEREST_SRC", "MAST_ROOT"]

EVEREST_DAT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                                 
EVEREST_SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAST_ROOT = 'https://archive.stsci.edu/missions/hlsp/everest/'

# In the future...
# EVEREST_DAT = os.path.expanduser(os.environ.get("EVEREST_DATA_DIR", os.path.join("~", ".everest")))