#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`__init__.py` - Initialization
--------------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .kepler import *

#: The string that identifies individual targets for this mission
IDSTRING = 'KIC'
#: The character abbreviation of the name given to an observing "season" for this mission
SEASONCHAR = 'Q'
#: The string representing the filter/band used in the mission 
MAGSTRING = r'K$_\mathrm{p}$'
#: The time units for the mission
TIMEUNITS = 'BJD'
#: The currently supported seasons
SEASONS = []
#: Returns :py:obj:`True` if argument is a valid `Kepler` target identifier (necessary but not sufficient)
ISTARGET = lambda x: ((type(x) is int) and (x > 1e6) and (x < 1e7))
#: The published light curve CSV file header
CSVHEADER = \
'''KIC %07d
============

'''