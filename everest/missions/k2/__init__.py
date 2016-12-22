#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`__init__.py` - Initialization
--------------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .k2 import *
from .sysrem import GetCBVs
from . import aux, pbs, pipelines, sysrem
from .pbs import Download, Run, Status, Publish

#: The string that identifies individual targets for this mission
IDSTRING = 'EPIC'
#: The character abbreviation of the name given to an observing "season" for this mission
SEASONCHAR = 'C'
#: The string representing the filter/band used in the mission 
MAGSTRING = r'K$_\mathrm{p}$'
#: The time units for the mission
TIMEUNITS = 'BJD - 2454833'