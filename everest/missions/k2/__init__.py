#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, unicode_literals
from .k2 import *
from . import aux, pbs

#: The string that identifies individual targets for this mission
IDSTRING = 'EPIC'
#: The character abbreviation of the name given to an observing "season" for this mission
SEASONCHAR = 'C'
#: The string representing the filter/band used in the mission 
MAGSTRING = r'K$_\mathrm{p}$'