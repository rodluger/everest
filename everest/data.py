#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`data.py` - Data access routines
----------------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .missions import Missions, MissionWrapper

#:
Breakpoint = MissionWrapper('Breakpoint')
#:
FITSFile = MissionWrapper('FITSFile')
#:
GetData = MissionWrapper('GetData')
#:
GetNeighbors = MissionWrapper('GetNeighbors')
#:
HasShortCadence = MissionWrapper('HasShortCadence')
#:
HDUCards = MissionWrapper('HDUCards')
#:
Season = MissionWrapper('Season')
#:
Statistics = MissionWrapper('Statistics')

