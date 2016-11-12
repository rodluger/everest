#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`pbs.py` - PBS cluster routines
---------------------------------------

Routines to run :py:mod:`everest` in batch mode on a PBS cluster.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .missions import Missions

def Download(mission = 'k2', **kwargs):
  '''
  
  '''
  
  if mission in Missions:
    return Missions[mission].pbs.Download(**kwargs)
  else:
    raise ValueError('Mission %s not supported.' % mission)

def Run(mission = 'k2', **kwargs):
  '''
  
  '''
  
  if mission in Missions:
    return Missions[mission].pbs.Run(**kwargs)
  else:
    raise ValueError('Mission %s not supported.' % mission)

def Status(mission = 'k2', **kwargs):
  '''
  
  '''
  
  if mission in Missions:
    return Missions[mission].pbs.Status(**kwargs)
  else:
    raise ValueError('Mission %s not supported.' % mission)
