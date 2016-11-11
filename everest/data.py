#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`data.py` - Data access routines
----------------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .missions import Missions

def Season(id, mission = 'k2', **kwargs):
  '''
  
  '''
  
  if mission in Missions:
    return Missions[mission].Season(id, **kwargs)
  else:
    raise ValueError('Mission %s not supported.' % mission)
  
def GetData(id, mission = 'k2', **kwargs):
  '''
  
  '''
  
  if mission in Missions:
    return Missions[mission].GetData(id, **kwargs)
  else:
    raise ValueError('Mission %s not supported.' % mission)

def GetNeighbors(id, mission = 'k2', **kwargs):
  '''
  
  '''

  if mission in Missions:
    return Missions[mission].GetNeighbors(id, **kwargs)
  else:
    raise ValueError('Mission %s not supported.' % mission)

def Breakpoint(id, mission = 'k2', **kwargs):
  '''
  
  '''
  
  if mission in Missions:
    return Missions[mission].Breakpoint(id, **kwargs)
  else:
    raise ValueError('Mission %s not supported.' % mission)

def Statistics(mission = 'k2', **kwargs):
  '''
  
  '''
  
  if mission in Missions:
    return Missions[mission].Statistics(**kwargs)
  else:
    raise ValueError('Mission %s not supported.' % mission)

def HasShortCadence(id, mission = 'k2', **kwargs):
  '''
  
  '''
  
  if mission in Missions:
    return Missions[mission].HasShortCadence(id, **kwargs)
  else:
    raise ValueError('Mission %s not supported.' % mission)