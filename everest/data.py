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
  A wrapper around the mission-specific :py:obj:`Season` function. Returns
  the number of the current observational season/quarter/campaign for the target.
  
  :param id: The target ID
  :param str mission: The mission name. Default `k2`

  '''
  
  if mission in Missions:
    return Missions[mission].Season(id, **kwargs)
  else:
    raise ValueError('Mission %s not supported.' % mission)
  
def GetData(id, mission = 'k2', **kwargs):
  '''
  A wrapper around the mission-specific :py:obj:`GetData` function. Returns
  a :py:class:`utils.DataContainer` instance with the raw light curve data.
  
  :param id: The target ID
  :param str mission: The mission name. Default `k2`  
  
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

def MakeFITS(model, mission = 'k2', **kwargs):
  '''
  
  '''
  
  if mission in Missions:
    return Missions[mission].MakeFITS(model, **kwargs)
  else:
    raise ValueError('Mission %s not supported.' % mission)