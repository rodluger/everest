#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, unicode_literals
import functools
from . import k2

#: The supported missions
Missions = {'k2': k2}
 
class MissionWrapper(object):
  '''
  A generic wrapper for mission-specific routines.
  
  '''
  
  def __init__(self, name):
    '''
    
    '''
    
    self.name = name

  def __call__(self, *args, mission = 'k2', **kwargs):
    '''
    
    '''
    
    if mission in Missions:
      return functools.reduce(getattr, [Missions[mission]] + self.name.split('.'))(*args, **kwargs)
    else:
      raise ValueError('Mission %s not supported.' % mission)