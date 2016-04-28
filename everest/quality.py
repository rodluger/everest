#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
quality.py
----------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np

def Crowding(k0, source, ny, nx, contour):
  '''
  Assess the severity of the aperture crowding.
  TODO: Do this more rigorously!
  
  '''

  # Get the distance from a source to the edge of the aperture
  # Not terribly robust, but sufficient for an estimate
  dist = np.inf
  m, n = int(np.round(source.y - source.y0)), int(np.round(source.x - source.x0))
  for dm in range(-3,4):
    for dn in range(-3,4):
      if m + dm >= 0 and m + dm < ny and n + dn >= 0 and n + dn < nx:
        if contour[m + dm, n + dn]:
          d = np.sqrt(dm ** 2 + dn ** 2)
          if d < dist:
            dist = d
  
  # The difference in magnitudes      
  dk = source.kepmag - k0
  
  # Compute the metric [0-5]
  if np.isnan(dk):
    return np.nan
  elif np.isinf(dist):
    return 0
  elif dist == 0 and dk < 0:
    return 5
  elif dist == 0 and dk < 3:
    return 4
  elif dist == 0 and dk < 5:
    return 3
  elif dist == 0 and dk > 5:
    return 1
  elif dist <= 1. and dk < 0:
    return 4
  elif dist <= 1. and dk < 3:
    return 2
  elif dist <= 1. and dk < 5:
    return 1
  elif dist <= 2. and dk < 0:
    return 3
  elif dist <= 2. and dk < 5:
    return 2
  elif dist <= 3. and dk < 0:
    return 2
  else:
    return 0

def Saturation(fpix):
  '''
  Assess the saturation severity. 
  TODO: Do this more rigorously!
  
  '''
  
  maxmed = np.median(fpix.T[0])
  if maxmed > 180000:
    satsev = 5
  elif maxmed > 170000:
    satsev = 4
  elif maxmed > 160000:
    satsev = 3
  elif maxmed > 150000:
    satsev = 2
  elif maxmed > 140000:
    satsev = 1
  else:
    satsev = 0
  return satsev

def Autocorrelation(kchisq):
  '''
  Gauge the acor fitting severity.
  
  '''
  
  if kchisq < 1.:
    acorsev = 0
  elif kchisq < 30.:
    acorsev = 1
  elif kchisq < 70.:
    acorsev = 2
  elif kchisq < 100.:
    acorsev = 3
  elif kchisq < 150.:
    acorsev = 4
  else:
    acorsev = 5
  return acorsev