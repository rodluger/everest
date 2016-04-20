#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
transit.py
----------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np
import matplotlib.pyplot as pl
import pysyzygy as ps
from scipy.optimize import fmin
import logging
log = logging.getLogger(__name__)

def Get_RpRs(d, **kwargs):
  '''
  Returns the value of Rp/Rs for a given depth, given
  the ``pysyzygy`` transit **kwargs.
  
  '''
  
  def Depth(RpRs, **kwargs):
    return 1 - ps.Transit(RpRs = RpRs, **kwargs)([kwargs.get('t0', 0.)])

  def DiffSq(r):
    return 1.e10 * (d - Depth(r, **kwargs)) ** 2  
  
  return fmin(DiffSq, [np.sqrt(d)], disp = False)

def Get_rhos(dur, **kwargs):
  '''
  Returns the value of rhos for a given transit duration, given
  the ``pysyzygy`` transit **kwargs.
  
  '''
  
  assert dur >= 0.05 and dur <= 0.5, "Invalid value for the duration."
  
  def Dur(rhos, **kwargs):
    t0 = kwargs.get('t0', 0.)
    time = np.linspace(-0.5, 0.5, 1000)
    try:
      t = time[np.where(ps.Transit(rhos = rhos, **kwargs)(time) < 1)]
    except:
      return 0.
    return t[-1] - t[0]

  def DiffSq(rhos):
    return (dur - Dur(rhos, **kwargs)) ** 2  
  
  return fmin(DiffSq, [0.2], disp = False)

def Transit(time, t0 = 0., dur = 0.1, per = 3.56789, depth = 0.001, **kwargs):
  '''
  A Mandel-Agol transit model, but with the depth and the duration
  as primary input variables.
  
  '''
  
  # Remove extra kwargs
  for k in ['mask', 'everest', 'k2sff']:
    kwargs.pop(k, None)
  
  # Note that rhos can affect RpRs, so we should really do this iteratively,
  # but the effect is pretty negligible!
  RpRs = Get_RpRs(depth, t0 = t0, per = per, **kwargs)
  rhos = Get_rhos(dur, t0 = t0, per = per, **kwargs)
  return ps.Transit(t0 = t0, per = per, RpRs = RpRs, rhos = rhos, **kwargs)(time)

def TopHat(time, t0 = 0., per = 3.56789, dur = 0.1, depth = 0.001, **kwargs):
  '''
  A simple top-hat transit model.
  
  '''
  
  inds = np.where( np.abs((time - t0 - per / 2.) % per - per / 2.) < dur / 2. )
  model = np.ones_like(time)
  model[inds] = 1. - depth
  return model