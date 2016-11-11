#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`transit.py` - Transit models
-------------------------------------
These are routines used to generate a transit model, primarily for
transit injection/recovery tests. These are wrappers around
:py:func:`pysyzygy.Transit`, with the added feature that
the transit :py:obj:`depth` and the transit :py:obj:`duration` can be specified
as input variables (as opposed to the planet-star radius ratio
and the stellar density, which :py:mod:`pysyzygy` expects).

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
  Returns the value of the planet radius over the stellar radius for a given depth :py:obj:`d`, given
  the :py:class:`everest.pysyzygy` transit :py:obj:`kwargs`.
  
  '''
  
  def Depth(RpRs, **kwargs):
    return 1 - ps.Transit(RpRs = RpRs, **kwargs)([kwargs.get('t0', 0.)])

  def DiffSq(r):
    return 1.e10 * (d - Depth(r, **kwargs)) ** 2  
  
  return fmin(DiffSq, [np.sqrt(d)], disp = False)

def Get_rhos(dur, **kwargs):
  '''
  Returns the value of the stellar density for a given transit duration :py:obj:`dur`, given
  the :py:class:`everest.pysyzygy` transit :py:obj:`kwargs`.
  
  '''
  
  assert dur >= 0.05 and dur <= 0.5, "Invalid value for the duration."
  
  def Dur(rhos, **kwargs):
    t0 = kwargs.get('t0', 0.)
    time = np.linspace(t0 - 0.5, t0 + 0.5, 1000)
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
  A `Mandel-Agol <http://adsabs.harvard.edu/abs/2002ApJ...580L.171M>`_ transit model, 
  but with the depth and the duration
  as primary input variables.
  
  :param numpy.ndarray time: The time array
  :param float t0: The time of first transit in units of :py:obj:`BJD` - 2454833.
  :param float dur: The transit duration in days. Don't go too crazy on this one -- very small \
                    or very large values will break the inverter. Default 0.1
  :param float per: The orbital period in days. Default 3.56789
  :param float depth: The fractional transit depth. Default 0.001
  :param dict kwargs: Any additional keyword arguments, passed directly to :py:func:`everest.pysyzygy.Transit`
  :returns tmod: The transit model evaluated at the same times as the :py:obj:`time` array
  
  '''
  
  # Note that rhos can affect RpRs, so we should really do this iteratively,
  # but the effect is pretty negligible!
  RpRs = Get_RpRs(depth, t0 = t0, per = per, **kwargs)
  rhos = Get_rhos(dur, t0 = t0, per = per, **kwargs)
  return ps.Transit(t0 = t0, per = per, RpRs = RpRs, rhos = rhos, **kwargs)(time)