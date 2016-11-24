#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`pld.py` - PLD models
-----------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .math import Interpolate
import os, sys
import numpy as np
from itertools import combinations_with_replacement as multichoose
import logging
log = logging.getLogger(__name__)

__all__ = ['rPLD', 'nPLD']

class rPLD(object):
  '''
  The standard PLD model.
  
  '''
  
  def X(self, i, j = slice(None, None, None)):
    '''
    Computes the design matrix at the given *PLD* order and the given indices. 
    The columns are the *PLD* vectors for the target at the
    corresponding order, computed as the product of the fractional pixel
    flux of all sets of :py:obj:`n` pixels, where :py:obj:`n` is the *PLD*
    order.
    
    '''

    if self.recursive:
      X1 = self.fpix[j] / self.flux[j].reshape(-1, 1)
    else:
      X1 = self.fpix[j] / self.fraw[j].reshape(-1, 1)
    return np.product(list(multichoose(X1.T, i + 1)), axis = 1).T
  
  def _setup(self, **kwargs):
    '''
    This is called during production de-trending, prior to
    calling the :py:pbj:`Detrender.run()` method.
    
    '''
    
    pass
    
class nPLD(object):
  '''
  The "neighboring stars" *PLD* model. This model uses the *PLD* vectors of neighboring
  stars to help in the de-trending and can lead to increased performance over the regular
  :py:class:`rPLD` model, particularly for dimmer stars.
  
  '''
        
  def X(self, i, j = slice(None, None, None)):
    '''
    Computes the design matrix at the given *PLD* order and the given indices. 
    The columns are the *PLD* vectors for the target at the
    corresponding order, computed as the product of the fractional pixel
    flux of all sets of :py:obj:`n` pixels, where :py:obj:`n` is the *PLD*
    order. At the end of each matrix, columns corresponding to the neighbor
    star *PLD* signals are appended. Note that for both speed and memory
    reasons, cross terms are **not** computed for the neighboring stars.
    
    '''
      
    if self.recursive:
      X1 = self.fpix[j] / self.flux[j].reshape(-1, 1)
    else:
      X1 = self.fpix[j] / self.fraw[j].reshape(-1, 1)
    
    X = np.product(list(multichoose(X1.T, i + 1)), axis = 1).T
    return np.hstack([X, self.X1N[j] ** (i + 1)])
  
  def _setup(self, **kwargs):
    '''
    This is called during production de-trending, prior to
    calling the :py:pbj:`Detrender.run()` method.
    
    :param tuple cdpp_range:  If :py:obj:`parent_model` is set, neighbors are selected only if \
                              their de-trended CDPPs fall within this range. Default `None`
    :param tuple mag_range:   Only select neighbors whose magnitudes are within this range. \
                              Default (11., 13.) 
    :param int neighbors:     The number of neighboring stars to use in the de-trending. The \
                              higher this number, the more signals there are and hence the more \
                              de-trending information there is. However, the neighboring star \
                              signals are regularized together with the target's signals, so adding \
                              too many neighbors will inevitably reduce the contribution of the \
                              target's own signals, which may reduce performance. Default `10`
    :param str parent_model:  By default, :py:class:`nPLD` is run in stand-alone mode. The neighbor \
                              signals are computed directly from their TPFs, so there is no need to \
                              have run *PLD* on them beforehand. However, if :py:obj:`parent_model` \
                              is set, :py:class:`nPLD` will use information from the \
                              :py:obj:`parent_model` model of each neighboring star when de-trending. \
                              This is particularly useful for identifying outliers in the neighbor \
                              signals and preventing them from polluting the current target. Setting \
                              :py:obj:`parent_model` to :py:class:`rPLD`, for instance, will use the \
                              outlier information in the :py:class:`rPLD` model of the neighbors \
                              (this must have been run ahead of time). Note, however, that tests with \
                              *K2* data show that including outliers in the neighbor signals actually \
                              *improves* the performance, since many of these outliers are associated \
                              with events such as thruster firings and are present in all light curves, \
                              and therefore *help* in the de-trending. Default `None`
    
    '''
    
    # Get neighbors
    self.parent_model = kwargs.get('parent_model', None)
    num_neighbors = kwargs.get('neighbors', 10)
    self.neighbors = self._mission.GetNeighbors(self.ID, 
                                   cadence = self.cadence,
                                   model = self.parent_model,
                                   neighbors = num_neighbors, 
                                   mag_range = kwargs.get('mag_range', (11., 13.)), 
                                   cdpp_range = kwargs.get('cdpp_range', None),
                                   aperture_name = self.aperture_name)
    if len(self.neighbors):
      if len(self.neighbors) < num_neighbors:
        log.warn("%d neighbors requested, but only %d found." % (num_neighbors, len(self.neighbors)))
    else:
      raise Exception("No neighbors found! Aborting.")
    
    for neighbor in self.neighbors:
      log.info("Loading data for neighboring target %d..." % neighbor)
      if self.parent_model is not None and self.cadence == 'lc':
        # We load the `parent` model. The advantage here is that outliers have
        # properly been identified and masked. I haven't tested this on short
        # cadence data, so I'm going to just forbid it...
        data = eval(self.parent_model)(neighbor, mission = self.mission, is_parent = True)
      else:
        # We load the data straight from the TPF. Much quicker, since no model must
        # be run in advance. Downside is we don't know where the outliers are. But based
        # on tests with K2 data, the de-trending is actually *better* if the outliers are
        # included! These are mostly thruster fire events and other artifacts common to
        # all the stars, so it makes sense that we might want to keep them in the design
        # matrix.
        data = self._mission.GetData(neighbor, season = self.season, clobber = self.clobber_tpf, 
                             cadence = self.cadence,
                             aperture_name = self.aperture_name, 
                             saturated_aperture_name = self.saturated_aperture_name, 
                             max_pixels = self.max_pixels,
                             saturation_tolerance = self.saturation_tolerance)
        data.mask = np.array(list(set(np.concatenate([data.badmask, data.nanmask]))), dtype = int)
        data.fraw = np.sum(data.fpix, axis = 1)
      
      # Compute the linear PLD vectors and interpolate over outliers, NaNs and bad timestamps
      X1 = data.fpix / data.fraw.reshape(-1, 1)
      X1 = Interpolate(data.time, data.mask, X1)
      if self.X1N is None:
        self.X1N = np.array(X1)
      else:
        self.X1N = np.hstack([self.X1N, X1])
      del X1
      del data