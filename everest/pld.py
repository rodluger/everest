#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`pld.py` - PLD models
-----------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .data import GetData, GetNeighbors
from .math import Interpolate
import os, sys
import numpy as np
from itertools import combinations_with_replacement as multichoose
import logging
log = logging.getLogger(__name__)

__all__ = ['rPLDBase', 'nPLDBase']

class rPLDBase(object):
  '''
  
  '''
   
  def get_X(self):
    '''
    Computes the design matrix at the current *PLD* order and stores it in
    :py:obj:`self.X`. This is a list of matrices, one for each *PLD* order.
    The columns in each matrix are the *PLD* vectors for the target at the
    corresponding order, computed as the product of the fractional pixel
    flux of all sets of :py:obj:`n` pixels, where :py:obj:`n` is the *PLD*
    order.
    
    '''
      
    if not self.is_parent:
      log.info("Computing the design matrix...")
    if self.recursive:
      X1 = self.fpix / self.flux.reshape(-1, 1)
    else:
      X1 = self.fpix / self.fraw.reshape(-1, 1)
    for n in range(self.pld_order): 
      if (self._X[n] is None) and ((n == self.lam_idx) or (self.lam[0][n] is not None)):
        self._X[n] = np.product(list(multichoose(X1.T, n + 1)), axis = 1).T
  
  def get_sc_model(self, order, weights, inds = None):
    '''
    Computes the short cadence model. The PLD model is always computed from
    the long cadence data; the weights are then dotted with the short cadence
    design matrix to compute the short cadence model. In principle, one could
    compute the model from the short cadence data, but (a) this is likely to
    take a very long time and lead to memory errors, and (b) it is likely to
    perform poorly, since the SNR ratio of each data point in short cadence is
    very low, making it difficult for PLD to pick out the instrumental component
    of the signal. See Deming et al. (2015) for a discussion on how computing the
    model on binned (i.e., long cadence) data is ideal.
    
    .. note:: The code below uses a :py:obj:`for` loop to dot each signal with its \
              corresponding weight, which is inefficient. However, matrix operations \
              on the short cadence design matrix take a **huge** amount of memory \
              and usually cause the script to crash. By computing the model this way, \
              the design matrix is never actually stored in memory, but processed one \
              column at a time. It actually works surprisingly fast.
    
    '''
    
    if inds is None:
      inds = range(self.sc_fpix.shape[0])
    
    if self.recursive:
      X1 = self.sc_fpix[inds] / self.sc_fraw[inds].reshape(-1, 1)
    else:
      X1 = self.sc_fpix[inds] / self.sc_flux[inds].reshape(-1, 1)
    
    model = np.zeros(len(inds))
    for ii, w in zip(multichoose(range(self.sc_fpix.shape[1]), order), weights):
      model += np.product([X1[:,i] for i in ii], axis = 0) * w
    
    return model
  
class nPLDBase(object):
  '''
  
  '''
        
  def get_X(self):
    '''
    Computes the design matrix at the current *PLD* order and stores it in
    :py:obj:`self.X`. This is a list of matrices, one for each *PLD* order.
    The columns in each matrix are the *PLD* vectors for the target at the
    corresponding order, computed as the product of the fractional pixel
    flux of all sets of :py:obj:`n` pixels, where :py:obj:`n` is the *PLD*
    order. At the end of each matrix, columns corresponding to the neighbor
    star *PLD* signals are appended. Note that for both speed and memory
    reasons, cross terms are **not** computed for the neighboring stars.
      
    '''
      
    if not self.is_parent:
      log.info("Computing the design matrix...")
    if self.recursive:
      X1 = self.fpix / self.flux.reshape(-1, 1)
    else:
      X1 = self.fpix / self.fraw.reshape(-1, 1)
    for n in range(self.pld_order): 
      if (self._X[n] is None) and ((n == self.lam_idx) or (self.lam[0][n] is not None)):
        self._X[n] = np.product(list(multichoose(X1.T, n + 1)), axis = 1).T
        self._X[n] = np.hstack([self._X[n], self.X1N ** (n + 1)])
  
  def get_sc_model(self, order, weights, inds = None):
    '''
    Computes the short cadence model, including the contribution of the
    neighboring stars. See :py:meth:`rPLD.get_sc_model` for more details.
    
    '''

    if self.recursive:
      X1 = self.sc_fpix[inds] / self.sc_fraw[inds].reshape(-1, 1)
    else:
      X1 = self.sc_fpix[inds] / self.sc_flux[inds].reshape(-1, 1)
    
    model = np.zeros(len(inds))
    for ii, w, n in zip(multichoose(range(self.sc_fpix.shape[1]), order), weights, range(len(weights))):
      model += np.product([X1[:,i] for i in ii], axis = 0) * w
    
    # Add the neighbors' contribution
    model += np.dot(self.sc_X1N[inds] ** (order), weights[n + 1:])
    
    return model