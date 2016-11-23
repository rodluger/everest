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