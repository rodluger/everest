#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`model.py` - De-trending
--------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .config import EVEREST_DAT
from .missions import Missions
from .utils import InitLog, Formatter, AP_SATURATED_PIXEL, AP_COLLAPSED_PIXEL
from .math import Chunks, RMS, CDPP6, SavGol, Interpolate
from .data import Season, GetData, GetNeighbors, Breakpoint, GetSimpleNeighbors
from .gp import GetCovariance, GetKernelParams
from .transit import Transit
from .dvs import DVS1, DVS2
from .model import Model
import os, sys, glob
import numpy as np
import george
import matplotlib.pyplot as pl
import matplotlib.image as mpimg
from matplotlib.ticker import MaxNLocator, FuncFormatter
from matplotlib.backends.backend_pdf import PdfPages
from scipy.ndimage import zoom
from itertools import combinations_with_replacement as multichoose
import traceback
import logging
log = logging.getLogger(__name__)

__all__ = ['cPLD', 'SimpleNeighborPLD', 'NoThrustersNeighborPLD']

class cPLD(Model):
  '''
  PLD for crowded stars. Uses the PLD vectors of 50
  bright, well-de-trended neighbors to de-trend each target. Unlike
  nPLD, the target's own pixels are not used.
  
  '''
        
  def __init__(self, *args, **kwargs):
    '''
    
    '''
    
    # Initialize
    super(cPLD, self).__init__(*args, **kwargs)
    
    # Check for saved model
    if self.load_model():
      return
    
    # Get neighbors
    num_neighbors = kwargs.get('neighbors', 50)
    self.neighbors = GetNeighbors(self.ID, mission = self.mission, 
                                  model = 'PLD',
                                  neighbors = num_neighbors, 
                                  mag_lo = kwargs.get('mag_lo', 11.), 
                                  mag_hi = kwargs.get('mag_hi', 13.),
                                  cdpp_lo = kwargs.get('cdpp_lo', 10.), 
                                  cdpp_hi = kwargs.get('cdpp_hi', 30.))
    if len(self.neighbors):
      if len(self.neighbors) < num_neighbors:
        log.warn("%d neighbors requested, but only %d found." % (num_neighbors, len(self.neighbors)))
    else:
      log.error("No neighbors found! Aborting.")
      return
      
    self._XNeighbors = [None for i in range(self.pld_order)]
    for neighbor in self.neighbors:
      log.info("Loading data for neighboring target %d..." % neighbor)
      data = PLD(neighbor, mission = self.mission, is_child = True)
      X1 = data.fpix / data.fraw.reshape(-1, 1)
      # Linearly interpolate over outliers, NaNs and bad timestamps
      X1 = Interpolate(data.time, data.mask, X1)
      for n in range(self.pld_order):
        if self._XNeighbors[n] is None:
          self._XNeighbors[n] = X1 ** (n + 1)
        else:
          self._XNeighbors[n] = np.hstack([self._XNeighbors[n], X1 ** (n + 1)])
      del data
      
    # Run
    self.run()
  
  def plot_weights(self):
    '''
    
    '''
    
    pass
     
  def get_X(self):
    '''
  
    '''
      
    if not self.is_child:
      log.info("Computing the design matrix...")
    X1 = self.fpix / self.flux.reshape(-1, 1)
    for n in range(self.pld_order): 
      if (self._X[n] is None) and ((n == self.lam_idx) or (self.lam[0][n] is not None)):
        self._X[n] = np.array(self._XNeighbors[n])
        
class SimpleNeighborPLD(Model):
  '''
  Neighbor PLD, but with no dependence on previous PLD
  runs. This version does not identify outliers, such as thruster
  firing events. Compare this to the output of PLD + nPLD.
  
  '''
        
  def __init__(self, *args, **kwargs):
    '''
    
    '''
    
    # Initialize
    super(SimpleNeighborPLD, self).__init__(*args, **kwargs)
    
    # Check for saved model
    if self.load_model():
      return
    
    # Get neighbors
    num_neighbors = kwargs.get('neighbors', 10)
    self.neighbors = GetSimpleNeighbors(self.ID, mission = self.mission, 
                                        neighbors = num_neighbors, 
                                        mag_lo = kwargs.get('mag_lo', 11.), 
                                        mag_hi = kwargs.get('mag_hi', 13.))
    if len(self.neighbors):
      if len(self.neighbors) < num_neighbors:
        log.warn("%d neighbors requested, but only %d found." % (num_neighbors, len(self.neighbors)))
    else:
      log.error("No neighbors found! Aborting.")
      return
      
    self._XNeighbors = [None for i in range(self.pld_order)]
    for neighbor in self.neighbors:
      log.info("Loading data for neighboring target %d..." % neighbor)
      # Get the neighbor's flux data
      data = GetData(neighbor, self.mission, season = self.season, clobber = self.clobber_tpf, 
                     aperture_name = self.aperture_name, 
                     saturated_aperture_name = self.saturated_aperture_name, 
                     max_pixels = self.max_pixels,
                     saturation_tolerance = self.saturation_tolerance)
      mask = np.array(list(set(np.concatenate([data.badmask, data.nanmask]))), dtype = int)
      fraw = np.sum(data.fpix, axis = 1)
      X1 = data.fpix / fraw.reshape(-1, 1)
      # Linearly interpolate over outliers, NaNs and bad timestamps
      X1 = Interpolate(data.time, mask, X1)
      for n in range(self.pld_order):
        if self._XNeighbors[n] is None:
          self._XNeighbors[n] = X1 ** (n + 1)
        else:
          self._XNeighbors[n] = np.hstack([self._XNeighbors[n], X1 ** (n + 1)])
      del data
      
    # Run
    self.run()
        
  def get_X(self):
    '''
  
    '''
      
    if not self.is_child:
      log.info("Computing the design matrix...")
    X1 = self.fpix / self.flux.reshape(-1, 1)
    for n in range(self.pld_order): 
      if (self._X[n] is None) and ((n == self.lam_idx) or (self.lam[0][n] is not None)):
        self._X[n] = np.product(list(multichoose(X1.T, n + 1)), axis = 1).T
        self._X[n] = np.hstack([self._X[n], self._XNeighbors[n]])

class NoThrustersNeighborPLD(Model):
  '''
  Neighbor PLD, but with no dependence on previous PLD
  runs. This version attempts to identify and remove outliers, such as thruster
  firing events. Compare this to the output of PLD + nPLD and SimpleNeighborPLD.
  
  '''
        
  def __init__(self, *args, **kwargs):
    '''
    
    '''
    
    # Initialize
    super(NoThrustersNeighborPLD, self).__init__(*args, **kwargs)
    
    # Check for saved model
    if self.load_model():
      return
    
    # Get neighbors
    num_neighbors = kwargs.get('neighbors', 10)
    self.neighbors = GetSimpleNeighbors(self.ID, mission = self.mission,
                                        neighbors = num_neighbors, 
                                        mag_lo = kwargs.get('mag_lo', 11.), 
                                        mag_hi = kwargs.get('mag_hi', 13.))
    if len(self.neighbors):
      if len(self.neighbors) < num_neighbors:
        log.warn("%d neighbors requested, but only %d found." % (num_neighbors, len(self.neighbors)))
    else:
      log.error("No neighbors found! Aborting.")
      return
      
    self._XNeighbors = [None for i in range(self.pld_order)]
    for neighbor in self.neighbors:
      log.info("Loading data for neighboring target %d..." % neighbor)
      
      # Get the neighbor's flux data
      data = GetData(neighbor, self.mission, season = self.season, clobber = self.clobber_tpf, 
                     aperture_name = self.aperture_name, 
                     saturated_aperture_name = self.saturated_aperture_name, 
                     max_pixels = self.max_pixels,
                     saturation_tolerance = self.saturation_tolerance)
      time = np.array(data.time)
      fpix = np.array(data.fpix)
      fraw = np.sum(fpix, axis = 1)
      flux = np.array(fraw)
      fraw_err = np.sqrt(np.sum(data.fpix_err ** 2, axis = 1))
      
      # Compute linear basis vectors
      mask = np.array(list(set(np.concatenate([data.badmask, data.nanmask]))), dtype = int)
      X1 = fpix / fraw.reshape(-1, 1)
      X1 = Interpolate(time, mask, X1)
      
      # Iterative sigma clipping
      i = 0
      outmask1 = np.array([], dtype = int)
      outmask2 = np.array([-1], dtype = int)
      
      # Loop as long as we're getting a different outlier mask
      while not np.array_equal(outmask1, outmask2):
                
        # Compute the linear PLD model
        gp = george.GP(np.nanstd(flux) ** 2 * george.kernels.Matern32Kernel(30. ** 2))
        gp.compute(time, fraw_err)
        A = np.dot(X1.T, gp.solver.apply_inverse(X1))
        B = np.dot(X1.T, gp.solver.apply_inverse(fraw))
        w = np.linalg.solve(A, B)
        model = np.dot(X1, w)
        flux = fraw - model
    
        # Get the outliers
        f = SavGol(flux)
        med = np.nanmedian(f)
        MAD = 1.4826 * np.nanmedian(np.abs(f - med))
        inds = np.where((f > med + self.osigma * MAD) | (f < med - self.osigma * MAD))[0]
      
        # Loop
        outmask2 = np.array(outmask1)
        outmask1 = np.array(inds, dtype = int)
        i += 1
        if i == self.oiter:
          break
      
      log.info("%d outliers removed." % len(outmask1))
      
      # Interpolate over the outliers
      X1 = Interpolate(time, outmask1, X1)
      
      # Compute higher order terms and append
      for n in range(self.pld_order):
        if self._XNeighbors[n] is None:
          self._XNeighbors[n] = X1 ** (n + 1)
        else:
          self._XNeighbors[n] = np.hstack([self._XNeighbors[n], X1 ** (n + 1)])
      del data
      
    # Run
    self.run()
        
  def get_X(self):
    '''
  
    '''
    
    if not self.is_child:
      log.info("Computing the design matrix...")
    X1 = self.fpix / self.flux.reshape(-1, 1)
    for n in range(self.pld_order): 
      if (self._X[n] is None) and ((n == self.lam_idx) or (self.lam[0][n] is not None)):
        self._X[n] = np.product(list(multichoose(X1.T, n + 1)), axis = 1).T
        self._X[n] = np.hstack([self._X[n], self._XNeighbors[n]])
 