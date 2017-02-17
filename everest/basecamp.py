#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`basecamp.py` - The Everest base class
----------------------------------------------

The :py:obj:`everest` engine. All :py:obj:`everest` models
inherit from :py:class:`Basecamp`.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from . import missions
from .utils import InitLog, Formatter, AP_SATURATED_PIXEL, AP_COLLAPSED_PIXEL
from .math import Chunks, Scatter, SavGol, Interpolate
from .gp import GetCovariance, GetKernelParams
from scipy.linalg import block_diag
import os, sys
import numpy as np
import george
import matplotlib.pyplot as pl
import matplotlib.image as mpimg
from matplotlib.ticker import MaxNLocator, FuncFormatter
from scipy.ndimage import zoom
from itertools import combinations_with_replacement as multichoose
import traceback
import logging
log = logging.getLogger(__name__)

__all__ = ['Basecamp']

class Basecamp(object):
  '''
  
  '''
  
  @property
  def _mission(self):
    '''
    
    '''
    
    return getattr(missions, self.mission)
  
  @_mission.setter
  def _mission(self, value):
    '''
    
    '''
    
    raise NotImplementedError("Can't set this property.") 
    
  @property
  def dir(self):
    '''
    Returns the directory where the raw data and output for the target is stored.
    
    '''
    
    return self._mission.TargetDirectory(self.ID, self.season)
      
  @dir.setter
  def dir(self, value):
    '''
    
    '''
    
    raise NotImplementedError("Can't set this property.") 
  
  @property
  def logfile(self):
    '''
    Returns the full path to the log file for the current run.
    
    '''

    return os.path.join(self.dir, '%s.log' % self.name)

  @logfile.setter
  def logfile(self, value):
    '''
    
    '''
    
    raise NotImplementedError("Can't set this property.") 
  
  @property
  def season(self):
    '''
    Returns the current observing season. For *K2*, this is the observing campaign,
    while for *Kepler*, it is the current quarter. 
    
    '''
    
    try:
      self._season
    except:
      self._season = self._mission.Season(self.ID)
      if hasattr(self._season, '__len__'):
        raise ValueError("Multiple seasons available for this target. Please specify the desired one with the `season` kwarg.")
    return self._season

  @season.setter
  def season(self, value):
    '''
    
    '''
    
    raise NotImplementedError("Can't set this property.") 
          
  @property
  def flux(self):
    '''
    The corrected/de-trended flux. This is computed by subtracting the linear
    model from the raw SAP flux.
    
    '''
    
    if self.nsub == 1:
      return self.fraw - self.model
    else:
      return [f - m for f, m in zip(self.fraw, self.model)]
  
  @flux.setter
  def flux(self, value):
    '''
    
    '''
    
    raise NotImplementedError("Can't set this property.") 

  @property
  def fcor(self):
    '''
    The CBV-corrected de-trended flux.
    
    '''
    
    if self.XCBV is None:
      return None
    else:
      if self.nsub == 1:
        return self.flux - self._mission.FitCBVs(self)
      else:
        return [f - m for f, m in zip(self.flux, self._mission.FitCBVs(self))]
  
  @fcor.setter
  def fcor(self, value):
    '''
    
    '''
    
    raise NotImplementedError("Can't set this property.") 
  
  @property
  def norm(self):
    '''
    The PLD normalization. Typically, this is just the simple aperture
    photometry flux (i.e., the sum of all the pixels in the aperture).
    
    '''
    
    return self._norm
    
  @norm.setter
  def norm(self, value):
    '''
    
    '''
    
    raise NotImplementedError("Can't set this property.") 

  @property
  def _breakpoints(self):
    '''
    A flat version of the breakpoints list.
    
    '''
    
    if self.nsub > 1:
      return [item for sublist in self.breakpoints for item in sublist]
    else:
      return list(self.breakpoints)
  
  @_breakpoints.setter
  def _breakpoints(self, value):
    '''
    
    '''
    
    raise NotImplementedError("Can't set this property.")

  @property
  def cdpps(self):
    '''
    The string version of the current value of the CDPP in *ppm*. This displays the CDPP for
    each segment of the light curve individually (if breakpoints are present).
    
    '''
    
    return " / ".join(["%.2f ppm" % c for c in self.cdpp_arr]) + (" (%.2f ppm)" % self.cdpp)
  
  @cdpps.setter
  def cdpps(self, value):
    '''
    
    '''
    
    raise NotImplementedError("Can't set this property.")
  
  @property
  def mask(self):
    '''
    The array of indices to be masked. This is the union of the sets of outliers, bad (flagged)
    cadences, transit cadences, and :py:obj:`NaN` cadences.
    
    '''

    if self.nsub > 1:
      return [np.array(list(set(np.concatenate([self.outmask[k], self.badmask[k], 
              self.transitmask[k], self.nanmask[k]]))), dtype = int) for k in range(self.nsub)]
    else:
      return np.array(list(set(np.concatenate([self.outmask, self.badmask, self.transitmask, self.nanmask]))), dtype = int)
  
  @mask.setter
  def mask(self, value):
    '''
    
    '''
    
    raise NotImplementedError("Can't set this property.")

  @property
  def weights(self):
    '''
    The PLD weights vector. The model may be computed by dotting the design matrix 
    :py:attr:`X` with this vector. Note that these are computed just for plotting
    purpoeses -- the actual weights are never explicitly computed during the de-trending,
    since it can be rather slow.
    
    '''
    
    if self._weights is None:
      self.get_weights()
    return self._weights
  
  @weights.setter
  def weights(self, value):
    '''
    
    '''
    
    raise NotImplementedError("Can't set this property.") 
  
  def get_norm(self):
    '''
    Computes the PLD normalization. In the base class, this is just
    the sum of all the pixel fluxes.
    
    '''
    
    self._norm = self.fraw
  
  def X(self, i, j = slice(None, None, None), k = slice(None, None, None)):
    '''
    Computes the design matrix at the given *PLD* order and the given indices. 
    The columns are the *PLD* vectors for the target at the
    corresponding order, computed as the product of the fractional pixel
    flux of all sets of :py:obj:`n` pixels, where :py:obj:`n` is the *PLD*
    order.
    
    :param int i: The PLD order
    :param j: The index or indices corresponding to the cadences to be returned. \
              Default :py:obj:`slice(None, None, None)`, i.e., all indices.
    :param k: The index of the sub-season (certain `K2` campaigns only). \
              Default :py:obj:`slice(None, None, None)`, i.e., no sub-season.
    
    '''

    X1 = self.fpix[k][j] / self.norm[k][j].reshape(-1, 1)
    X = np.product(list(multichoose(X1.T, i + 1)), axis = 1).T
    if self.X1N is not None:
      return np.hstack([X, self.X1N[k][j] ** (i + 1)])
    else:
      return X
  
  def subseason(self, b):
    '''
    
    '''
    
    if self.nsub > 1:
      ksz = np.cumsum([len(k) for k in self.breakpoints])
      return np.argmax(b < ksz)
    else:
      return slice(None, None, None)
  
  def compute(self):
    '''
    Compute the model for the current value of lambda.

    '''

    log.info('Computing the model...')
        
    # Loop over all chunks
    model = [None for b in self._breakpoints]
    for b, brkpt in enumerate(self._breakpoints):
      
      # Get the sub-season
      k = self.subseason(b)
      
      # Masks for current chunk
      m = self.get_masked_chunk(b)
      c = self.get_chunk(b)
      
      # This block of the masked covariance matrix
      mK = GetCovariance(self.kernel_params, self.time[k][m], self.fraw_err[k][m])
      
      # Get median
      med = np.nanmedian(self.fraw[k][m])
      
      # Normalize the flux
      f = self.fraw[k][m] - med
      
      # The X^2 matrices
      A = np.zeros((len(m), len(m)))
      B = np.zeros((len(c), len(m)))
      
      # Loop over all orders
      for n in range(self.pld_order):

        # Only compute up to the current PLD order
        if (self.lam_idx >= n) and (self.lam[b][n] is not None):
          XM = self.X(n,m,k)
          XC = self.X(n,c,k)
          A += self.lam[b][n] * np.dot(XM, XM.T)
          B += self.lam[b][n] * np.dot(XC, XM.T)
          del XM, XC
      
      W = np.linalg.solve(mK + A, f)
      model[b] = np.dot(B, W)
      del A, B, W

    # The standard case first
    if self.nsub == 1:

      # Join the chunks after applying the correct offset
      if len(model) > 1:

        # First chunk
        self.model = model[0][:-self.bpad]
  
        # Center chunks
        for m in model[1:-1]:
          # Join the chunks at the first non-outlier cadence
          if len(self.model) - 1 in self.mask:
            i = np.argmax(-np.diff(self.mask[::-1]) > 1) + 2
          else:
            i = 1
          if i <= self.bpad:
            offset = self.model[-i] - m[self.bpad - i]
          else:
            # All data points in the padded region are bad;
            # We have no choice, so let's align the chunks anyways
            offset = self.model[-1] - m[self.bpad - 1]
          self.model = np.concatenate([self.model, m[self.bpad:-self.bpad] + offset])
  
        # Last chunk
        if len(self.model) - 1 in self.mask:
          i = np.argmax(-np.diff(self.mask[::-1]) > 1) + 2
        else:
          i = 1
        if i <= self.bpad:
          offset = self.model[-i] - m[self.bpad - i]
        else:
          # All data points in the padded region are bad;
          # We have no choice, so let's align the chunks anyways
          offset = self.model[-1] - m[self.bpad - 1]
        self.model = np.concatenate([self.model, model[-1][self.bpad:] + offset])      
  
      else:
        
        # Cake!
        self.model = model[0]
  
      # Subtract the global median
      self.model -= np.nanmedian(self.model)
    
    # Now the case with multiple sub-seasons
    else:

      # Loop over all sub-seasons
      self.model = [None for i in range(self.nsub)]
      n = 0
      for k in range(self.nsub):
        
        # Loop over all the segments and assign
        # them to the correct sub-season
        nseg = len(self.breakpoints[k])
        for j, m in enumerate(model[n:n+nseg]):
          
          # Only chunk?
          if nseg == 1:
            self.model[k] = m
            
          # First chunk?
          elif j == 0:
            self.model[k] = m[:-self.bpad]
                  
          # Last chunk?
          elif j == nseg - 1:
            if len(self.model[k]) - 1 in self.mask[k]:
              i = np.argmax(-np.diff(self.mask[k][::-1]) > 1) + 2
            else:
              i = 1
            if i <= self.bpad:
              offset = self.model[k][-i] - m[self.bpad - i]
            else:
              # All data points in the padded region are bad;
              # We have no choice, so let's align the chunks anyways
              offset = self.model[k][-1] - m[self.bpad - 1]
            self.model[k] = np.concatenate([self.model[k], m[self.bpad:] + offset])
          
          # Center chunks?
          else:
            # Join the chunks at the first non-outlier cadence
            if len(self.model[k]) - 1 in self.mask[k]:
              i = np.argmax(-np.diff(self.mask[k][::-1]) > 1) + 2
            else:
              i = 1
            if i <= self.bpad:
              offset = self.model[k][-i] - m[self.bpad - i]
            else:
              # All data points in the padded region are bad;
              # We have no choice, so let's align the chunks anyways
              offset = self.model[k][-1] - m[self.bpad - 1]
            self.model[k] = np.concatenate([self.model[k], m[self.bpad:-self.bpad] + offset])
        
        n += nseg
          
        # Subtract the global median
        self.model[k] -= np.nanmedian(self.model[k])
        
    # Get the CDPP and reset the weights
    self.cdpp_arr = self.get_cdpp_arr()
    self.cdpp = self.get_cdpp()
    self._weights = None

  def compute_joint(self):
    '''
    Compute the model in a single step, allowing for covariance
    between cadences in different chunks. This should in principle
    help remove kinks at the breakpoints, but is more expensive to
    compute. 
    
    .. warning:: Everest does not cross-validate with this sort of \
                 model (it's too expensive), so care is needed when using \
                 light curves de-trended with this method, as they may not \
                 be optimized against over-fitting/under-fitting. I suspect \
                 it doesn't matter much, though.
    
    '''

    # This could be implemented eventually...
    if self.nsub > 1:
      raise NotImplementedError("Joint model evaluation not implemented when sub-seasons are present.")
  
    log.info('Computing the model...')
    A = [None for b in self._breakpoints]
    B = [None for b in self._breakpoints]
    
    # Loop over all chunks
    for b, brkpt in enumerate(self._breakpoints):

      # Masks for current chunk
      m = self.get_masked_chunk(b, pad = False)
      c = self.get_chunk(b, pad = False)
      
      # The X^2 matrices
      A[b] = np.zeros((len(m), len(m)))
      B[b] = np.zeros((len(c), len(m)))
      
      # Loop over all orders
      for n in range(self.pld_order):

        # Only compute up to the current PLD order
        if (self.lam_idx >= n) and (self.lam[b][n] is not None):
          XM = self.X(n,m)
          XC = self.X(n,c)
          A[b] += self.lam[b][n] * np.dot(XM, XM.T)
          B[b] += self.lam[b][n] * np.dot(XC, XM.T)
          del XM, XC
    
    # Merge chunks. BIGA and BIGB are sparse, but unfortunately
    # scipy.sparse doesn't handle sparse matrix inversion all that
    # well when the *result* is not itself sparse. So we're sticking
    # with regular np.linalg.
    BIGA = block_diag(*A)
    del A
    BIGB = block_diag(*B)
    del B
    
    # Compute the model
    mK = GetCovariance(self.kernel_params, self.apply_mask(self.time), self.apply_mask(self.fraw_err))
    f = self.apply_mask(self.fraw)
    f -= np.nanmedian(f)
    W = np.linalg.solve(mK + BIGA, f)
    self.model = np.dot(BIGB, W)

    # Subtract the global median
    self.model -= np.nanmedian(self.model)
    
    # Get the CDPP and reset the weights
    self.cdpp_arr = self.get_cdpp_arr()
    self.cdpp = self.get_cdpp()
    self._weights = None

  def apply_mask(self, x = None, k = slice(None, None, None)):
    '''
    Returns the outlier mask, an array of indices corresponding to the non-outliers.
    
    :param numpy.ndarray x: If specified, returns the masked version of :py:obj:`x` instead. Default :py:obj:`None`
    
    '''
    
    if x is None:
      return np.delete(np.arange(len(self.time[k])), self.mask[k])
    else:
      return np.delete(x, self.mask[k], axis = 0)

  def get_chunk(self, b, x = None, pad = True):
    '''
    Returns the indices corresponding to a given light curve chunk.
    
    :param int b: The index of the chunk to return
    :param numpy.ndarray x: If specified, applies the mask to array :py:obj:`x`. Default :py:obj:`None`

    '''
    
    # Get the subseason
    k = self.subseason(b)
    M = np.arange(len(self.time[k]))
    if b > 0 and self.subseason(b - 1) == k:
      res = M[(M > self._breakpoints[b - 1] - int(pad) * self.bpad) & (M <= self._breakpoints[b] + int(pad) * self.bpad)]
    else:
      res = M[M <= self._breakpoints[b] + int(pad) * self.bpad]
    if x is None:
      return res
    else:
      return x[res]
    
  def get_masked_chunk(self, b, x = None, pad = True):
    '''
    Same as :py:meth:`get_chunk`, but first removes the outlier indices.
    :param int b: The index of the chunk to return
    :param numpy.ndarray x: If specified, applies the mask to array :py:obj:`x`. Default :py:obj:`None`
    
    '''
    
    # Get the subseason
    k = self.subseason(b)
    M = self.apply_mask(np.arange(len(self.time[k])), k = k)
    if b > 0 and self.subseason(b - 1) == k:
      res = M[(M > self._breakpoints[b - 1] - int(pad) * self.bpad) & (M <= self._breakpoints[b] + int(pad) * self.bpad)]
    else:
      res = M[M <= self._breakpoints[b] + int(pad) * self.bpad]
    if x is None:
      return res
    else:
      return x[res]
    
  def get_weights(self):
    '''
    Computes the PLD weights vector :py:obj:`w`.
    Not currently used in the code.
    
    '''
    
    # This could be implemented eventually...
    if self.nsub > 1:
      raise NotImplementedError("This function is not implemented for light curves with sub-seasons.")
    
    log.info("Computing PLD weights...")
    
    # Loop over all chunks
    weights = [None for i in range(len(self.breakpoints))]
    for b, brkpt in enumerate(self.breakpoints):

      # Masks for current chunk
      m = self.get_masked_chunk(b)
      c = self.get_chunk(b)
      
      # This block of the masked covariance matrix
      _mK = GetCovariance(self.kernel_params, self.time[m], self.fraw_err[m])
      
      # This chunk of the normalized flux
      f = self.fraw[m] - np.nanmedian(self.fraw)  
      
      # Loop over all orders
      _A = [None for i in range(self.pld_order)]
      for n in range(self.pld_order):
        if self.lam_idx >= n:
          X = self.X(n,m)
          _A[n] = np.dot(X, X.T)
          del X
          
      # Compute the weights
      A = np.sum([l * a for l, a in zip(self.lam[b], _A) if l is not None], axis = 0)
      W = np.linalg.solve(_mK + A, f)
      weights[b] = [l * np.dot(self.X(n,m).T, W) for n, l in enumerate(self.lam[b]) if l is not None]
    
    self._weights = weights
  
  def get_cdpp_arr(self, flux = None):
    '''
    Returns the CDPP value in *ppm* for each of the chunks in the light curve.
    
    '''
    
    if flux is None:
      flux = self.flux
    return np.array([self._mission.CDPP(flux[self.subseason(b)][self.get_masked_chunk(b)], cadence = self.cadence) for b in range(self.nseg)])
  
  def get_cdpp(self, flux = None):
    '''
    Returns the scalar CDPP for the light curve.
    
    '''
    
    if flux is None:
      flux = self.flux
      if self.nsub > 1:
        for k in range(self.nsub):
          flux[k] = self.apply_mask(flux[k], k = k)
        flux = np.concatenate(flux)
      else:
        flux = self.apply_mask(flux)
    return self._mission.CDPP(flux, cadence = self.cadence)
  
  def plot_aperture(self, axes, labelsize = 8):
    '''
    Plots the aperture and the pixel images at the beginning, middle, and end of 
    the time series. Also plots a high resolution image of the target, if available.
    
    '''
    
    log.info('Plotting the aperture...')
        
    # Get colormap
    try:
      plasma = pl.get_cmap('plasma')
    except ValueError:
      plasma = pl.get_cmap('Greys')
    plasma.set_bad(alpha = 0)

    # Get aperture contour
    def PadWithZeros(vector, pad_width, iaxis, kwargs):
      vector[:pad_width[0]] = 0
      vector[-pad_width[1]:] = 0
      return vector
    
    # Plot first, mid, and last TPF image
    title = ['start', 'mid', 'end']
    for i, image in enumerate(self.pixel_images):
      
      # Sub-seasons?
      if self.nsub > 1:
        if i == 0:
          # First sub-season
          k = 0
        elif i == 2:
          # Last sub-season
          k = self.nsub - 1
        else:
          # Middle sub-season, round down
          k = (self.nsub - 1) // 2
      else:
        k = slice(None, None, None)
      
      # Grab the image
      ny, nx = image.shape
      contour = np.zeros((ny,nx))
      contour[np.where(self.aperture[k])] = 1
      contour = np.lib.pad(contour, 1, PadWithZeros)
      highres = zoom(contour, 100, order = 0, mode='nearest') 
      extent = np.array([-1, nx, -1, ny])
      
      # Plot it
      ax = axes[i]
      if np.all(np.isnan(self.flux[k])):
        alpha = 0.3
      else:
        alpha = 1
      ax.imshow(image, aspect = 'auto', interpolation = 'nearest', cmap = plasma, alpha = alpha)
      ax.contour(highres, levels=[0.5], extent=extent, origin='lower', colors='r', linewidths=1)
      
      # Check for saturated columns
      for x in range(self.aperture[k].shape[0]):
        for y in range(self.aperture[k].shape[1]):
          if self.aperture[k][x][y] == AP_SATURATED_PIXEL:
            ax.fill([y - 0.5, y + 0.5, y + 0.5, y - 0.5], 
                    [x - 0.5, x - 0.5, x + 0.5, x + 0.5], fill = False, hatch='xxxxx', color = 'r', lw = 0)      
      
      ax.axis('off')       
      ax.set_xlim(-0.7, nx - 0.3)
      ax.set_ylim(-0.7, ny - 0.3)
      ax.annotate(title[i], xy = (0.5, 0.975), xycoords = 'axes fraction',
                  ha = 'center', va = 'top', size = labelsize, color = 'w')
      if i == 1:
        for source in self.nearby:
          ax.annotate('%.1f' % source['mag'], 
                      xy = (source['x'] - source['x0'], source['y'] - source['y0']), 
                      ha = 'center', va = 'center', size = labelsize - 2, color = 'w', fontweight = 'bold')    
      
    # Plot hi res image
    if self.hires is not None:
      ax = axes[-1]
      ax.imshow(self.hires, aspect = 'auto', extent = (-0.5, nx - 0.5, -0.5, ny - 0.5), interpolation = 'bicubic', cmap = plasma)
      ax.contour(highres, levels=[0.5], extent=extent, origin='lower', colors='r', linewidths=1)
      ax.axis('off')
      ax.set_xlim(-0.7, nx - 0.3)
      ax.set_ylim(-0.7, ny - 0.3)
      ax.annotate('hires', xy = (0.5, 0.975), xycoords = 'axes fraction',
                  ha = 'center', va = 'top', size = labelsize, color = 'w')