#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`basecamp.py` - The Everest base class
----------------------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from . import missions
from .utils import InitLog, Formatter, AP_SATURATED_PIXEL, AP_COLLAPSED_PIXEL
from .math import Chunks, RMS, CDPP6, SavGol, Interpolate
from .gp import GetCovariance, GetKernelParams
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
    
    return self.fraw - self.model
  
  @flux.setter
  def flux(self, value):
    '''
    
    '''
    
    raise NotImplementedError("Can't set this property.") 
  
  @property
  def cdpps(self):
    '''
    The string version of the current value of the CDPP in *ppm*. This displays the CDPP for
    each segment of the light curve individually (if breakpoints are present).
    
    '''
    
    return " / ".join(["%.2f ppm" % c for c in self.cdpp6_arr]) + (" (%.2f ppm)" % self.cdpp6)
  
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
    
    return np.array(list(set(np.concatenate([self.outmask, self.badmask, self.transitmask, self.nanmask]))), dtype = int)
  
  @mask.setter
  def mask(self, value):
    '''
    
    '''
    
    raise NotImplementedError("Can't set this property.")
    
  @property
  def X(self):
    '''
    The current *PLD* design matrix.
    
    '''
    
    return self._DesignMatrix
  
  @X.setter
  def X(self, value):
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
  
  def compute(self):
    '''
    Compute the model for the current value of lambda.

    '''

    log.info('Computing the model...')

    # Loop over all chunks
    model = [None for b in self.breakpoints]
    for b, brkpt in enumerate(self.breakpoints):
    
      # Masks for current chunk
      m = self.get_masked_chunk(b)
      c = self.get_chunk(b)
      
      # This block of the masked covariance matrix
      mK = GetCovariance(self.kernel_params, self.time[m], self.fraw_err[m])
      
      # The normalized flux
      f = self.fraw[m] - np.nanmedian(self.fraw[m]) # ***
      
      # The X^2 matrices
      A = np.zeros((len(m), len(m)))
      B = np.zeros((len(c), len(m)))
      
      # Loop over all orders
      for n in range(self.pld_order):

        # Only compute up to the current PLD order
        if (self.lam_idx >= n) and (self.lam[b][n] is not None):
          XM = self.X(n,m)
          XC = self.X(n,c)
          A += self.lam[b][n] * np.dot(XM, XM.T)
          B += self.lam[b][n] * np.dot(XC, XM.T)
          del XM, XC
      
      W = np.linalg.solve(mK + A, f)
      model[b] = np.dot(B, W)
      del A, B, W

    # Join the chunks after applying the correct offset
    if len(model) > 1:

      # First chunk
      self.model = model[0][:-self.bpad]
  
      # Center chunks
      for m in model[1:-1]:
        offset = self.model[-1] - m[self.bpad - 1]
        self.model = np.concatenate([self.model, m[self.bpad:-self.bpad] + offset])
  
      # Last chunk
      offset = self.model[-1] - model[-1][self.bpad - 1]
      self.model = np.concatenate([self.model, model[-1][self.bpad:] + offset])      
  
    else:

      self.model = model[0]
  
    # Subtract the global median
    self.model -= np.nanmedian(self.model)
    
    # Get the CDPP and reset the weights
    self.cdpp6_arr = self.get_cdpp_arr()
    self.cdpp6 = self.get_cdpp()
    self._weights = None

  def apply_mask(self, x = None):
    '''
    Returns the outlier mask, an array of indices corresponding to the non-outliers.
    
    :param numpy.ndarray x: If specified, returns the masked version of :py:obj:`x` instead. Default :py:obj:`None`
    
    '''
    
    if x is None:
      return np.delete(np.arange(len(self.time)), self.mask)
    else:
      return np.delete(x, self.mask, axis = 0)

  def get_chunk(self, b, x = None):
    '''
    Returns the indices corresponding to a given light curve chunk.
    
    :param int b: The index of the chunk to return
    :param numpy.ndarray x: If specified, applies the mask to array :py:obj:`x`. Default :py:obj:`None`

    '''

    M = np.arange(len(self.time))
    if b > 0:
      res = M[(M > self.breakpoints[b - 1] - self.bpad) & (M <= self.breakpoints[b] + self.bpad)]
    else:
      res = M[M <= self.breakpoints[b] + self.bpad]
    if x is None:
      return res
    else:
      return x[res]
    
  def get_masked_chunk(self, b, x = None):
    '''
    Same as :py:meth:`get_chunk`, but first removes the outlier indices.
    :param int b: The index of the chunk to return
    :param numpy.ndarray x: If specified, applies the mask to array :py:obj:`x`. Default :py:obj:`None`
    
    '''
    
    M = self.apply_mask(np.arange(len(self.time)))
    if b > 0:
      res = M[(M > self.breakpoints[b - 1] - self.bpad) & (M <= self.breakpoints[b] + self.bpad)]
    else:
      res = M[M <= self.breakpoints[b] + self.bpad]
    if x is None:
      return res
    else:
      return x[res]
    
  def get_weights(self):
    '''
    Computes the PLD weights vector :py:obj:`w` (for plotting and/or computing the short
    cadence model).
    
    '''
    
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
    Returns the 6-hr CDPP value in *ppm* for each of the chunks in the light curve.
    
    '''
    
    if flux is None:
      flux = self.flux
    return np.array([CDPP6(flux[self.get_masked_chunk(b)], cadence = self.cadence) for b, _ in enumerate(self.breakpoints)])
  
  def get_cdpp(self, flux = None):
    '''
    Returns the scalar 6-hr CDPP for the light curve.
    
    '''
    
    if flux is None:
      flux = self.flux
    return CDPP6(self.apply_mask(flux), cadence = self.cadence)
  
  def plot_weights(self, ax, cax):
    '''
    Plots the *PLD* weights on the CCD for each of the *PLD* orders.
    
    .. note:: Only works for light curves with zero or one breakpoints.
    
    '''
    
    # Check number of segments
    if len(self.breakpoints) > 2:
      return
    
    # Loop over all PLD orders and over all chunks
    npix = len(self.fpix[1])
    ap = self.aperture.flatten()
    ncol = 1 + 2 * (len(self.weights[0]) - 1)
    raw_weights = np.zeros((len(self.breakpoints), ncol, self.aperture.shape[0], self.aperture.shape[1]), dtype = float)
    scaled_weights = np.zeros((len(self.breakpoints), ncol, self.aperture.shape[0], self.aperture.shape[1]), dtype = float)
    
    # Loop over orders
    for o in range(len(self.weights[0])):
      if o == 0:
        oi = 0
      else:
        oi = 1 + 2 * (o - 1)
        
      # Loop over chunks
      for b in range(len(self.weights)):
      
        c = self.get_chunk(b)
        rw_ii = np.zeros(npix); rw_ij = np.zeros(npix)
        sw_ii = np.zeros(npix); sw_ij = np.zeros(npix)
        X = np.nanmedian(self.X(o,c), axis = 0)
      
        # Compute all sets of pixels at this PLD order, then
        # loop over them and assign the weights to the correct pixels
        sets = np.array(list(multichoose(np.arange(npix).T, o + 1)))
        for i, s in enumerate(sets):
          if (o == 0) or (s[0] == s[1]):
            # Not the cross-terms
            j = s[0]
            rw_ii[j] += self.weights[b][o][i]
            sw_ii[j] += X[i] * self.weights[b][o][i]
          else:
            # Cross-terms
            for j in s:
              rw_ij[j] += self.weights[b][o][i]
              sw_ij[j] += X[i] * self.weights[b][o][i]
          
        # Make the array 2D and plot it
        rw = np.zeros_like(ap, dtype = float)
        sw = np.zeros_like(ap, dtype = float)
        n = 0
        for i, a in enumerate(ap):
          if (a & 1):
            rw[i] = rw_ii[n]
            sw[i] = sw_ii[n]
            n += 1
        raw_weights[b][oi] = rw.reshape(*self.aperture.shape)
        scaled_weights[b][oi] = sw.reshape(*self.aperture.shape)

        if o > 0:
          # Make the array 2D and plot it
          rw = np.zeros_like(ap, dtype = float)
          sw = np.zeros_like(ap, dtype = float)
          n = 0
          for i, a in enumerate(ap):
            if (a & 1):
              rw[i] = rw_ij[n]
              sw[i] = sw_ij[n]
              n += 1
          raw_weights[b][oi + 1] = rw.reshape(*self.aperture.shape)
          scaled_weights[b][oi + 1] = sw.reshape(*self.aperture.shape)

    # Plot the images
    log.info('Plotting the PLD weights...')
    rdbu = pl.get_cmap('RdBu_r')
    rdbu.set_bad('k')
    for b in range(len(self.weights)):
      rmax = max([-raw_weights[b][o].min() for o in range(ncol)] +
                 [raw_weights[b][o].max() for o in range(ncol)])
      smax = max([-scaled_weights[b][o].min() for o in range(ncol)] +
                 [scaled_weights[b][o].max() for o in range(ncol)])
      for o in range(ncol):
        imr = ax[2 * b, o].imshow(raw_weights[b][o], aspect = 'auto', interpolation = 'nearest', cmap = rdbu, origin = 'lower', vmin = -rmax, vmax = rmax)
        ims = ax[2 * b + 1, o].imshow(scaled_weights[b][o], aspect = 'auto', interpolation = 'nearest', cmap = rdbu, origin = 'lower', vmin=-smax, vmax=smax)
    
      # Colorbars
      def fmt(x, pos):
        a, b = '{:.0e}'.format(x).split('e')
        b = int(b)
        if float(a) > 0:
          a = r'+' + a
        elif float(a) == 0:
          return ''
        return r'${} \times 10^{{{}}}$'.format(a, b) 
      cbr = pl.colorbar(imr, cax = cax[2 * b], format = FuncFormatter(fmt))
      cbr.ax.tick_params(labelsize = 8) 
      cbs = pl.colorbar(ims, cax = cax[2 * b + 1], format = FuncFormatter(fmt))
      cbs.ax.tick_params(labelsize = 8) 
  
    # Plot aperture contours
    def PadWithZeros(vector, pad_width, iaxis, kwargs):
      vector[:pad_width[0]] = 0
      vector[-pad_width[1]:] = 0
      return vector
    ny, nx = self.aperture.shape
    contour = np.zeros((ny,nx))
    contour[np.where(self.aperture)] = 1
    contour = np.lib.pad(contour, 1, PadWithZeros)
    highres = zoom(contour, 100, order = 0, mode='nearest') 
    extent = np.array([-1, nx, -1, ny])
    for axis in ax.flatten():
      axis.contour(highres, levels=[0.5], extent=extent, origin='lower', colors='r', linewidths=1)
      
      # Check for saturated columns
      for x in range(self.aperture.shape[0]):
        for y in range(self.aperture.shape[1]):
          if self.aperture[x][y] == AP_SATURATED_PIXEL:
            axis.fill([y - 0.5, y + 0.5, y + 0.5, y - 0.5], 
                      [x - 0.5, x - 0.5, x + 0.5, x + 0.5], fill = False, hatch='xxxxx', color = 'r', lw = 0)
      
      axis.set_xlim(-0.5, nx - 0.5)
      axis.set_ylim(-0.5, ny - 0.5)
      axis.set_xticks([]) 
      axis.set_yticks([])
  
    # Labels
    titles = [r'$1^{\mathrm{st}}$', 
              r'$2^{\mathrm{nd}}\ (i = j)$',
              r'$2^{\mathrm{nd}}\ (i \neq j)$',
              r'$3^{\mathrm{rd}}\ (i = j)$',
              r'$3^{\mathrm{rd}}\ (i \neq j)$'] + ['' for i in range(10)]
    for i, axis in enumerate(ax[0]):
      axis.set_title(titles[i], fontsize = 12)
    for j in range(len(self.weights)):
      ax[2 * j, 0].text(-0.55, -0.15, r'$%d$' % (j + 1), fontsize = 16, transform = ax[2 * j, 0].transAxes)
      ax[2 * j, 0].set_ylabel(r'$w_{ij}$', fontsize = 18)
      ax[2 * j + 1, 0].set_ylabel(r'$\bar{X}_{ij} \cdot w_{ij}$', fontsize = 18)

  def plot_aperture(self, axes, labelsize = 8):
    '''
    Plots the aperture and the pixel images at the beginning, middle, and end of 
    the time series. Also plots a high resolution image of the target, if available.
    
    '''
    
    log.info('Plotting the aperture...')
        
    # Get colormap
    plasma = pl.get_cmap('plasma')
    plasma.set_bad(alpha = 0)
    
    # Get aperture contour
    def PadWithZeros(vector, pad_width, iaxis, kwargs):
      vector[:pad_width[0]] = 0
      vector[-pad_width[1]:] = 0
      return vector
    ny, nx = self.pixel_images[0].shape
    contour = np.zeros((ny,nx))
    contour[np.where(self.aperture)] = 1
    contour = np.lib.pad(contour, 1, PadWithZeros)
    highres = zoom(contour, 100, order = 0, mode='nearest') 
    extent = np.array([-1, nx, -1, ny])
    
    # Plot first, mid, and last TPF image
    title = ['start', 'mid', 'end']
    for i, image in enumerate(self.pixel_images):
      ax = axes[i]
      ax.imshow(image, aspect = 'auto', interpolation = 'nearest', cmap = plasma)
      ax.contour(highres, levels=[0.5], extent=extent, origin='lower', colors='r', linewidths=1)
      
      # Check for saturated columns
      for x in range(self.aperture.shape[0]):
        for y in range(self.aperture.shape[1]):
          if self.aperture[x][y] == AP_SATURATED_PIXEL:
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