#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`basecamp.py` - The Everest base class
----------------------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .missions import Missions
from .utils import InitLog, Formatter, AP_SATURATED_PIXEL, AP_COLLAPSED_PIXEL
from .math import Chunks, RMS, CDPP6, SavGol, Interpolate
from .data import Season, Breakpoint
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
  def dir(self):
    '''
    Returns the directory where the raw data and output for the target is stored.
    
    '''
    
    return Missions[self.mission].TargetDirectory(self.ID, self.season)
      
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
      self._season = Season(self.ID, mission = self.mission)
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
  def sc_flux(self):
    '''
    The corrected/de-trended short cadence flux, computed by subtracting the
    short cadence linear model from the raw short cadence SAP flux. Only available
    if short cadence data exists for the target.
    
    '''
    
    if self.sc_model is not None:
      return self.sc_fraw - self.sc_model
    else:
      return None
  
  @sc_flux.setter
  def sc_flux(self, value):
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
  def sc_mask(self):
    '''
    Similar to :py:attr:`mask`, but for the short cadence data (if available).
    
    '''
    
    return np.array(list(set(np.concatenate([self.sc_badmask, self.sc_nanmask, self.sc_transitmask, self.sc_outmask]))), dtype = int)
  
  @sc_mask.setter
  def sc_mask(self, value):
    '''
    
    '''
    
    raise NotImplementedError("Can't set this property.")
  
  @property
  def X(self):
    '''
    The current *PLD* design matrix.
    
    '''
    
    return self._X
  
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
  
  def precompute(self):
    '''
    Pre-compute some expensive matrices used during the regularization step 
    and store them in memory.
    
    '''
    
    # Loop over all chunks
    for b, brkpt in enumerate(self.breakpoints):
    
      # Masks for current chunk
      m = self.get_masked_chunk(b)
      c = self.get_chunk(b)
      
      # Loop over all orders
      for n in range(self.pld_order):
        if self.X[n] is not None:
          self._A[b][n] = np.dot(self.X[n][m], self.X[n][m].T)
          self._B[b][n] = np.dot(self.X[n][c], self.X[n][m].T)
      
      # This block of the masked covariance matrix
      self._mK[b] = self.K[m][:,m]
      
      # Normalized flux
      self._f[b] = self.fraw[m] - np.nanmedian(self.fraw)  
  
  def compute(self, precompute = True, sc = False):
    '''
    Compute the model for the current value of lambda.
    
    :param bool precompute: Precompute the expensive :py:obj:`A` and :py:obj:`B` matrices? Default :py:obj:`True`
    :param bool sc: Compute the short cadence model? Default :py:obj:`False`
    
    '''
    
    if not sc:
      
      log.info('Computing the model...')
      if precompute:
        self.precompute()
  
      # Loop over all chunks
      model = [None for b in self.breakpoints]
      for b, brkpt in enumerate(self.breakpoints):
        A = np.sum([l * a for l, a in zip(self.lam[b], self._A[b]) if l is not None], axis = 0)
        B = np.sum([l * b for l, b in zip(self.lam[b], self._B[b]) if l is not None], axis = 0)
        W = np.linalg.solve(self._mK[b] + A, self._f[b])
        model[b] = np.dot(B, W)
  
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
      
    else:
      
      log.info('Computing the short cadence model...')
      # Get the PLD weights
      norders = len(self.weights[0])
      nchunks = len(self.weights)
      model_arr = [None for b in range(nchunks)]
    
      # Compute the dot product order by order, chunk by chunk
      # We load each segment of the design matrix individually,
      # since it can be quite large.
      for o in range(norders):
        for b in range(nchunks):
          log.info('Computing short cadence model (%d/%d)...' % (b + o * nchunks + 1, norders * nchunks))
          c = self.get_chunk(b, sc = True)
          m = self.get_sc_model(o + 1, inds = c, weights = self.weights[b][o])
          if model_arr[b] is None:
            model_arr[b] = m
          else:
            model_arr[b] += m
    
      # Join the chunks after applying the correct offset
      if len(model_arr) > 1:
        model = model_arr[0][:-self.sc_bpad]
        for m in model_arr[1:-1]:
          offset = model[-1] - m[self.sc_bpad - 1]
          model = np.concatenate([model, m[self.sc_bpad:-self.sc_bpad] + offset])
        offset = model[-1] - model_arr[-1][self.sc_bpad - 1]
        model = np.concatenate([model, model_arr[-1][self.sc_bpad:] + offset])
      else:
        model = model_arr[0]
      model -= np.nanmedian(model)
      self.sc_model = model

  def apply_mask(self, x = None, sc = False):
    '''
    Returns the outlier mask, an array of indices corresponding to the non-outliers.
    
    :param numpy.ndarray x: If specified, returns the masked version of :py:obj:`x` instead. Default :py:obj:`None`
    :param bool sc: Return the short cadence mask? Default :py:obj:`False`
    
    '''
    
    if x is None:
      if sc:
        return np.delete(np.arange(len(self.sc_time)), self.sc_mask)
      else:
        return np.delete(np.arange(len(self.time)), self.mask)
    else:
      if sc:
        return np.delete(x, self.sc_mask, axis = 0)
      else:
        return np.delete(x, self.mask, axis = 0)

  def get_chunk(self, b, x = None, sc = False):
    '''
    Returns the indices corresponding to a given light curve chunk.
    
    :param int b: The index of the chunk to return
    :param numpy.ndarray x: If specified, applies the mask to array :py:obj:`x`. Default :py:obj:`None`
    :param bool sc: Get the indices of the short cadence chunk? Default :py:obj:`False`
    
    '''
    
    if not sc:
      M = np.arange(len(self.time))
      if b > 0:
        res = M[(M > self.breakpoints[b - 1] - self.bpad) & (M <= self.breakpoints[b] + self.bpad)]
      else:
        res = M[M <= self.breakpoints[b] + self.bpad]
      if x is None:
        return res
      else:
        return x[res]
    else:
      sc_breakpoints = [np.nanargmin(np.abs(self.sc_time - self.time[b])) for b in self.breakpoints]
      M = np.arange(len(self.sc_time))
      if b > 0:
        res = M[(M > sc_breakpoints[b - 1] - self.sc_bpad) & (M <= sc_breakpoints[b] + self.sc_bpad)]
      else:
        res = M[M <= sc_breakpoints[b] + self.sc_bpad]
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
  
  def get_X(self):
    '''
    *Implemented by subclasses*
    
    '''
  
    raise NotImplementedError('This method must be implemented in a subclass.')
  
  def get_sc_model(self):
    '''
    *Implemented by subclasses*
    
    '''
    
    raise NotImplementedError('This method must be implemented in a subclass.')
  
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
      _mK = self.K[m][:,m]
      
      # This chunk of the normalized flux
      f = self.fraw[m] - np.nanmedian(self.fraw)  
      
      # Loop over all orders
      _A = [None for i in range(self.pld_order)]
      for n in range(self.pld_order):
        if self.X[n] is not None:
          _A[n] = np.dot(self.X[n][m], self.X[n][m].T)
      
      # Compute the weights
      A = np.sum([l * a for l, a in zip(self.lam[b], _A) if l is not None], axis = 0)
      W = np.linalg.solve(_mK + A, f)
      weights[b] = [l * np.dot(self.X[n][m].T, W) for n, l in enumerate(self.lam[b]) if l is not None]
    
    self._weights = weights
  
  def get_cdpp_arr(self):
    '''
    Returns the 6-hr CDPP value in *ppm* for each of the chunks in the light curve.
    
    '''
    
    return np.array([CDPP6(self.flux[self.get_masked_chunk(b)]) for b, _ in enumerate(self.breakpoints)])
  
  def get_cdpp(self):
    '''
    Returns the scalar 6-hr CDPP for the light curve.
    
    '''
    
    return CDPP6(self.apply_mask(self.flux))
  
  def plot_weights(self, ax = None, cax = None):
    '''
    Plots the *PLD* weights on the CCD for each of the *PLD* orders.
    
    '''
    
    # Set up the axes?
    if (ax is None) or (cax is None):
      fig = pl.figure(figsize = (12, 12))
      fig.subplots_adjust(top = 0.95, bottom = 0.025, left = 0.1, right = 0.92)
      fig.canvas.set_window_title('%s %d' % (Missions[self.mission].IDSTRING, self.ID))
      ax = [pl.subplot2grid((80, 130), (20 * j, 25 * i), colspan = 23, rowspan = 18) 
            for j in range(len(self.breakpoints) * 2) for i in range(1 + 2 * (self.pld_order - 1))]
      cax = [pl.subplot2grid((80, 130), (20 * j, 25 * (1 + 2 * (self.pld_order - 1))), 
             colspan = 4, rowspan = 18) for j in range(len(self.breakpoints) * 2)]
      ax = np.array(ax).reshape(2 * len(self.breakpoints), -1)
      cax = np.array(cax)
      show = True
    else:
      show = False
    
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
        X = np.nanmedian(self.X[o][c], axis = 0)
      
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
  
    if show:
      pl.show()

  def plot_aperture(self, axes = None):
    '''
    Plots the aperture and the pixel images at the beginning, middle, and end of 
    the time series. Also plots a high resolution image of the target, if available.
    
    '''
    
    log.info('Plotting the aperture...')
    
    # Set up the axes?
    if axes is None:
      fig, axes = pl.subplots(2,2, figsize = (6, 8))
      fig.subplots_adjust(top = 0.975, bottom = 0.025, left = 0.05, 
                          right = 0.95, hspace = 0.05, wspace = 0.05)
      axes = axes.flatten()
      fig.canvas.set_window_title('%s %d' % (Missions[self.mission].IDSTRING, self.ID))
      show = True
      labelsize = 12
    else:
      show = False
      labelsize = 8
    
    # Get colormap
    plasma = pl.get_cmap('plasma')
    plasma.set_bad('w')
    
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
    
    if show:
      pl.show()

  def plot(self, ax = None):
    '''
    Plots the final de-trended light curve.
    
    '''
    
    log.info('Plotting the light curve...')
    
    # Set up axes?
    if ax is None:
      fig, ax = pl.subplots(1, figsize = (13, 6))
      fig.canvas.set_window_title('EVEREST Light curve')
      show = True
    else:
      show = False

    # Plot the good data points
    ax.plot(self.apply_mask(self.time), self.apply_mask(self.flux), ls = 'none', marker = '.', color = 'k', markersize = 4, alpha = 0.5)
    
    # Plot the outliers
    bnmask = np.array(list(set(np.concatenate([self.badmask, self.nanmask]))), dtype = int)
    O1 = lambda x: x[self.outmask]
    O2 = lambda x: x[bnmask]
    M = lambda x: np.delete(x, bnmask)
    ax.plot(O1(self.time), O1(self.flux), ls = 'none', color = "#777777", marker = '.', markersize = 4, alpha = 0.5)
    ax.plot(O2(self.time), O2(self.flux), 'r.', markersize = 4, alpha = 0.25)

    # Plot the GP
    _, amp, tau = self.kernel_params
    gp = george.GP(amp ** 2 * george.kernels.Matern32Kernel(tau ** 2))
    gp.compute(self.apply_mask(self.time), self.apply_mask(self.fraw_err))
    med = np.nanmedian(self.apply_mask(self.flux))
    y, _ = gp.predict(self.apply_mask(self.flux) - med, self.time)
    y += med
    ax.plot(M(self.time), M(y), 'r-', lw = 0.5, alpha = 0.5)

    # Appearance
    ax.set_title('%s %d' % (Missions[self.mission].IDSTRING, self.ID), fontsize = 18)
    ax.set_xlabel('Time (%s)' % Missions[self.mission].TIMEUNITS, fontsize = 18)
    ax.set_ylabel('Flux', fontsize = 18)
    for brkpt in self.breakpoints[:-1]:
      ax.axvline(self.time[brkpt], color = 'r', ls = '--', alpha = 0.25)
    if len(self.cdpp6_arr) == 2:
      ax.annotate('%.2f ppm' % self.cdpp6_arr[0], xy = (0.02, 0.975), xycoords = 'axes fraction', 
                  ha = 'left', va = 'top', fontsize = 12, color = 'r')
      ax.annotate('%.2f ppm' % self.cdpp6_arr[1], xy = (0.98, 0.975), xycoords = 'axes fraction', 
                  ha = 'right', va = 'top', fontsize = 12, color = 'r')
    else:
      ax.annotate('%.2f ppm' % self.cdpp6, xy = (0.02, 0.975), xycoords = 'axes fraction', 
                  ha = 'left', va = 'top', fontsize = 12, color = 'r')
    ax.margins(0.01, 0.1)          
    
    # Get y lims that bound 99% of the flux
    flux = np.delete(self.flux, bnmask)
    N = int(0.995 * len(flux))
    hi, lo = flux[np.argsort(flux)][[N,-N]]
    fsort = flux[np.argsort(flux)]
    pad = (hi - lo) * 0.1
    ylim = (lo - pad, hi + pad)
    ax.set_ylim(ylim)   
    ax.get_yaxis().set_major_formatter(Formatter.Flux)
    
    # Indicate off-axis outliers
    for i in np.where(self.flux < ylim[0])[0]:
      if i in bnmask:
        color = "#ffcccc"
      elif i in self.outmask:
        color = "#cccccc"
      else:
        color = "#ccccff"
      ax.annotate('', xy=(self.time[i], ylim[0]), xycoords = 'data',
                  xytext = (0, 15), textcoords = 'offset points',
                  arrowprops=dict(arrowstyle = "-|>", color = color))
    for i in np.where(self.flux > ylim[1])[0]:
      if i in bnmask:
        color = "#ffcccc"
      elif i in self.outmask:
        color = "#cccccc"
      else:
        color = "#ccccff"
      ax.annotate('', xy=(self.time[i], ylim[1]), xycoords = 'data',
                  xytext = (0, -15), textcoords = 'offset points',
                  arrowprops=dict(arrowstyle = "-|>", color = color))

    if show:
      pl.show()