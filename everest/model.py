#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`model.py` - De-trending models
---------------------------------------

This module contains the generic models used to de-trend light curves for the various
supported missions. Most of the functionality is implemented in :py:class:`Model`, and
specific de-trending methods are implemented as subclasses.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .config import EVEREST_DAT
from .missions import Missions
from .utils import InitLog, Formatter, AP_SATURATED_PIXEL, AP_COLLAPSED_PIXEL
from .math import Chunks, RMS, CDPP6, SavGol, Interpolate
from .data import Season, GetData, GetNeighbors, Breakpoint, HasShortCadence, MakeFITS
from .gp import GetCovariance, GetKernelParams
from .transit import Transit
from .dvs import DVS1, DVS2
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

__all__ = ['Model', 'Inject', 'MaskedInject', 'rPLD', 'nPLD']

class Model(object):
  '''
  A generic *PLD* model with scalar matrix *L2* regularization. Includes functionality
  for loading pixel-level light curves, identifying outliers, generating the data
  covariance matrix, computing the regularized pixel model, and plotting the results.
  Specific models are implemented as subclasses.
  
  **General:**
  
  :param ID: The target star ID (*EPIC*, *KIC*, or *TIC* number, for instance)
  :param bool clobber: Overwrite existing :py:obj:`everest` models? Default :py:obj:`False`
  :param bool clobber_tpf: Download and overwrite the saved raw TPF data? Default :py:obj:`False`
  :param bool debug: De-trend in debug mode? If :py:obj:`True`, prints all output to screen and \
                     enters :py:obj:`pdb` post-mortem mode for debugging when an error is raised.
                     Default :py:obj:`False`
  :param bool make_fits: Generate a *FITS* file at the end of rthe run? Default :py:obj:`False`
  :param str mission: The name of the mission. Default `k2`
  
  **Model:**
  
  :param str aperture_name: The name of the aperture to use. These are defined in the datasets and are \
                            mission specific. Defaults to the mission default
  :param int bpad: When light curve breakpoints are set, the light curve chunks must be stitched together \
                   at the end. To prevent kinks and/or discontinuities, the chunks are made to overlap by \
                   :py:obj:`bpad` cadences on either end. The chunks are then mended and the overlap is \
                   discarded. Default 100
  :param bool breakpoint: Add a light curve breakpoint when de-trending? If :py:obj:`True`, splits the \
                          light curve into two chunks and de-trends each one separately, then stitches them \
                          back and the end. This is useful for missions like *K2*, where the light curve noise \
                          properties are very different at the beginning and end of each campaign. The cadences \
                          at which breakpoints are inserted are specified in the :py:func:`Breakpoint` function \
                          of each mission. Default :py:obj:`True`
  :param int cdivs: The number of light curve subdivisions when cross-validating. During each iteration, \
                    one of these subdivisions will be masked and used as the validation set. Default 3
  :param int giter: The number of iterations when optimizing the GP. During each iteration, the minimizer \
                    is initialized with a perturbed guess; after :py:obj:`giter` iterations, the step with \
                    the highest likelihood is kept. Default 3
  :param float gp_factor: When computing the initial kernel parameters, the red noise amplitude is set to \
                          the standard deviation of the data times this factor. Larger values generally \
                          help with convergence, particularly for very variable stars. Default 100
  :param array_like kernel_params: The initial value of the :py:obj:`Matern-3/2` kernel parameters \
                                   (white noise amplitude in flux units, red noise amplitude in flux units, \
                                   and timescale in days). Default :py:obj:`None` (determined from the data)
  :param array_like lambda_arr: The array of :math:`\Lambda` values to iterate over during the \
                                cross-validation step. :math:`\Lambda` is the regularization parameter,
                                or the standard deviation of \
                                the Gaussian prior on the weights for each order of PLD. \
                                Default ``10 ** np.arange(0,18,0.5)``
  :param float leps: The fractional tolerance when optimizing :math:`\Lambda`. The chosen value of \
                     :math:`\Lambda` will be within this amount of the minimum of the CDPP curve. \
                     Default 0.05
  :param int max_pixels: The maximum number of pixels. Very large apertures are likely to cause memory \
                         errors, particularly for high order PLD. If the chosen aperture exceeds this many \
                         pixels, a different aperture is chosen from the dataset. If no apertures with fewer \
                         than this many pixels are available, an error is thrown. Default 75
  :param bool optimize_gp: Perform the GP optimization steps? Default :py:obj:`True`
  :param float osigma: The outlier standard deviation threshold. Default 5
  :param int oiter: The maximum number of steps taken during iterative sigma clipping. Default 10
  :param int pld_order: The pixel level decorrelation order. Default `3`. Higher orders may cause memory errors
  :param bool recursive: Calculate the fractional pixel flux recursively? If :py:obj:`False`, \
                         always computes the fractional pixel flux as :math:`f_{ij}/\sum_i{f_{ij}}`. \
                         If :py:obj:`True`, uses the current de-trended flux as the divisor in an attempt \
                         to minimize the amount of instrumental information that is divided out in this \
                         step: :math:`f_{ij}/(\sum_i{f_{ij}} - M)`. Default :py:obj:`True`
  :param str saturated_aperture_name: If the target is found to be saturated, de-trending is performed \
                                      on this aperture instead. Defaults to the mission default
  :param float saturation_tolerance: The tolerance when determining whether or not to collapse a column \
                                     in the aperture. The column collapsing is implemented in the individual \
                                     mission modules. Default -0.1, i.e., if a target is 10% shy of the \
                                     nominal saturation level, it is considered to be saturated.
  :param int sc_bpad: Same as :py:obj:`bpad`, but for the short cadence data. Default 3000
  
  '''

  def __init__(self, ID, **kwargs):
    '''
    
    '''
        
    # Initialize logging
    self.ID = ID
    self.recursive = kwargs.get('recursive', True)
    self.make_fits = kwargs.get('make_fits', False)
    self.mission = kwargs.get('mission', 'k2')
    self.has_sc = HasShortCadence(self.ID, season = self.season)
    self.clobber = kwargs.get('clobber', False)
    self.debug = kwargs.get('debug', False)
    self.is_parent = kwargs.get('is_parent', False)
    if not self.is_parent:
      screen_level = kwargs.get('screen_level', logging.CRITICAL)
      log_level = kwargs.get('log_level', logging.DEBUG)
      InitLog(self.logfile, log_level, screen_level, self.debug)
      log.info("Initializing %s model for %d." % (self.name, self.ID))
    
    # Read general model kwargs
    self.lambda_arr = kwargs.get('lambda_arr', 10 ** np.arange(0,18,0.5))
    if self.lambda_arr[0] != 0:
      self.lambda_arr = np.append(0, self.lambda_arr)
    self.leps = kwargs.get('leps', 0.05)
    self.osigma = kwargs.get('osigma', 5)
    self.oiter = kwargs.get('oiter', 10)
    self.cdivs = kwargs.get('cdivs', 3)
    self.giter = kwargs.get('giter', 3)
    self.optimize_gp = kwargs.get('optimize_gp', True)
    self.kernel_params = kwargs.get('kernel_params', None)    
    self.clobber_tpf = kwargs.get('clobber_tpf', False)
    self.bpad = kwargs.get('bpad', 100)
    self.sc_bpad = kwargs.get('sc_bpad', 3000)
    self.aperture_name = kwargs.get('aperture', None)
    self.saturated_aperture_name = kwargs.get('saturated_aperture', None)
    self.max_pixels = kwargs.get('max_pixels', 75)
    self.saturation_tolerance = kwargs.get('saturation_tolerance', -0.1)
    self.gp_factor = kwargs.get('gp_factor', 100.)
    
    # Handle breakpointing. The breakpoint is the *last* index of the first 
    # light curve chunk. The code is (for the most part) written to allow for
    # multiple breakpoints, but the plotting can currently only handle one.
    if kwargs.get('breakpoint', True):
      self.breakpoints = [Breakpoint(self.ID), 999999]
      if self.breakpoints[0] is None:
        raise ValueError('Invalid breakpoint set in `%d.Breakpoint()`.' % (self.mission))
    else:
      self.breakpoints = [999999]
    nseg = len(self.breakpoints)

    # Get the pld order
    pld_order = kwargs.get('pld_order', 3)
    assert (pld_order > 0), "Invalid value for the de-trending order."
    self.pld_order = pld_order

    # Initialize model params 
    self.lam_idx = -1
    self.lam = [[1e5] + [None for i in range(self.pld_order - 1)] for b in range(nseg)]
    self._A = [[None for i in range(self.pld_order)] for b in range(nseg)]
    self._B = [[None for i in range(self.pld_order)] for b in range(nseg)]
    self._mK = [None for b in range(nseg)]
    self._f = [None for b in range(nseg)]
    self._X = [None for i in range(self.pld_order)]
    self.XNeighbors = [None for i in range(self.pld_order)]
    self.sc_X1N = None
    self.cdpp6_arr = np.array([np.nan for b in range(nseg)])
    self.cdppr_arr = np.array([np.nan for b in range(nseg)])
    self.cdppv_arr = np.array([np.nan for b in range(nseg)])
    self.cdpp6 = np.nan
    self.cdppr = np.nan
    self.cdppv = np.nan
    self.gppp = np.nan
    self.neighbors = []
    self.loaded = False
    self._weights = None
    
    # Initialize plotting
    self.dvs1 = DVS1(nseg > 1, pld_order = self.pld_order)
    self.dvs2 = DVS2(nseg > 1)
  
  @property
  def name(self):
    '''
    Returns the name of the current :py:class:`Model` subclass.
    
    '''
    
    return self.__class__.__name__
      
  @name.setter
  def name(self, value):
    '''
    
    '''
    
    raise NotImplementedError("Can't set this property.") 
  
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
      self._season = Season(self.ID, self.mission)
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
    
    else:
      
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
    
  def _precompute(self, mask, b):
    '''
    Pre-compute the matrices :py:obj:`A` and :py:obj:`B` (cross-validation step only)
    for chunk :py:obj:`b`.
    
    '''
    
    # Get current chunk and mask outliers
    m1 = self.get_masked_chunk(b)
    flux = self.fraw[m1]
    K = self.K[m1][:,m1]
    med = np.nanmedian(self.fraw)
    
    # Now mask the validation set
    M = lambda x, axis = 0: np.delete(x, mask, axis = axis)
    m2 = M(m1)
    mK = M(M(K, axis = 0), axis = 1)
    f = M(flux) - med
    
    # Pre-compute the matrices
    A = [None for i in range(self.pld_order)]
    B = [None for i in range(self.pld_order)] 
    for n in range(self.pld_order):
      if self.X[n] is not None:
        A[n] = np.dot(self.X[n][m2], self.X[n][m2].T)
        B[n] = np.dot(self.X[n][m1], self.X[n][m2].T)
    
    return A, B, mK, f
    
  def _compute(self, b, A, B, mK, f):
    '''
    Compute the model (cross-validation step only) for chunk :py:obj:`b`.
    
    '''

    A = np.sum([l * a for l, a in zip(self.lam[b], A) if l is not None], axis = 0)
    B = np.sum([l * b for l, b in zip(self.lam[b], B) if l is not None], axis = 0)
    W = np.linalg.solve(mK + A, f)
    model = np.dot(B, W)
    model -= np.nanmedian(model)
    
    return model
  
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
  
  def get_outliers(self):
    '''
    Performs iterative sigma clipping to get outliers.
    
    '''
            
    log.info("Clipping outliers...")
    log.info('Iter %d/%d: %d outliers' % (0, self.oiter, len(self.outmask)))
    M = lambda x: np.delete(x, np.concatenate([self.nanmask, self.badmask]), axis = 0)
    t = M(self.time)
    outmask = [np.array([-1]), np.array(self.outmask)]
    
    # Loop as long as the last two outlier arrays aren't equal
    while not np.array_equal(outmask[-2], outmask[-1]):

      # Check if we've done this too many times
      if len(outmask) - 1 > self.oiter:
        log.error('Maximum number of iterations in ``get_outliers()`` exceeded. Skipping...')
        break
    
      # Check if we're going in circles
      if np.any([np.array_equal(outmask[-1], i) for i in outmask[:-1]]):
        log.error('Function ``get_outliers()`` is going in circles. Skipping...')
        break
      
      # Compute the model to get the flux
      self.compute()
    
      # Get the outliers
      f = SavGol(M(self.flux))
      med = np.nanmedian(f)
      MAD = 1.4826 * np.nanmedian(np.abs(f - med))
      inds = np.where((f > med + self.osigma * MAD) | (f < med - self.osigma * MAD))[0]
      
      # Project onto unmasked time array
      inds = np.array([np.argmax(self.time == t[i]) for i in inds])
      self.outmask = np.array(inds, dtype = int)
      
      # Add them to the running list
      outmask.append(np.array(inds))
      
      # Log
      log.info('Iter %d/%d: %d outliers' % (len(outmask) - 2, self.oiter, len(self.outmask)))

  def optimize_lambda(self, validation):
    '''
    Returns the index of :py:attr:`self.lambda_arr` that minimizes the validation scatter
    in the segment with minimum at the lowest value of :py:obj:`lambda`, with
    fractional tolerance :py:attr:`self.leps`.
    
    :param numpy.ndarray validation: The scatter in the validation set as a function of :py:obj:`lambda`
    
    '''
    
    maxm = 0
    minr = len(validation)
    for n in range(validation.shape[1]):
      m = np.nanargmin(validation[:,n])
      if m > maxm:
        maxm = m
      r = np.where((validation[:,n] - validation[m,n]) / 
                    validation[m,n] <= self.leps)[0][-1]
      if r < minr:
        minr = r
    return min(maxm, minr)

  def cross_validate(self, ax, info = ''):
    '''
    Cross-validate to find the optimal value of :py:obj:`lambda`.
    
    :param ax: The current :py:obj:`matplotlib.pyplot` axis instance to plot the \
               cross-validation results.
    :param str info: The label to show in the bottom right-hand corner of the plot. Default `''`
    
    '''
    
    # Loop over all chunks
    ax = np.atleast_1d(ax)
    for b, brkpt in enumerate(self.breakpoints):
    
      log.info("Cross-validating chunk %d/%d..." % (b + 1, len(self.breakpoints)))
      med_training = np.zeros_like(self.lambda_arr)
      med_validation = np.zeros_like(self.lambda_arr)
        
      # Mask for current chunk 
      m = self.get_masked_chunk(b)
        
      # Mask transits and outliers
      time = self.time[m]
      flux = self.fraw[m]
      ferr = self.fraw_err[m]
      med = np.nanmedian(self.fraw)
    
      # The precision in the validation set
      validation = [[] for k, _ in enumerate(self.lambda_arr)]
    
      # The precision in the training set
      training = [[] for k, _ in enumerate(self.lambda_arr)]
    
      # Setup the GP
      _, amp, tau = self.kernel_params
      gp = george.GP(amp ** 2 * george.kernels.Matern32Kernel(tau ** 2))
      gp.compute(time, ferr)
    
      # The masks
      masks = list(Chunks(np.arange(0, len(time)), len(time) // self.cdivs))
    
      # Loop over the different masks
      for i, mask in enumerate(masks):
      
        log.info("Section %d/%d..." % (i + 1, len(masks)))
      
        # Pre-compute (training set)
        pre_t = self._precompute([], b)

        # Pre-compute (validation set)
        pre_v = self._precompute(mask, b)
    
        # Iterate over lambda
        for k, lam in enumerate(self.lambda_arr):
      
          # Update the lambda matrix
          self.lam[b][self.lam_idx] = lam
      
          # Training set. Note that we're computing the MAD, not the
          # standard deviation, as this handles extremely variable
          # stars much better!
          model = self._compute(b, *pre_t)
          gpm, _ = gp.predict(flux - model - med, time[mask])
          fdet = (flux - model)[mask] - gpm
          scatter = 1.e6 * (1.4826 * np.nanmedian(np.abs(fdet / med - 
                                                  np.nanmedian(fdet / med))) /
                                                  np.sqrt(len(mask)))
          training[k].append(scatter)
      
          # Validation set
          model = self._compute(b, *pre_v)
          gpm, _ = gp.predict(flux - model - med, time[mask])
          fdet = (flux - model)[mask] - gpm
          scatter = 1.e6 * (1.4826 * np.nanmedian(np.abs(fdet / med - 
                                                  np.nanmedian(fdet / med))) /
                                                  np.sqrt(len(mask)))
          validation[k].append(scatter)
      
      # Finalize
      training = np.array(training)
      validation = np.array(validation)
      for k, _ in enumerate(self.lambda_arr):

        # Take the mean
        med_validation[k] = np.nanmean(validation[k])
        med_training[k] = np.nanmean(training[k])
            
      # Compute best model
      i = self.optimize_lambda(validation)
      v_best = med_validation[i]
      t_best = med_training[i]
      self.cdppv_arr[b] = v_best / t_best
      self.lam[b][self.lam_idx] = self.lambda_arr[i]
      log.info("Found optimum solution at log(lambda) = %.1f." % np.log10(self.lam[b][self.lam_idx]))
      self.compute()

      # Plotting hack: first x tick will be -infty
      lambda_arr = np.array(self.lambda_arr)
      lambda_arr[0] = 10 ** (np.log10(lambda_arr[1]) - 3)
    
      # Plot cross-val
      for n in range(len(masks)):
        ax[b].plot(np.log10(lambda_arr), validation[:,n], 'r-', alpha = 0.3)
      ax[b].plot(np.log10(lambda_arr), med_training, 'b-', lw = 1., alpha = 1)
      ax[b].plot(np.log10(lambda_arr), med_validation, 'r-', lw = 1., alpha = 1)            
      ax[b].axvline(np.log10(self.lam[b][self.lam_idx]), color = 'k', ls = '--', lw = 0.75, alpha = 0.75)
      ax[b].axhline(v_best, color = 'k', ls = '--', lw = 0.75, alpha = 0.75)
      ax[b].set_ylabel(r'Scatter (ppm)', fontsize = 5)
      hi = np.max(validation[0])
      lo = np.min(training)
      rng = (hi - lo)
      ax[b].set_ylim(lo - 0.15 * rng, hi + 0.15 * rng)
      if rng > 2:
        ax[b].get_yaxis().set_major_formatter(Formatter.CDPP)
        ax[b].get_yaxis().set_major_locator(MaxNLocator(integer = True))
      elif rng > 0.2:
        ax[b].get_yaxis().set_major_formatter(Formatter.CDPP1F)
      else:
        ax[b].get_yaxis().set_major_formatter(Formatter.CDPP2F)
        
      # Fix the x ticks
      xticks = [np.log10(lambda_arr[0])] + list(np.linspace(np.log10(lambda_arr[1]), np.log10(lambda_arr[-1]), 6))
      ax[b].set_xticks(xticks)
      ax[b].set_xticklabels(['' for x in xticks])
      pad = 0.01 * (np.log10(lambda_arr[-1]) - np.log10(lambda_arr[0]))
      ax[b].set_xlim(np.log10(lambda_arr[0]) - pad, np.log10(lambda_arr[-1]) + pad)
      ax[b].annotate('%s.%d' % (info, b), xy = (0.02, 0.025), xycoords = 'axes fraction', 
                     ha = 'left', va = 'bottom', fontsize = 8, alpha = 0.5, 
                     fontweight = 'bold')
    
    # Tidy up
    if len(ax) > 1:
      ax[0].xaxis.set_ticks_position('top')
    for axis in ax[1:]:
      axis.spines['top'].set_visible(False)
      axis.xaxis.set_ticks_position('bottom')
    
    # A hack to mark the first xtick as -infty
    labels = ['%.1f' % x for x in xticks]
    labels[0] = r'$-\infty$'
    ax[-1].set_xticklabels(labels) 
    ax[-1].set_xlabel(r'Log $\Lambda$', fontsize = 5)
  
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

  def finalize(self):
    '''
    This method is called at the end of the de-trending, prior to plotting the final results.
    Subclass it to add custom functionality to individual models.
    
    '''
    
    pass
  
  def get_ylim(self):
    '''
    Computes the ideal y-axis limits for the light curve plot. Attempts to set
    the limits equal to those of the raw light curve, but if more than 1% of the
    flux lies either above or below these limits, auto-expands to include those
    points. At the end, adds 5% padding to both the top and the bottom.
    
    '''
    
    bn = np.array(list(set(np.concatenate([self.badmask, self.nanmask]))), dtype = int)
    fraw = np.delete(self.fraw, bn)
    lo, hi = fraw[np.argsort(fraw)][[3,-3]]
    flux = np.delete(self.flux, bn)
    fsort = flux[np.argsort(flux)]
    if fsort[int(0.01 * len(fsort))] < lo:
      lo = fsort[int(0.01 * len(fsort))]
    if fsort[int(0.99 * len(fsort))] > hi:
      hi = fsort[int(0.99 * len(fsort))]
    pad = (hi - lo) * 0.05
    ylim = (lo - pad, hi + pad)
    return ylim
    
  def plot(self, ax, info_left = '', info_right = '', color = 'b'):
    '''
    Plots the current light curve. This is called at several stages to plot the
    de-trending progress as a function of the different *PLD* orders.
    
    :param ax: The current :py:obj:`matplotlib.pyplot` axis instance
    :param str info_left: Information to display at the left of the plot. Default `''`
    :param str info_right: Information to display at the right of the plot. Default `''`
    :param str color: The color of the data points. Default `'b'`
    
    '''

    # Plot
    ax.plot(self.apply_mask(self.time), self.apply_mask(self.flux), ls = 'none', marker = '.', color = color, markersize = 2, alpha = 0.5)
    ylim = self.get_ylim()
    
    # Plot the outliers
    bnmask = np.array(list(set(np.concatenate([self.badmask, self.nanmask]))), dtype = int)
    O1 = lambda x: x[self.outmask]
    O2 = lambda x: x[bnmask]
    ax.plot(O1(self.time), O1(self.flux), ls = 'none', color = "#777777", marker = '.', markersize = 2, alpha = 0.5)
    ax.plot(O2(self.time), O2(self.flux), 'r.', markersize = 2, alpha = 0.25)
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
    
    # Plot the breakpoints
    for brkpt in self.breakpoints[:-1]:
      ax.axvline(self.time[brkpt], color = 'r', ls = '--', alpha = 0.5)
    
    # Appearance
    ax.annotate('%.2f ppm' % self.cdpp6_arr[0], xy = (0.02, 0.975), xycoords = 'axes fraction', 
                ha = 'left', va = 'top', fontsize = 8)
    if len(self.cdpp6_arr) == 2:
      ax.annotate('%.2f ppm' % self.cdpp6_arr[1], xy = (0.98, 0.975), xycoords = 'axes fraction', 
                  ha = 'right', va = 'top', fontsize = 8)
    ax.annotate(info_right, xy = (0.98, 0.025), xycoords = 'axes fraction', 
                ha = 'right', va = 'bottom', fontsize = 10, alpha = 0.5, 
                fontweight = 'bold')            
    ax.annotate(info_left, xy = (0.02, 0.025), xycoords = 'axes fraction', 
                ha = 'left', va = 'bottom', fontsize = 8)     
    ax.set_xlabel(r'Time (BJD - 2454833)', fontsize = 5)
    ax.margins(0.01, 0.1)
    ax.set_ylim(*ylim)
    ax.get_yaxis().set_major_formatter(Formatter.Flux)
  
  def plot_info(self, dvs):
    '''
    Plots miscellaneous de-trending information on the data validation summary figure.
    
    :param dvs: A :py:class:`dvs.DVS1` or :py:class:`dvs.DVS2` figure instance
    
    '''
    
    axl, axc, axr = dvs.title()
    axc.annotate("%s %d" % (Missions[self.mission].IDSTRING, self.ID),
                 xy = (0.5, 0.5), xycoords = 'axes fraction', 
                 ha = 'center', va = 'center', fontsize = 18)
    
    axc.annotate(r"%.2f ppm $\rightarrow$ %.2f ppm" % (self.cdppr, self.cdpp6),
                 xy = (0.5, 0.2), xycoords = 'axes fraction',
                 ha = 'center', va = 'center', fontsize = 8, color = 'k',
                 fontstyle = 'italic')
    
    axl.annotate("%s %s%02d: %s" % (self.mission.upper(), 
                 Missions[self.mission].SEASONCHAR, self.season, self.name),
                 xy = (0.5, 0.5), xycoords = 'axes fraction', 
                 ha = 'center', va = 'center', fontsize = 12,
                 color = 'k')
    
    axl.annotate(self.aperture_name if len(self.neighbors) == 0 else "%s, %d neighbors" % (self.aperture_name, len(self.neighbors)),
                 xy = (0.5, 0.2), xycoords = 'axes fraction',
                 ha = 'center', va = 'center', fontsize = 8, color = 'k',
                 fontstyle = 'italic')
    
    axr.annotate("%s %.3f" % (Missions[self.mission].MAGSTRING, self.mag),
                 xy = (0.5, 0.5), xycoords = 'axes fraction', 
                 ha = 'center', va = 'center', fontsize = 12,
                 color = 'k')
    
    axr.annotate(r"GP %.3f ppm" % (self.gppp),
                 xy = (0.5, 0.2), xycoords = 'axes fraction',
                 ha = 'center', va = 'center', fontsize = 8, color = 'k',
                 fontstyle = 'italic')
  
  def plot_final(self):
    '''
    Plots the final de-trended light curve.
    
    '''
    
    # Plot the light curve
    ax = self.dvs1.top_left()
    bnmask = np.array(list(set(np.concatenate([self.badmask, self.nanmask]))), dtype = int)
    M = lambda x: np.delete(x, bnmask)
    ax.plot(M(self.time), M(self.flux), ls = 'none', marker = '.', color = 'k', markersize = 2, alpha = 0.3)

    # Plot the GP
    _, amp, tau = self.kernel_params
    gp = george.GP(amp ** 2 * george.kernels.Matern32Kernel(tau ** 2))
    gp.compute(self.apply_mask(self.time), self.apply_mask(self.fraw_err))
    med = np.nanmedian(self.apply_mask(self.flux))
    y, _ = gp.predict(self.apply_mask(self.flux) - med, self.time)
    y += med
    ax.plot(M(self.time), M(y), 'r-', lw = 0.5, alpha = 0.5)
    
    # Appearance
    ax.annotate('Final', xy = (0.98, 0.025), xycoords = 'axes fraction', 
                ha = 'right', va = 'bottom', fontsize = 10, alpha = 0.5, 
                fontweight = 'bold') 
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
    
    # Compute the CDPP of the GP-detrended flux
    self.gppp = CDPP6(self.apply_mask(self.flux - y + med))
  
  def plot_aperture(self, dvs):
    '''
    Plots the aperture and the pixel images at the beginning, middle, and end of 
    the time series. Also plots a high resolution image of the target, if available.
    
    :param dvs: A :py:class:`dvs.DVS1` or :py:class:`dvs.DVS2` figure instance
    
    '''
    
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
      ax = dvs.top_right()
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
                  ha = 'center', va = 'top', size = 8, color = 'w')
      if i == 1:
        for source in self.nearby:
          ax.annotate('%.1f' % source['mag'], 
                      xy = (source['x'] - source['x0'], source['y'] - source['y0']), 
                      ha = 'center', va = 'center', size = 6, color = 'w', fontweight = 'bold')    
      
    # Plot hi res image
    if self.hires is not None:
      ax = dvs.top_right()
      ax.imshow(self.hires, aspect = 'auto', extent = (-0.5, nx - 0.5, -0.5, ny - 0.5), interpolation = 'bicubic', cmap = plasma)
      ax.contour(highres, levels=[0.5], extent=extent, origin='lower', colors='r', linewidths=1)
      ax.axis('off')
      ax.set_xlim(-0.7, nx - 0.3)
      ax.set_ylim(-0.7, ny - 0.3)
      ax.annotate('hires', xy = (0.5, 0.975), xycoords = 'axes fraction',
                  ha = 'center', va = 'top', size = 8, color = 'w')
  
  def plot_weights(self):
    '''
    Plots the *PLD* weights on the CCD for each of the *PLD* orders.
    
    '''
        
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
    ax, cax = self.dvs2.weights_grid()
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
  
  def plot_page2(self):
    '''
    Plots the second page of the data validation summary.
    
    '''
    
    # Plot the short cadence light curve (if available). Otherwise,
    # plot the raw light curve.
    ax1 = self.dvs2.lc1()
    
    if self.has_sc:
      bnmask = np.array(list(set(np.concatenate([self.sc_badmask, self.sc_nanmask]))), dtype = int)
      M = lambda x: np.delete(x, bnmask)
      ax1.plot(M(self.sc_time), M(self.sc_flux), ls = 'none', marker = '.', color = 'k', markersize = 2, alpha = 0.03, zorder = -1)
      ax1.set_rasterization_zorder(0)
      ax1.annotate('SC', xy = (0.98, 0.025), xycoords = 'axes fraction', 
                  ha = 'right', va = 'bottom', fontsize = 10, alpha = 0.5, 
                  fontweight = 'bold') 
      ax1.margins(0.01, 0.1)  
      flux = np.delete(self.sc_flux, bnmask)
      N = int(0.995 * len(flux))
      hi, lo = flux[np.argsort(flux)][[N,-N]]
      fsort = flux[np.argsort(flux)]
      pad = (hi - lo) * 0.1
      ylim = (lo - pad, hi + pad)   
      ax1.set_ylim(ylim)   
      ax1.get_yaxis().set_major_formatter(Formatter.Flux)  
    else:
      ax1.plot(self.time, self.fraw, ls = 'none', marker = '.', color = 'k', markersize = 2, alpha = 0.5)
      ax1.annotate('Raw', xy = (0.98, 0.025), xycoords = 'axes fraction', 
                  ha = 'right', va = 'bottom', fontsize = 10, alpha = 0.5, 
                  fontweight = 'bold') 
      ax1.margins(0.01, 0.1)  
      bnmask = np.array(list(set(np.concatenate([self.badmask, self.nanmask]))), dtype = int)
      flux = np.delete(self.flux, bnmask)
      N = int(0.995 * len(flux))
      hi, lo = flux[np.argsort(flux)][[N,-N]]
      fsort = flux[np.argsort(flux)]
      pad = (hi - lo) * 0.1
      ylim = (lo - pad, hi + pad)   
      ax1.set_ylim(ylim)   
      ax1.get_yaxis().set_major_formatter(Formatter.Flux) 

    # Plot the long cadence light curve
    ax2 = self.dvs2.lc2()
    bnmask = np.array(list(set(np.concatenate([self.badmask, self.nanmask]))), dtype = int)
    M = lambda x: np.delete(x, bnmask)
    ax2.plot(M(self.time), M(self.flux), ls = 'none', marker = '.', color = 'k', markersize = 2, alpha = 0.3)
    ax2.annotate('LC', xy = (0.98, 0.025), xycoords = 'axes fraction', 
                ha = 'right', va = 'bottom', fontsize = 10, alpha = 0.5, 
                fontweight = 'bold') 
    ax2.margins(0.01, 0.1)          
    ax2.set_ylim(ylim)   
    ax2.get_yaxis().set_major_formatter(Formatter.Flux) 
    
    # Plot the PLD weights
    self.plot_weights()
    
  def load_tpf(self):
    '''
    Loads the target pixel file.
    
    '''
    
    if not self.loaded:
      data = GetData(self.ID, self.mission, season = self.season, clobber = self.clobber_tpf, 
                     aperture_name = self.aperture_name, 
                     saturated_aperture_name = self.saturated_aperture_name, 
                     max_pixels = self.max_pixels,
                     saturation_tolerance = self.saturation_tolerance)
      self.cadn = data.cadn
      self.time = data.time
      self.model = np.zeros_like(self.time)
      self.fpix = data.fpix
      self.fraw = np.sum(self.fpix, axis = 1)
      self.fpix_err = data.fpix_err
      self.fraw_err = np.sqrt(np.sum(self.fpix_err ** 2, axis = 1))
      self.nanmask = data.nanmask
      self.badmask = data.badmask
      self.transitmask = np.array([], dtype = int)
      self.outmask = np.array([], dtype = int)
      self.aperture = data.aperture
      self.aperture_name = data.aperture_name
      self.apertures = data.apertures
      self.quality = data.quality
      self.Xpos = data.Xpos
      self.Ypos = data.Ypos
      self.mag = data.mag
      self.pixel_images = data.pixel_images
      self.nearby = data.nearby
      self.hires = data.hires
      self.saturated = data.saturated
      self.meta = data.meta
      self.bkg = data.bkg
      if self.has_sc:
        self.sc_cadn = data.sc_cadn
        self.sc_time = data.sc_time
        self.sc_fpix = data.sc_fpix
        self.sc_fraw = np.sum(self.sc_fpix, axis = 1)
        self.sc_fpix_err = data.sc_fpix_err
        self.sc_fraw_err = np.sqrt(np.sum(self.sc_fpix_err ** 2, axis = 1))
        self.sc_nanmask = data.sc_nanmask
        self.sc_badmask = data.sc_badmask
        self.sc_transitmask = np.array([], dtype = int)
        self.sc_outmask = np.array([], dtype = int)
        self.sc_Xpos = data.sc_Xpos
        self.sc_Ypos = data.sc_Ypos
        self.sc_meta = data.sc_meta
        self.sc_model = np.zeros_like(self.sc_time)
        self.sc_quality = data.sc_quality
        self.sc_bkg = data.sc_bkg
      else:
        self.sc_cadn = None
        self.sc_time = None
        self.sc_fpix = None
        self.sc_fraw = None
        self.sc_fpix_err = None
        self.sc_fraw_err = None
        self.sc_nanmask = None
        self.sc_badmask = None
        self.sc_transitmask = None
        self.sc_outmask = None
        self.sc_Xpos = None
        self.sc_Ypos = None
        self.sc_meta = None
        self.sc_model = None
        self.sc_quality = None
        self.sc_bkg = None
      
      # Update the last breakpoint to the correct value
      self.breakpoints[-1] = len(self.time) - 1
      self.loaded = True
  
  def load_model(self, name = None):
    '''
    Loads a saved version of the model.
    
    '''
    
    if self.clobber:
      return False
    
    if name is None:
      name = self.name    
    file = os.path.join(self.dir, '%s.npz' % name)
    if os.path.exists(file):
      if not self.is_parent: 
        log.info("Loading '%s.npz'..." % name)
      try:
        data = np.load(file)
        for key in data.keys():
          try:
            setattr(self, key, data[key][()])
          except NotImplementedError:
            pass
        self.K = GetCovariance(self.kernel_params, self.time, self.fraw_err)
        self.get_X()
        pl.close()
        return True
      except:
        log.warn("Error loading '%s.npz'." % name)
        os.rename(file, file + '.bad')
    
    if self.is_parent:
      raise Exception('Unable to load `%s` model for target %d.' % (self.name, self.ID))
    
    return False

  def save_model(self):
    '''
    Saves all of the de-trending information to disk in an `npz` file
    and saves the DVS as a `pdf`.
    
    '''
    
    # Save the data
    log.info("Saving data to '%s.npz'..." % self.name)
    d = dict(self.__dict__)
    d.pop('_X', None)
    d.pop('_weights', None)
    d.pop('_A', None)
    d.pop('_B', None)
    d.pop('_f', None)
    d.pop('_mK', None)
    d.pop('K', None)
    d.pop('dvs1', None)
    d.pop('dvs2', None)
    d.pop('clobber', None)
    d.pop('clobber_tpf', None)
    d.pop('debug', None)
    np.savez(os.path.join(self.dir, self.name + '.npz'), **d)
    
    # Save the DVS
    pdf = PdfPages(os.path.join(self.dir, self.name + '.pdf'))
    pdf.savefig(self.dvs1.fig)
    pl.close(self.dvs1.fig)
    pdf.savefig(self.dvs2.fig)
    pl.close(self.dvs2.fig)
    d = pdf.infodict()
    d['Title'] = 'EVEREST: %s de-trending of %s %d' % (self.name, Missions[self.mission].IDSTRING, self.ID)
    d['Author'] = 'Rodrigo Luger'
    pdf.close()
    
  def exception_handler(self, pdb):
    '''
    A custom exception handler.
    
    :param pdb: If :py:obj:`True`, enters PDB post-mortem mode for debugging.
    
    '''
    
    # Grab the exception
    exctype, value, tb = sys.exc_info()
    
    # Log the error and create a .err file
    errfile = os.path.join(self.dir, self.name + '.err')
    with open(errfile, 'w') as f:
      for line in traceback.format_exception_only(exctype, value):
        l = line.replace('\n', '')
        log.error(l)
        print(l, file = f)
      for line in traceback.format_tb(tb):
        l = line.replace('\n', '')
        log.error(l)
        print(l, file = f)
    
    # Re-raise?
    if pdb:
      raise
  
  def update_gp(self):
    '''
    Calls :py:func:`gp.GetKernelParams` to optimize the GP and obtain the
    covariance matrix for the regression.
    
    '''
    
    self.kernel_params = GetKernelParams(self.time, self.flux, self.fraw_err, 
                                         mask = self.mask, guess = self.kernel_params, 
                                         giter = self.giter)
    self.K = GetCovariance(self.kernel_params, self.time, self.fraw_err)
  
  def init_kernel(self):
    '''
    Initializes the covariance matrix with a guess at the GP kernel parameters.
    
    '''
    
    if self.kernel_params is None:
      X = self.apply_mask(self.fpix / self.flux.reshape(-1, 1))
      y = self.apply_mask(self.flux) - np.dot(X, np.linalg.solve(np.dot(X.T, X), np.dot(X.T, self.apply_mask(self.flux))))      
      white = np.nanmedian([np.nanstd(c) for c in Chunks(y, 13)])
      amp = self.gp_factor * np.nanstd(y)
      tau = 30.0
      self.kernel_params = [white, amp, tau]
    self.K = GetCovariance(self.kernel_params, self.time, self.fraw_err)
  
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
  
  def run(self):
    '''
    Runs the de-trending step.
    
    '''
    
    try:
          
      # Load raw data
      log.info("Loading target data...")
      self.load_tpf()
      self.plot_aperture(self.dvs1)     
      self.plot_aperture(self.dvs2) 
      self.init_kernel()
      M = self.apply_mask(np.arange(len(self.time)))
      self.cdppr_arr = self.get_cdpp_arr()
      self.cdpp6_arr = np.array(self.cdppr_arr)
      self.cdppv_arr = np.array(self.cdppr_arr)
      self.cdppr = self.get_cdpp()
      self.cdpp6 = self.cdppr
      self.cdppv = self.cdppv

      log.info("%s (Raw): CDPP6 = %s" % (self.name, self.cdpps))
      self.plot(self.dvs1.left(), info_right = 'Raw', color = 'k')
      
      # Loop
      for n in range(self.pld_order):
        self.lam_idx += 1
        self.get_X()
        self.get_outliers()
        if n > 0 and self.optimize_gp:
          self.update_gp()
        self.cross_validate(self.dvs1.right(), info = 'CV%d' % n)
        self.cdpp6_arr = self.get_cdpp_arr()
        self.cdppv_arr *= self.cdpp6_arr
        self.cdpp6 = self.get_cdpp()
        self.cdppv = np.mean(self.cdppv_arr)
        log.info("%s (%d/%d): CDPP = %s" % (self.name, n + 1, self.pld_order, self.cdpps))
        self.plot(self.dvs1.left(), info_right= 'LC%d' % (n + 1), info_left = '%d outliers' % len(self.outmask))
      
      # Compute short cadence model, if available
      if self.has_sc:
        self.compute(sc = True)
        
      # Save
      self.finalize()
      self.plot_final()
      self.plot_page2()
      self.plot_info(self.dvs1)
      self.plot_info(self.dvs2)
      self.save_model()
      
      if self.make_fits:
        log.info('Generating FITS file...')
        MakeFITS(self)
        
    except:
    
      self.exception_handler(self.debug)
  
def Inject(ID, model = 'nPLD', t0 = None, per = None, dur = 0.1, depth = 0.001,
           mask = False, trn_win = 5, poly_order = 3, **kwargs):
  '''
  Run one of the :py:obj:`everest` models with injected transits and attempt to recover the
  transit depth at the end with a simple linear regression with a polynomial baseline. 
  The depth is stored in the :py:obj:`inject` attribute of the model (a dictionary) as 
  :py:obj:`rec_depth`. A control injection is also performed, in which the transits are injected 
  into the de-trended data; the recovered depth in the control run is stored in :py:obj:`inject` 
  as :py:obj:`rec_depth_control`. The results are plotted on page 2 of the data validation summary.
  
  :param int ID: The target id
  :param str model: The name of the :py:obj:`everest` model to run. Default `"PLD"`
  :param float t0: The transit ephemeris in days. Default is to draw from the uniform distributon [0., :py:obj:`per`)
  :param float per: The injected planet period in days. Default is to draw from the uniform distribution [2, 10]
  :param float dur: The transit duration in days. Must be in the range [0.05, 0.5]. Default 0.1
  :param float depth: The fractional transit depth. Default 0.001
  :param bool mask: Explicitly mask the in-transit cadences when computing the PLD model? Default :py:obj:`False`
  :param float trn_win: The size of the transit window in units of the transit duration
  :param int poly_order: The order of the polynomial used to fit the continuum
  
  '''
  
  # Randomize the planet params
  if per is None:
    a = 3.
    b = 10.
    per = a + (b - a) * np.random.random()
  if t0 is None:
    t0 = per * np.random.random()
  
  # Get the actual class
  _model = eval(model)
  inject = {'t0': t0, 'per': per, 'dur': dur, 'depth': depth, 'mask': mask,
            'poly_order': poly_order, 'trn_win': trn_win}
  
  # Define the injection class
  class Injection(_model):
    '''
    The :py:obj:`Injection` class is a special subclass of a user-selected :py:obj:`everest` model.
    See :py:func:`Inject` for more details.
    
    '''
    
    def __init__(self, *args, inject = None, parent_class = None, **kwargs):
      '''
      
      '''
      
      self.inject = inject
      self.parent_class = parent_class
      self.kwargs = kwargs
      super(Injection, self).__init__(*args, **kwargs)

    @property
    def name(self):
      '''
      
      '''
      
      if self.inject['mask']:
        maskchar = 'M'
      else:
        maskchar = 'U'
      return '%s_Inject_%s%g' % (self.parent_class, maskchar, self.inject['depth'])
    
    def load_tpf(self):
      '''
      Loads the target pixel files and injects transits at the pixel level.
      
      '''
      
      # Load the TPF
      super(Injection, self).load_tpf()
      log.info("Injecting transits...")
    
      # Inject the transits into the regular data
      transit_model = Transit(self.time, t0 = self.inject['t0'], per = self.inject['per'], dur = self.inject['dur'], depth = self.inject['depth'])
      for i in range(self.fpix.shape[1]):
        self.fpix[:,i] *= transit_model 
      self.fraw = np.sum(self.fpix, axis = 1)
      if self.inject['mask']:
        self.transitmask = np.array(list(set(np.concatenate([self.transitmask, np.where(transit_model < 1.)[0]]))), dtype = int)

      # Now inject into the short cadence, if available
      if self.has_sc:
        transit_model = Transit(self.sc_time, t0 = self.inject['t0'], per = self.inject['per'], dur = self.inject['dur'], depth = self.inject['depth'])
        for i in range(self.sc_fpix.shape[1]):
          self.sc_fpix[:,i] *= transit_model 
        self.sc_fraw = np.sum(self.sc_fpix, axis = 1)
        if self.inject['mask']:
          self.sc_transitmask = np.array(list(set(np.concatenate([self.sc_transitmask, np.where(transit_model < 1.)[0]]))), dtype = int)
 
    def recover_depth(self):
      '''
      Recovers the injected transit depth from the long cadence data with a simple LLS solver.
      The results are all stored in the :py:obj:`inject` attribute of the model.

      '''
      
      # Control run
      transit_model = Transit(self.time, t0 = self.inject['t0'], per = self.inject['per'], dur = self.inject['dur'], depth = self.inject['depth'])
      kwargs = dict(self.kwargs)
      kwargs.update({'clobber': False})
      control = eval(self.parent_class)(self.ID, is_parent = True, **kwargs)
      control.fraw *= transit_model 
      
      # Get params
      log.info("Recovering transit depth...")
      t0 = self.inject['t0']
      per = self.inject['per']
      dur = self.inject['dur']
      depth = self.inject['depth']
      trn_win = self.inject['trn_win']
      poly_order = self.inject['poly_order']
      
      for run, tag in zip([self, control], ['', '_control']):
      
        # Compute the model
        mask = np.array(list(set(np.concatenate([run.badmask, run.nanmask]))), dtype = int)
        flux = np.delete(run.flux / np.nanmedian(run.flux), mask)  
        time = np.delete(run.time, mask)
        transit_model = (Transit(time, t0 = t0, per = per, dur = dur, depth = depth) - 1) / depth
      
        # Count the transits
        t0 += np.ceil((time[0] - dur - t0) / per) * per
        ttimes0 = np.arange(t0, time[-1] + dur, per)
        tinds = []
        for tt in ttimes0:
          # Get indices for this chunk
          inds = np.where(np.abs(time - tt) < trn_win * dur / 2.)[0]
          # Ensure there's a transit in this chunk, and that
          # there are enough points for the polynomial fit
          if np.any(transit_model[inds] < 0.) and len(inds) > poly_order:
            tinds.append(inds)

        # Our design matrix
        sz = (poly_order + 1) * len(tinds)
        X = np.empty((0, 1 + sz), dtype = float)
        Y = np.array([], dtype = float)
        T = np.array([], dtype = float)
      
        # Loop over all transits
        for i, inds in enumerate(tinds):
          # Get the transit model
          trnvec = transit_model[inds].reshape(-1, 1)
          # Normalize the time array
          t = time[inds]
          t = (t - t[0]) / (t[-1] - t[0])
          # Cumulative arrays
          T = np.append(T, time[inds])
          Y = np.append(Y, flux[inds])
          # Polynomial vector
          polyvec = np.array([t ** o for o in range(0, poly_order + 1)]).T
          # Update the design matrix with this chunk
          lzeros = np.zeros((len(t), i * (poly_order + 1)))
          rzeros = np.zeros((len(t), sz - (i + 1) * (poly_order + 1)))
          chunk = np.hstack((trnvec, lzeros, polyvec, rzeros))
          X = np.vstack((X, chunk))

        # Get the relative depth
        A = np.dot(X.T, X)
        B = np.dot(X.T, Y)
        C = np.linalg.solve(A, B)
        rec_depth = C[0]
  
        # Get the uncertainties
        sig = 1.4826 * np.nanmedian(np.abs(flux - np.nanmedian(flux))) / np.nanmedian(flux)
        cov = sig ** 2 * np.linalg.solve(A, np.eye(A.shape[0]))
        err = np.sqrt(np.diag(cov))
        rec_depth_err = err[0]
      
        # Store the results
        self.inject.update({'rec_depth%s' % tag: rec_depth, 'rec_depth_err%s' % tag: rec_depth_err})
      
        # Store the detrended, folded data
        D = (Y - np.dot(C[1:], X[:,1:].T) + np.nanmedian(Y)) / np.nanmedian(Y)
        T = (T - t0 - per / 2.) % per - per / 2.  
        self.inject.update({'fold_time%s' % tag: T, 'fold_flux%s' % tag: D})

    def plot_page2(self):
      '''
      Plots the injection recovery results on page 2 of the DVS.
      
      '''
    
      # Plot the recovered folded transits
      ax1 = self.dvs2.lc1()      
      ax1.plot(self.inject['fold_time'], self.inject['fold_flux'], 'k.', alpha = 0.3)
      x = np.linspace(np.min(self.inject['fold_time']), np.max(self.inject['fold_time']), 500)
      try:
        y = Transit(x, t0 = 0., per = self.inject['per'], dur = self.inject['dur'], depth = self.inject['rec_depth'])
      except:
        # Log the error, and carry on
        exctype, value, tb = sys.exc_info()
        for line in traceback.format_exception_only(exctype, value):
          l = line.replace('\n', '')
          log.error(l)
        y = np.ones_like(x) * np.nan
      ax1.plot(x, y, 'r-')
      ax1.annotate('INJECTED', xy = (0.98, 0.025), xycoords = 'axes fraction', 
                   ha = 'right', va = 'bottom', fontsize = 10, alpha = 0.5, 
                   fontweight = 'bold')
      ax1.annotate('True depth:\nRecovered depth:',
                   xy = (0.02, 0.025), 
                   xycoords = 'axes fraction', 
                   ha = 'left', va = 'bottom', fontsize = 8, color = 'r')
      ax1.annotate('%.6f\n%.6f' % (self.inject['depth'], self.inject['rec_depth']),
                   xy = (0.25, 0.025), 
                   xycoords = 'axes fraction', 
                   ha = 'left', va = 'bottom', fontsize = 8, color = 'r')
      ax1.margins(0, None)
      ax1.ticklabel_format(useOffset=False)

      # Plot the recovered folded transits (control)
      ax2 = self.dvs2.lc2()      
      ax2.plot(self.inject['fold_time_control'], self.inject['fold_flux_control'], 'k.', alpha = 0.3)
      x = np.linspace(np.min(self.inject['fold_time_control']), np.max(self.inject['fold_time_control']), 500)
      try:
        y = Transit(x, t0 = 0., per = self.inject['per'], dur = self.inject['dur'], depth = self.inject['rec_depth_control'])
      except:
        # Log the error, and carry on
        exctype, value, tb = sys.exc_info()
        for line in traceback.format_exception_only(exctype, value):
          l = line.replace('\n', '')
          log.error(l)
        y = np.ones_like(x) * np.nan
      ax2.plot(x, y, 'r-')
      ax2.annotate('CONTROL', xy = (0.98, 0.025), xycoords = 'axes fraction', 
                   ha = 'right', va = 'bottom', fontsize = 10, alpha = 0.5, 
                   fontweight = 'bold')
      ax2.annotate('True depth:\nRecovered depth:',
                   xy = (0.02, 0.025), 
                   xycoords = 'axes fraction', 
                   ha = 'left', va = 'bottom', fontsize = 8, color = 'r')
      ax2.annotate('%.6f\n%.6f' % (self.inject['depth'], self.inject['rec_depth_control']),
                   xy = (0.25, 0.025), 
                   xycoords = 'axes fraction', 
                   ha = 'left', va = 'bottom', fontsize = 8, color = 'r')
      ax2.margins(0, None)
      ax2.ticklabel_format(useOffset=False)
      N = int(0.995 * len(self.inject['fold_flux_control']))
      hi, lo = self.inject['fold_flux_control'][np.argsort(self.inject['fold_flux_control'])][[N,-N]]
      fsort = self.inject['fold_flux_control'][np.argsort(self.inject['fold_flux_control'])]
      pad = (hi - lo) * 0.1
      ylim = (lo - pad, hi + pad)   
      ax2.set_ylim(ylim) 
      ax1.set_ylim(ylim)
    
      # Plot the PLD weights
      self.plot_weights()
    
    def finalize(self):
      '''
      Calls the depth recovery routine at the end of the de-trending step.
      
      '''
      
      super(Injection, self).finalize()
      self.recover_depth()
    
  return Injection(ID, inject = inject, parent_class = model, **kwargs)

def MaskedInject(ID, **kwargs):
  '''
  
  '''
  
  return Inject(ID, mask = True, **kwargs)

class rPLD(Model):
  '''
  The standard PLD model, inheriting all of its features from :py:class:`Model`.
  
  '''
        
  def __init__(self, *args, **kwargs):
    '''
    
    '''
    
    # Initialize
    super(rPLD, self).__init__(*args, **kwargs)
    
    # Check for saved model
    if self.load_model():
      return

    # Run
    self.run()
      
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

class nPLD(Model):
  '''
  The "neighboring stars" *PLD* model. This model uses the *PLD* vectors of neighboring
  stars to help in the de-trending and can lead to increased performance over the regular
  :py:class:`rPLD` model, particularly for dimmer stars.
  
  :param tuple cdpp_range: If :py:obj:`parent_model` is set, neighbors are selected only if \
                           their de-trended CDPPs fall within this range. Default `None`
  :param tuple mag_range: Only select neighbors whose magnitudes are within this range. \
                          Default (11., 13.) 
  :param int neighbors: The number of neighboring stars to use in the de-trending. The \
                        higher this number, the more signals there are and hence the more \
                        de-trending information there is. However, the neighboring star \
                        signals are regularized together with the target's signals, so adding \
                        too many neighbors will inevitably reduce the contribution of the \
                        target's own signals, which may reduce performance. Default `10`
  :param str parent_model: By default, :py:class:`nPLD` is run in stand-alone mode. The neighbor \
                           signals are computed directly from their TPFs, so there is no need to \
                           have run *PLD* on them beforehand. However, if :py:obj:`parent_model` is set, \
                           :py:class:`nPLD` will use information from the :py:obj:`parent_model` model of
                           each neighboring star when de-trending. This is particularly useful for \
                           identifying outliers in the neighbor signals and preventing them from polluting \
                           the current target. Setting :py:obj:`parent_model` to :py:class:`rPLD`, for instance, \
                           will use the outlier information in the :py:class:`rPLD` model of the neighbors \
                           (this must have been run ahead of time). Note, however, that tests with *K2* data \
                           show that including outliers in the neighbor signals actually *improves* the performance, \
                           since many of these outliers are associated with events such as thruster firings and are \
                           present in all light curves, and therefore *help* in the de-trending. Default `None`
  
  '''
        
  def __init__(self, *args, **kwargs):
    '''
    
    '''
    
    # Initialize
    super(nPLD, self).__init__(*args, **kwargs)
    
    # Check for saved model
    if self.load_model():
      return
    
    # Get neighbors
    self.parent_model = kwargs.get('parent_model', None)
    num_neighbors = kwargs.get('neighbors', 10)
    self.neighbors = GetNeighbors(self.ID, mission = self.mission, 
                                  model = self.parent_model,
                                  neighbors = num_neighbors, 
                                  mag_range = kwargs.get('mag_range', (11., 13.)), 
                                  cdpp_range = kwargs.get('cdpp_range', None),
                                  aperture_name = self.aperture_name)
    if len(self.neighbors):
      if len(self.neighbors) < num_neighbors:
        log.warn("%d neighbors requested, but only %d found." % (num_neighbors, len(self.neighbors)))
    else:
      log.error("No neighbors found! Aborting.")
      return
    
    for neighbor in self.neighbors:
      log.info("Loading data for neighboring target %d..." % neighbor)
      if self.parent_model is not None:
        # We load the `parent` model. The advantage here is that outliers have
        # properly been identified and masked
        data = eval(self.parent_model)(neighbor, mission = self.mission, is_parent = True)
      else:
        # We load the data straight from the TPF. Much quicker, since no model must
        # be run in advance. Downside is we don't know where the outliers are. But based
        # on tests with K2 data, the de-trending is actually *better* if the outliers are
        # included! These are mostly thruster fire events and other artifacts common to
        # all the stars, so it makes sense that we might want to keep them in the design
        # matrix.
        data = GetData(neighbor, self.mission, season = self.season, clobber = self.clobber_tpf, 
                       aperture_name = self.aperture_name, 
                       saturated_aperture_name = self.saturated_aperture_name, 
                       max_pixels = self.max_pixels,
                       saturation_tolerance = self.saturation_tolerance)
        data.mask = np.array(list(set(np.concatenate([data.badmask, data.nanmask]))), dtype = int)
        data.fraw = np.sum(data.fpix, axis = 1)
        if self.has_sc:
          data.sc_mask = np.array(list(set(np.concatenate([data.sc_badmask, data.sc_nanmask]))), dtype = int)
          data.sc_fraw = np.sum(data.sc_fpix, axis = 1)
      
      # Compute the linear PLD vectors and interpolate over outliers, NaNs and bad timestamps
      X1 = data.fpix / data.fraw.reshape(-1, 1)
      X1 = Interpolate(data.time, data.mask, X1)
      if self.has_sc:
        _scX1N = data.sc_fpix / data.sc_fraw.reshape(-1, 1)
        _scX1N = Interpolate(data.sc_time, data.sc_mask, _scX1N)
        if self.sc_X1N is None:
          self.sc_X1N = np.array(_scX1N)
        else:
          self.sc_X1N = np.hstack([self.sc_X1N, _scX1N])
        del _scX1N
      for n in range(self.pld_order):
        if self.XNeighbors[n] is None:
          self.XNeighbors[n] = X1 ** (n + 1)
        else:
          self.XNeighbors[n] = np.hstack([self.XNeighbors[n], X1 ** (n + 1)])
      del data
      del X1

    # Run
    self.run()
        
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
        self._X[n] = np.hstack([self._X[n], self.XNeighbors[n]])
  
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