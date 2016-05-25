#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`usertools.py` - User tools
-----------------------------------

User tools to download, process, and plot the :py:class:`everest` light curves.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .config import EVEREST_DAT, EVEREST_SRC, MAST_ROOT
from .kernels import KernelModels
from .data import Campaign
from .utils import PlotBounds, MADOutliers
import george
import numpy as np
import matplotlib.pyplot as pl
import glob
import os
try:
  import pyfits
except ImportError:
  try:
    import astropy.io.fits as pyfits
  except ImportError:
    raise Exception('Please install the `pyfits` package.')
from tempfile import NamedTemporaryFile
import six
from six.moves import urllib
import shutil

def _QtFigureLayout(n, nsub):
  '''
  Adapted from http://stackoverflow.com/a/29039755

  '''
  
  # Figure manager
  mgr = pl.get_current_fig_manager()
  
  # Toggle full screen to get size
  mgr.full_screen_toggle()
  py = mgr.canvas.height()
  px = mgr.canvas.width()
  mgr.full_screen_toggle()

  # Width of the window border in pixels
  d = 10  
  if n == 0:
    # x-top-left-corner, y-top-left-corner, x-width, y-width (in pixels)
    mgr.window.setGeometry(d, 2 * d, 2 * px / 3 - d, py - 2 * d)
  else:

    
    height = py / nsub - 4 * d
    top = 2 * d + (n - 1) * height + 3 * d * n
    
    mgr.window.setGeometry(2 * px / 3 + d, top, px / 3 - 2 * d, height)

  # Bring to front
  mgr.window.raise_()
  
def _EverestVersion():
  '''
  Returns the current :py:mod:`everest` version on MAST.
  
  '''
  
  url = MAST_ROOT + 'version.txt'
  r = urllib.request.Request(url)
  handler = urllib.request.urlopen(r)
  code = handler.getcode()
  if int(code) != 200:
    raise Exception("Error code {0} for URL '{1}'".format(code, url))
  data = handler.read().decode('utf-8')
  if data.endswith('\n'):
    data = data[:-1]
  return data

def _DownloadFITSFile(EPIC, clobber = False):
  '''
  Download a given :py:mod:`everest` FITS file from MAST.
  Returns the full local path to the FITS file.
  
  :param bool clobber: If `True`, download and overwrite existing files. Default `False`
  
  '''
    
  # Get the campaign
  campaign = Campaign(EPIC)  
    
  # Get the url
  mast_version = _EverestVersion()
  url = MAST_ROOT + 'c%02d/' % campaign + ('%09d' % EPIC)[:4] + '00000/' + \
        'hlsp_everest_k2_llc_%d-c%02d_kepler_v%s.fits' % (EPIC, campaign, mast_version)
  
  # Get the local file name
  filename = os.path.join(EVEREST_DAT, 'fits', 'c%02d' % campaign, 
                          ('%09d' % EPIC)[:4] + '00000',
                          'hlsp_everest_k2_llc_%d-c%02d_kepler_v%s.fits' % (EPIC, campaign, mast_version))
  if not os.path.exists(os.path.dirname(filename)):
    os.makedirs(os.path.dirname(filename))
  
  # Download the data
  r = urllib.request.Request(url)
  handler = urllib.request.urlopen(r)
  code = handler.getcode()
  if int(code) != 200:
    raise Exception("Error code {0} for URL '{1}'".format(code, url))
  data = handler.read()

  # Atomically save to disk
  f = NamedTemporaryFile("wb", delete=False)
  f.write(data)
  f.flush()
  os.fsync(f.fileno())
  f.close()
  shutil.move(f.name, filename)
  
  return filename
  
def GetFITSFile(EPIC, clobber = False):
  '''
  Returns the path to a given :py:mod:`everest` FITS file.
  In case of multiple versions, returns the file corresponding to the
  latest :py:mod:`everest` version. In case a local copy does not exist,
  downloads the most recent version from MAST.
  
  :param bool clobber: If `True`, download and overwrite existing files. Default `False`
  
  '''
  
  campaign = Campaign(EPIC)
  path = os.path.join(EVEREST_DAT, 'fits', 'c%02d' % campaign, 
                     ('%09d' % EPIC)[:4] + '00000',
                      'hlsp_everest_k2_llc_%d-c%02d_kepler_v*.fits' % (EPIC, campaign))
  files = glob.glob(path)
  if len(files):
    return files[-1]
  else:
    return _DownloadFITSFile(EPIC, clobber = clobber)

class Mask(object):
  '''
  An object containing information about which portions of the light curve to mask
  when computing the PLD coefficients.
  
  :param list ranges: A list of 2-tuples in the form `(low, high)`, where \
                      `low` and `high` are the lower and upper bounds, respectively, \
                      of the regions to be masked. These should be in units of \
                      `BJD - 245833` (same as the `time` array). Default `[]`
  :param list transits: A list of 3-tuples in the form `(per, t0, dur)` and/or \
                        lists of 2-tuples in the form `(t, dur)` for any \
                        transits/eclipses to be masked. In the first case, \
                        the parameter `per` is the \
                        transit period in days, `t0` is the time of first transit,
                        and `dur` is the total transit duration in days. (For eclipsing \
                        binaries with secondary eclipses, specify two \
                        separate tuples, one for primary eclipse and one for secondary \
                        eclipse.) The second case is for non-periodic transits and allows \
                        the user to specify the time `t` of *each* transit and its duration \
                        `dur`.
  
  :param dict kwargs: Extra keyword arguments for this class (for internal use only)                
  
  :returns: When called with an array `x` as an argument, a :py:class:`Mask` instance \
            returns the array `x` with the mask applied
              
  '''
  
  def __init__(self, ranges = [], transits = [], **kwargs):
    
    # Check the ranges param
    if len(ranges):
      if np.all([hasattr(r, '__len__') for r in ranges]):
        if np.any([len(r) != 2 for r in ranges]):
          raise Exception('Param `ranges` must be a list of `(low, high)` pairs.')
      elif np.all([not hasattr(r, '__len__') for r in ranges]):
        ranges = [ranges]
      else:
        raise Exception('Param `ranges` must be a list of `(low, high)` pairs.')
    self._ranges = ranges
    
    # Check the transits param
    if len(transits):
      if np.all([hasattr(r, '__len__') for r in transits]):
        if np.all([len(r) == 2 for r in transits]):
          # User provided a single planet with ttvs
          transits = [transits]
        else:
          for transit in transits:
            if np.any([hasattr(r, '__len__') for r in ranges]) and not np.all([hasattr(r, '__len__') for r in ranges]):
              # User provided some scalars and some lists
              raise Exception('Invalid value for the `transits` param. Please consult the docs.')
            elif np.all([hasattr(r, '__len__') for r in ranges]):
              if np.any([len(r) != 2 for r in ranges]):
                # User provided all lists, but they're not the right length
                raise Exception('Invalid value for the `transits` param. Please consult the docs.')
      elif np.all([not hasattr(r, '__len__') for r in transits]):
        if len(transits) == 3 or len(transits) == 2:
          # User provided a single planet or a single transit time
          transits = [transits]
        else:
          raise Exception('Invalid value for the `transits` param. Please consult the docs.')
      else:
        raise Exception('Invalid value for the `transits` param. Please consult the docs.')
    self._transits = transits
    
    # Internal
    self.inds = kwargs.get('inds', [])
    self.time = kwargs.get('time', None)
  
  @property
  def all_inds(self):
    '''
    Returns the indices of all points that will be masked.
    
    '''
    
    # Check that the user provided `time`
    if self.time is None:
      raise Exception("Please set the `time` property!")
    
    # Calculate time indices from the ranges
    rinds = []
    for lo, hi in self._ranges:
      rinds.extend(np.where((self.time >= lo) & (self.time <= hi))[0])
    
    # Calculate time indices from the transits
    tinds = []
    for transit in self._transits:
      if np.all([not hasattr(r, '__len__') for r in transit]) and len(transit) == 3:
        # This is a 3-tuple for a single planet
        per, t0, dur = transit
        t0 += np.ceil((self.time[0] - dur - t0) / per) * per
        for t in np.arange(t0, self.time[-1] + dur, per):
          tinds.extend(np.where(np.abs(self.time - t) < dur / 2.)[0])
      else:
        # This is a list of 2-tuples corresponding to transit times
        for tup in transit:
          # This is a 2-tuple for a single transit 
          t, dur = tup
          tinds.extend(np.where(np.abs(self.time - t) < dur / 2.)[0])
        
    # Add in the explicit indices
    inds = list(set(self.inds + list(rinds) + list(tinds)))
    
    return inds
  
  def inv(self, x):
    '''
    Returns the actual masked portion of the array `x`;
    this is the inverse of what a call to a `Mask` object returns.
    
    '''
    
    return x[self.all_inds]
  
  def __call__(self, x):
    '''
    Returns the masked version of `x`.
    
    '''
    
    return np.delete(x, self.all_inds, axis = 0)

def _Plot(EPIC, time, flux, fpld, fwhite, mask):
  '''
  Plot the raw and de-trended light curves.
  
  '''
  
  fig, ax = pl.subplots(2, figsize = (12,8), sharex = True)
  fig.subplots_adjust(left = 0.1, right = 0.95, 
                      top = 0.925, bottom = 0.1,
                      hspace = 0.05)
  # Plot                   
  ax[0].plot(mask(time), mask(flux), 'k.', markersize = 3, alpha = 0.5)
  ax[0].plot(mask.inv(time), mask.inv(flux), 'r.', markersize = 3, alpha = 0.5,
             label = 'Masked')
  ax[0].legend(loc = 'upper left', fontsize = 9, numpoints = 3)
  ax[1].plot(mask(time), mask(fpld), 'b.', markersize = 3, alpha = 0.5)
  ax[1].plot(mask.inv(time), mask.inv(fpld), 'r.', markersize = 3, alpha = 0.5)

  # Labels
  ax[0].set_ylabel('SAP Flux', fontsize = 18)
  ax[1].set_ylabel('EVEREST Flux', fontsize = 18)
  ax[1].set_xlabel('Time (BJD - 2454833)', fontsize = 18)
  pl.suptitle('EPIC %d' % EPIC, fontsize = 22)
  fig.canvas.set_window_title('EPIC %d' % EPIC)
  
  # Set the appropriate bounds (same for both plots)
  bounds = min(PlotBounds(flux), PlotBounds(fpld))
  ax[0].set_ylim(*bounds)
  ax[1].set_ylim(*bounds)
  try:
    _QtFigureLayout(0, len(mask._transits))
  except:
    pass
  
  # Plot folded transits?
  if fwhite is not None:
    figs = _PlotFolded(EPIC, time, flux, fpld, fwhite, mask)
  
  pl.show()

def _PlotFolded(EPIC, time, flux, fpld, fwhite, mask):
  '''
  Plot the raw and de-trended folded light curves.
  
  '''
  
  figs = []
  
  # Calculate time indices from the transits
  for n, transit in enumerate(mask._transits):
    
    # Instantiate
    fig, ax = pl.subplots(1, figsize = (4,6), sharex = True)
    fig.subplots_adjust(left = 0.15, right = 0.95, 
                      top = 0.95, bottom = 0.1)
    
    if np.all([not hasattr(r, '__len__') for r in transit]) and len(transit) == 3:
      # This is a 3-tuple for a single planet
      per, t0, dur = transit
      t0 += np.ceil((time[0] - dur - t0) / per) * per
      fold = lambda t: (t - t0 - per / 2.) % per - per / 2.
      ax.plot(fold(time), fwhite, 'b.', alpha = 0.4)
      ax.set_xlim(-dur * 0.75, dur * 0.75)

    else:
      # This is a list of 2-tuples corresponding to transit times
      ftime = np.array(time)
      for i, t in enumerate(time):
        ftime[i] -= transit[np.argmin([np.abs(t - tup[0]) for tup in transit])][0]
      ax.plot(ftime, fwhite, 'b.', alpha = 0.4)
      ax.set_xlim(-transit[0][1] * 0.75, transit[0][1] * 0.75)
    
    # Labels
    fig.canvas.set_window_title('EPIC %d.%02d' % (EPIC, n + 1))
    
    # Layout on screen
    fmg = pl.get_current_fig_manager()
    try:
      _QtFigureLayout(n + 1, len(mask._transits))
    except:
      pass
    
    figs.append(fig)
  
  return figs
  
def Detrend(EPIC, mask = None, clobber = False, plot = False):
  '''
  Detrends a given EPIC target with custom user options. If a local copy does not
  exist, automatically downloads the :py:mod:`everest` FITS file from MAST.
  
  :param everest.tools.Mask mask: A :py:class:`Mask` instance containing information \
                                  on which portions of the light curve to mask.
  :param bool clobber: If `True`, download and overwrite existing files. Default `False`
  :param bool plot: Plot the light curve? Default `False`
  
  :returns: `(time, flux, mask_inds)`, the time and de-trended flux arrays, as well as the \
            indices corresponding to the points that were masked
   
  '''
  
  file = GetFITSFile(EPIC, clobber = clobber)
  with pyfits.open(file) as hdulist:
    
    # Get the original arrays
    time = hdulist[1].data['TIME']
    flux = hdulist[1].data['RAW_FLUX']
    ferr = hdulist[1].data['RAW_FERR']
    
    # Get the outliers
    outliers = hdulist[1].data['OUTLIER']
    oinds = np.where(outliers)[0]
    
    # Add the outliers to the mask object
    if mask is None:
      mask = Mask()
    mask.inds = list(set(mask.inds + list(oinds)))
    mask.time = time

    # Get the gaussian process kernel
    knum = hdulist[3].header['KNUM']
    kpars = [hdulist[3].header['KPAR%02d' % n] for n in range(10)]
    kpars = [k for k in kpars if k != '']
    kernel = KernelModels[knum]
    kernel[:] = kpars
    gp = george.GP(kernel.george_kernel())
    
    # Get the design matrix
    X = hdulist[3].data['X']
    
    # Do some linear algebra to get the PLD coefficients `C`
    # and the PLD `model`
    gp.compute(mask(time), mask(ferr))
    A = np.dot(mask(X).T, gp.solver.apply_inverse(mask(X)))
    B = np.dot(mask(X).T, gp.solver.apply_inverse(mask(flux)))
    C = np.linalg.solve(A, B)
    model = np.dot(C, X.T)
    
    # Subtract the model and add the median back in to get
    # our final de-trended flux
    fpld = flux - model + np.median(flux)
        
    # Plot?
    if plot:
      # Compute the fully whitened flux (only if we're plotting folded transits)
      if len(mask._transits):
        fwhite = flux - model
        fwhite += np.median(flux)
        med = np.median(mask(fwhite))
        outliers = MADOutliers(mask(time), mask(fwhite))  
        gp.compute(np.delete(mask(time), outliers), np.delete(mask(ferr), outliers))
        mu, _ = gp.predict(np.delete(mask(fwhite), outliers) - med, time)
        fwhite = (fwhite - mu) / med
      else:
        fwhite = None
        frw = None
      _Plot(EPIC, time, flux, fpld, fwhite, mask)
    
    return time, fpld, mask.all_inds