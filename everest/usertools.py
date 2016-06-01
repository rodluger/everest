#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`usertools.py` - User tools
-----------------------------------

User tools to download, process, and plot the :py:class:`everest` light curves.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from . import __published__, __version__
from .config import EVEREST_DAT, EVEREST_SRC, MAST_ROOT, HYAK_ROOT
from .kernels import KernelModels
from .data import Campaign
from .utils import PlotBounds, MADOutliers
import k2plr as kplr
import george
import numpy as np
import matplotlib.pyplot as pl
import glob
import os, subprocess
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
  
  if __published__:
  
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
    
  else:

    # Get the url
    url = HYAK_ROOT + 'c%02d/' % campaign + ('%09d' % EPIC)[:4] + '00000/' + \
          'hlsp_everest_k2_llc_%d-c%02d_kepler_v%s.fits' % (EPIC, campaign, __version__)
  
    # Get the local file name
    filename = os.path.join(EVEREST_DAT, 'fits', 'c%02d' % campaign, 
                            ('%09d' % EPIC)[:4] + '00000',
                            'hlsp_everest_k2_llc_%d-c%02d_kepler_v%s.fits' % (EPIC, campaign, __version__))
    if not os.path.exists(os.path.dirname(filename)):
      os.makedirs(os.path.dirname(filename))
  
    # Download the data
    subprocess.call(['scp', url, filename])
  
  return filename

def _GetFITSFile(EPIC, clobber = False):
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

def _FlagColor(flag):
  '''
  
  '''
  
  if flag <= 1:
    return "g"
  elif flag <= 3:
    return "y"
  else:
    return "r"

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
  
  :returns: When called with an array `x` as an argument, a :py:class:`Mask` instance \
            returns the array `x` with the mask applied
              
  '''
  
  def __init__(self, ranges = [], transits = []):
    
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
    self.inds = []
    self.time = None
  
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

class Everest(object):
  '''
  
  '''
  
  def __init__(self, EPIC, clobber = False):
    '''
    
    '''
    
    self.EPIC = EPIC
    self.file = _GetFITSFile(self.EPIC, clobber = clobber)
    self.mask = Mask()
    self.detrended = False
    
  def set_mask(self, **kwargs):
    '''
    
    '''
    
    self.mask = Mask(**kwargs)
    self.detrended = False

  def detrend(self):
    '''
    Detrends with custom user options. If a local copy does not
    exist, automatically downloads the :py:mod:`everest` FITS file from MAST.
   
    '''
    
    # Check if we've done this already
    if self.detrended:
      return
    
    # Open the FITS file
    with pyfits.open(self.file) as hdulist:
    
      # Get the original arrays
      self.time = hdulist[1].data['TIME']
      # Hack to prevent george from complaining below
      self.time[np.isnan(self.time)] = 0 
      self.raw_flux = hdulist[1].data['RAW_FLUX']
      self.raw_ferr = hdulist[1].data['RAW_FERR']
    
      # Get the outliers
      outliers = hdulist[1].data['OUTLIER']
      oinds = np.where(outliers)[0]
    
      # Add the outliers to the mask object
      self.mask.inds = list(set(self.mask.inds + list(oinds)))
      self.mask.time = self.time
    
      # Get the gaussian process kernel
      knum = hdulist[3].header['KNUM']
      kpars = [hdulist[3].header['KPAR%02d' % n] for n in range(10)]
      kpars = [k for k in kpars if k != '']
      kernel = KernelModels[knum]
      kernel[:] = kpars
      self.gp = george.GP(kernel.george_kernel())
    
      # Get the design matrix
      self.X = hdulist[3].data['X']
    
      # Do some linear algebra to get the PLD coefficients `C`
      # and the PLD `model`
      self.gp.compute(self.mask(self.time), self.mask(self.raw_ferr))
      A = np.dot(self.mask(self.X).T, self.gp.solver.apply_inverse(self.mask(self.X)))
      B = np.dot(self.mask(self.X).T, self.gp.solver.apply_inverse(self.mask(self.raw_flux)))
      self.C = np.linalg.solve(A, B)
      self.model = np.dot(self.C, self.X.T)
    
      # Subtract the model and add the median back in to get
      # our final de-trended flux
      self.flux = self.raw_flux - self.model + np.nanmedian(self.raw_flux)
      
      # Warning flags
      self.crwdflag = hdulist[1].header['CRWDFLAG']  
      self.satflag = hdulist[1].header['SATFLAG']
      
      # Set flag
      self.detrended = True

  def plot(self, pipeline = 'everest'):
    '''
    Plot the raw and de-trended light curves.
  
    '''
    
    if not self.detrended:
      self.detrend()
    
    fig, ax = pl.subplots(2, figsize = (12,8), sharex = True)
    fig.subplots_adjust(left = 0.1, right = 0.95, 
                        top = 0.925, bottom = 0.1,
                        hspace = 0.05)
    # Plot                   
    ax[0].plot(self.mask(self.time), self.mask(self.raw_flux), 'k.', markersize = 3, alpha = 0.5)
    ax[0].plot(self.mask.inv(self.time), self.mask.inv(self.raw_flux), 'r.', markersize = 3, alpha = 0.5,
               label = 'Masked')
    ax[0].legend(loc = 'upper left', fontsize = 9, numpoints = 3)
    
    if pipeline.lower() == 'everest':
      ax[1].plot(self.mask(self.time), self.mask(self.flux), 'b.', markersize = 3, alpha = 0.5)
      ax[1].plot(self.mask.inv(self.time), self.mask.inv(self.flux), 'r.', markersize = 3, alpha = 0.5)
      ax[1].set_ylabel('EVEREST Flux', fontsize = 18)
      bounds = min(PlotBounds(self.raw_flux), PlotBounds(self.flux))
    elif pipeline.lower() == 'k2sff':
      try:
        k2sff = kplr.K2SFF(self.EPIC)
      except:
        return
      k2sff.fcor *= np.nanmedian(self.flux)
      ax[1].plot(k2sff.time, k2sff.fcor, 'b.', markersize = 3, alpha = 0.5)
      ax[1].set_ylabel('K2SFF Flux', fontsize = 18)
      bounds = min(PlotBounds(self.raw_flux), PlotBounds(k2sff.fcor))      
    elif pipeline.lower() == 'k2sff':
      try:
        k2sc = kplr.K2SC(self.EPIC)
      except:
        return
      ax[1].plot(k2sc.time, k2sc.pdcflux, 'b.', markersize = 3, alpha = 0.5)
      ax[1].set_ylabel('K2SC Flux', fontsize = 18)
      bounds = min(PlotBounds(self.raw_flux), PlotBounds(k2sc.pdcflux))
    elif pipeline.lower() == 'k2varcat':
      try:
        k2varcat = kplr.K2VARCAT(self.EPIC)
      except:
        return
      k2varcat.flux *= np.nanmedian(self.flux)
      ax[1].plot(k2varcat.time, k2varcat.flux, 'b.', markersize = 3, alpha = 0.5)
      ax[1].set_ylabel('K2VARCAT Flux', fontsize = 18)
      bounds = min(PlotBounds(self.raw_flux), PlotBounds(k2varcat.flux))    
    else:
      return
      
    # Labels
    ax[0].set_ylabel('SAP Flux', fontsize = 18)
    ax[1].set_xlabel('Time (BJD - 2454833)', fontsize = 18)
    pl.suptitle('EPIC %d' % self.EPIC, fontsize = 22)
    fig.canvas.set_window_title('EPIC %d' % self.EPIC)
  
    # Set the appropriate bounds (same for both plots)
    ax[0].set_ylim(*bounds)
    ax[1].set_ylim(*bounds)
  
    # Warning flags
    ax[1].annotate("Crowding:   %d/5" % self.crwdflag, xy = (0.02, 0.95), xycoords = "axes fraction", ha="left", va="top", fontsize=12, color=_FlagColor(self.crwdflag))
    ax[1].annotate("Saturation: %d/5" % self.satflag, xy = (0.02, 0.885), xycoords = "axes fraction", ha="left", va="top", fontsize=12, color=_FlagColor(self.satflag))
  
    return fig, ax

  def plot_folded(self):
    '''
    Plot the raw and de-trended folded light curves.
  
    '''
    
    if not len(self.mask._transits):
      return None, None
    
    if not self.detrended:
      self.detrend()
  
    # Compute whitened flux
    med = np.nanmedian(self.mask(self.flux))
    outliers = MADOutliers(self.mask(self.time), self.mask(self.flux))  
    self.gp.compute(np.delete(self.mask(self.time), outliers), np.delete(self.mask(self.raw_ferr), outliers))
    mu, _ = self.gp.predict(np.delete(self.mask(self.flux), outliers) - med, self.time)
    fwhite = (self.flux - mu) / med
     
    # Instantiate
    fig, ax = pl.subplots(len(self.mask._transits), figsize = (6, 5 + 1.5 * len(self.mask._transits) - 1), sharex = True)
    ax = np.atleast_1d(ax)
    fig.subplots_adjust(left = 0.2, right = 0.95, top = 0.95, bottom = 0.1)

    # Calculate time indices from the transits
    for n, transit in enumerate(self.mask._transits):
      if np.all([not hasattr(r, '__len__') for r in transit]) and len(transit) == 3:
        # This is a 3-tuple for a single planet
        per, t0, dur = transit
        t0 += np.ceil((self.time[0] - dur - t0) / per) * per
        fold = lambda t: (t - t0 - per / 2.) % per - per / 2.
        ax[n].plot(fold(self.time), fwhite, 'b.', alpha = 0.4)
        ax[n].set_xlim(-dur * 0.75, dur * 0.75)

      else:
        # This is a list of 2-tuples corresponding to transit times
        ftime = np.array(self.time)
        for i, t in enumerate(self.time):
          ftime[i] -= transit[np.argmin([np.abs(t - tup[0]) for tup in transit])][0]
        ax[n].plot(ftime, fwhite, 'b.', alpha = 0.4)
        ax[n].set_xlim(-transit[0][1] * 0.75, transit[0][1] * 0.75)
    
      # Labels
      ax[n].set_title('EPIC %d.%02d' % (self.EPIC, n + 1))
      ax[n].set_ylabel('Norm. Flux', fontsize = 14)
  
    ax[-1].set_xlabel('Time (days)', fontsize = 18)
    fig.canvas.set_window_title('EPIC %d' % self.EPIC)
  
    return fig, ax