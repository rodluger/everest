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
  
  return data

def _DownloadFITSFile(EPIC, clobber = False):
  '''
  Download a given :py:mod:`everest` FITS file from MAST.
  Returns the full local path to the FITS file.
  
  '''
    
  # Get the url
  mast_version = _EverestVersion()
  url = MAST_ROOT + 'c%02d' % campaign + ('%09d' % EPIC)[:4] + '00000' + \
        'hlsp_everest_k2_llc_%d-c%02d_kepler_v%s.fits' % (EPIC, campaign, mast_version)
  
  # Get the local file name
  campaign = Campaign(EPIC)
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
  :param list transits: Either a list of 3-tuples in the form `(per, t0, dur)` or a \
                        list of 2-tuples in the form `(t, dur)` for any \
                        transits/eclipses to be masked. In the first case, \
                        the parameter `per` is the \
                        transit period in days, `t0` is the time of first transit,
                        and `dur` is the total transit duration in days. (For eclipsing \
                        binaries with secondary eclipses, consider specifying two \
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
        if np.any([len(r) != 2 for r in transits]) and np.any([len(r) != 3 for r in transits]):
          raise Exception('Invalid value for the `transits` param. Please consult the docs.')
      elif np.all([not hasattr(r, '__len__') for r in transits]):
        transits = [transits]
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
      if len(transit) == 3:
        per, t0, dur = transit
        t0 += np.ceil((self.time[0] - dur - t0) / per) * per
        for t in np.arange(t0, self.time[-1] + dur, per):
          tinds.extend(np.where(np.abs(self.time - t) < dur / 2.)[0])
      elif len(transit) == 2:
        t, dur = transit
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
  
def Detrend(EPIC, mask = None, clobber = False):
  '''
  Detrends a given EPIC target with custom user options. If a local copy does not
  exist, automatically downloads the :py:mod:`everest` FITS file from MAST.
  
  :param everest.tools.Mask mask: A :py:class:`Mask` instance containing information \
                                  on which portions of the light curve to mask.
  
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
    
    return time, fpld, mask.all_inds