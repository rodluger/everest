#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`tools.py` - User tools
-------------------------------

User tools.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os
EVEREST_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
from .kernels import KernelModels
import george
import numpy as np
import matplotlib.pyplot as pl
import glob
try:
  import pyfits
except ImportError:
  try:
    import astropy.io.fits as pyfits
  except ImportError:
    raise Exception('Please install the `pyfits` package.')

def _DownloadFITSFile(EPIC):
  '''
  Download a given :py:mod:`everest` FITS file from MAST.
  Returns the full local path to the FITS file.
  
  '''
  
  raise NotImplementedError('Downloading functionality coming soon!')

def _FITSFile(EPIC):
  '''
  Returns the path to a given :py:mod:`everest` FITS file.
  In case of multiple versions, returns the file corresponding to the
  latest :py:mod:`everest` version.
  
  '''
  
  path = os.path.join(EVEREST_ROOT, 'fits', 'c*', 
                     ('%09d' % EPIC)[:4] + '00000',
                      'hlsp_everest_k2_llc_%d-c*_kepler_v*.fits' % EPIC)
  files = glob.glob(path)
  if len(files):
    return files[-1]
  else:
    return _DownloadFITSFile(EPIC)

class Mask(object):
  '''
  
  '''
  
  def __init__(self, ranges = [], axis = 0, inds = [], time = None):
    self._axis = axis
    self._ranges = ranges
    self.inds = inds
    self.time = time
  
  def __call__(self, x):
    
    # Check that the user provided `time`
    if self.time is None:
      raise Exception("Please set the `time` property!")
    
    # Calculate time indices
    inds = []
    for lo, hi in self._ranges:
      inds.extend(np.where((self.time >= lo) & (self.time <= hi))[0])
    
    # Add in the explicit indices
    inds = list(set(self.inds + list(inds)))
  
    # Return the masked array
    return np.delete(x, inds, axis = self._axis)
  
def Detrend(EPIC, mask = None):
  '''
  
  '''
  
  file = _FITSFile(EPIC)
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

    pl.plot(time, fpld)
    pl.show()