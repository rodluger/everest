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
  
def DetrendWithMask(EPIC):
  '''
  
  '''
  
  file = _FITSFile(EPIC)
  with pyfits.open(file) as hdulist:
    
    # Get the original arrays
    time = hdulist[1].data['TIME']
    flux = hdulist[1].data['RAW_FLUX']
    
    # Get the outlier mask
    outliers = hdulist[1].data['OUTLIER']
    M = Mask(np.where(outliers))
    
    # Get the gaussian process kernel
    knum = hdulist[3].header['KNUM']
    kpars = [hdulist[3].header['KPAR%02d' % n] for n in range(10)]
    
    import pdb; pdb.set_trace()
    
    kernel = KernelModels[knum]
    kernel[:] = kpars
    gp = george.GP(kernel.george_kernel())
    
    # Get the design matrix
    X = hdulist[3].data['X']
    
    gp.compute(M(time), M(ferr))
    A = np.dot(M(X).T, gp.solver.apply_inverse(M(X)))
    B = np.dot(M(X).T, gp.solver.apply_inverse(M(flux)))
    C = np.linalg.solve(A, B)
    model = np.dot(C, X.T)
    fpld = flux - model
    fpld += np.median(flux)
    fpld_norm = fpld / np.median(M(fpld))
