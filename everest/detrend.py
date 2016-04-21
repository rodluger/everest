#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
detrend.py
----------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .utils import Mask, RMS, Chunks
import numpy as np
from scipy.misc import comb as nchoosek
from itertools import combinations_with_replacement as multichoose
from sklearn.decomposition import PCA
import logging
log = logging.getLogger(__name__)
      
def PLDBasis(fpix, time = None, breakpoints = None, pld_order = 1, cross_terms = True, 
             max_components = 300):
  '''
  Returns the basis vectors ``X`` we're going to use for the PLD model.
  ``X`` has shape (``npts``, ``n_components``). The first column ([:,0])
  is a vector of ones. The following columns are the principal components,
  in decreasing order of explained variance.
  
  '''
  
  # Our basis vector matrix: one set of basis vectors per chunk
  X = None
  
  # Fractional pixel fluxes
  npts, npix = fpix.shape
  frac = fpix / np.sum(fpix, axis = 1).reshape(-1, 1)
  
  # The number of signals (equation 6 in the paper)
  nsignals = int(np.sum([nchoosek(npix + k - 1, k) for k in range(1, pld_order + 1)]))
  
  # Add the breakpoints
  if breakpoints:
    assert time is not None, "Missing ``time`` array in ``PLDBasis()``."
    f = []
    m = 0
    for b in list(breakpoints):
      n = np.argmin(np.abs(time - b))
      f.append(frac[m:n])
      m = n
    f.append(frac[m:])
    frac = f
  else:
    frac = [frac]
  
  # Adjust the maximum number of PCA components; must be less than
  # the number of points in a chunk, otherwise they won't be
  # linearly independent! Must also, of course, be less than or
  # equal to the number of signals in the first place
  n_components = max_components
  for f in frac:
    n_components = min(n_components, f.shape[0], nsignals)
  
  # Setup the design matrix
  sz = n_components * len(frac)
  X = np.empty((0, sz), dtype = float)
  for i, f in enumerate(frac):
  
    # Get the PLD arrays
    if not cross_terms:
      x = np.hstack([f ** (n + 1) for n in range(pld_order)])
    else:
      x = np.empty(shape = (f.shape[0], 0), dtype = 'float64')
      for n in range(1, pld_order + 1):  
        xn = np.product(list(multichoose(f.T, n)), axis = 1).T
        x = np.hstack([x, xn])
        
    # Perform PCA on them
    pca = PCA(n_components = n_components - 1)
    xpca = pca.fit_transform(x)

    # Prepend a column vector of ones, since PCA transform removes 
    # the property that the basis vectors all sum to one.
    x = np.hstack([np.ones(xpca.shape[0]).reshape(-1, 1), xpca])
  
    # Pad with zeros on the left and right so that the chunks are
    # all independent of each other
    lzeros = np.zeros((x.shape[0], i * x.shape[1]))
    rzeros = np.zeros((x.shape[0], sz - (i + 1) * x.shape[1]))
    chunk = np.hstack((lzeros, x, rzeros))
    X = np.vstack((X, chunk))
  
  return X, n_components

def PLDModel(C, X):
  '''
  The PLD model. It's really simple!
  
  '''

  return np.dot(C, X.T)

def PLDCoeffs(X, Y, time, errors, gp, mask = []):
  '''
  Get the PLD coefficients by GLS regression with a GP.
  
  '''
  
  mask = Mask(mask)
  mT = mask(time)
  mE = mask(errors)
  mX = mask(X)
  mY = mask(Y)
  gp.compute(mT, mE)
  A = np.dot(mX.T, gp.solver.apply_inverse(mX))
  B = np.dot(mX.T, gp.solver.apply_inverse(mY))
  return np.linalg.solve(A, B)

def SliceX(X, n, npc):
  '''
  Reduce the dimensionality of X from ``npc`` components to ``n`` components.
  Trivial in the case where there's no breakpoints, slightly tricky when we have
  submatrices in ``X``.
  
  '''
  
  # Number of chunks in X
  nchunks = X.shape[1] // npc
  inds = np.concatenate([np.arange(i * npc + n, (i + 1) * npc) for i in range(nchunks)])
  return np.delete(X, inds, 1)

def ComputeScatter(X, Y, time, errors, gp, mask = [], niter = 30, nmasks = 10):
  '''
  Compute the median scatter in the de-trended light curve and the 
  median scatter in small masked portions of the light curve to
  gauge overfitting.
  
  '''
  
  # Setup some variables
  masked_scatter = []
  mask_orig = np.array(mask, dtype = int)
  mask_new = np.array(mask_orig)
  
  # Mask transits and outliers
  M = Mask(mask_orig)
  mT = M(time)
  mE = M(errors)
  mY = M(Y)
  mX = M(X)
  med = np.median(mY)
  
  # The precision in the unmasked light curve
  C = PLDCoeffs(X, Y, time, errors, gp, mask_orig)
  M = PLDModel(C, X)
  unmasked_scatter = RMS((Y - M + med) / med)

  # Get the indices of all contiguous 13-cadence segments of data,
  # separated into ``nmasks`` chunks
  sz = len(mT) // nmasks
  chunks = list(Chunks(mT, sz))
  inds = [[] for i in chunks]
  for c, chunk in enumerate(chunks):
    for i, t in enumerate(chunk[:-13]):
      if chunk[i + 13] - t <= 0.28:
        inds[c].append(i + c * sz)
  
  # Compute the precision several times and take the median
  for n in range(niter):
    
    # Get all our masks
    masks = [np.arange(s, s + 13) for s in [np.random.choice(i) for i in inds]]
    
    # Redefine the mask function
    mask_new = list(np.append(mask_orig, np.concatenate(masks)))
    
    # Get coefficients based on the masked data
    C = PLDCoeffs(X, Y, time, errors, gp, mask_new)
    
    # Predict the model in that interval
    M = [PLDModel(C, X[m]) for m in masks]
    
    # The masked de-trended interval and its precision in ppm
    masked_scatter.append(np.median([1.e6 * np.std((Y[m] - M[i] + med) / med) / np.sqrt(13) for i, m in enumerate(masks)]))
    
  # Take the median and return
  masked_scatter = np.median(masked_scatter)
  return masked_scatter, unmasked_scatter