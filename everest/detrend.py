#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
detrend.py
----------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .utils import Mask, RMS, Chunks
import numpy as np
import george
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

def Outliers(time, flux, fpix, ferr, mask = [], sigma = 5):
  '''
  Return the indices of outliers we should remove from the light curve when computing
  the PLD coeffs.
  
  '''

  # Mask the arrays right off the bat
  time_orig = np.array(time)
  time = np.delete(time, mask)
  flux = np.delete(flux, mask)
  fpix = np.delete(fpix, mask, axis = 0)
  ferr = np.delete(ferr, mask)

  # Set up a generic GP
  amp = np.median([np.std(y) for y in Chunks(flux, int(2. / np.median(time[1:] - time [:-1])))])
  gp = george.GP(amp ** 2 * george.kernels.Matern32Kernel(2. ** 2))
  
  # Compute the basis vectors for 1st order PLD w/ 5 chunks
  nchunks = 5
  brkpts = [time[int(k)] for k in np.linspace(0, len(time), nchunks + 1)[1:-1]]
  X, _ = PLDBasis(fpix, time = time, pld_order = 1, max_components = 50, breakpoints = brkpts)
  
  # First we (tentatively) clip outliers from the raw flux.
  med = np.median(flux)
  MAD = 1.4826 * np.median(np.abs(flux - med))
  i = np.where((flux > med + sigma * MAD) | (flux < med - sigma * MAD))[0]
  log.info('Iteration #00: %d outliers.' % len(i))
  
  # Now do iterative sigma clipping. 
  j = [-1]
  count = 0
  while not np.array_equal(i, j):
    
    # Reset; the loop ends when we get the same outliers twice in a row
    j = i
    count += 1
    if count > 25:
      # We will continue, even though there may be issues
      log.error('Maximum number of iterations in ``Outliers()`` exceeded.')
      break
    
    # Remove both the PLD component and the GP component
    C = PLDCoeffs(X, flux, time, ferr, gp, mask = i)
    M = PLDModel(C, X)
    mu, _ = gp.predict(np.delete(flux - M, i), time)
    fdet = flux - M - mu
  
    # Clip!
    med = np.median(fdet)
    MAD = 1.4826 * np.median(np.abs(fdet - med))
    i = np.where((fdet > med + sigma * MAD) | (fdet < med - sigma * MAD))[0]
  
    # Log
    log.info('Iteration #%02d: %d outliers.' % (count, len(i)))
  
  return [np.argmax(time_orig == time[j]) for j in i]