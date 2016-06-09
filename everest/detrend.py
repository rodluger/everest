#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`detrend.py` - Linear algebra
-------------------------------------

This is the heart of :py:mod:`everest`. These routines compute the principal components
of the fractional pixel flux functions and solve the GLS problem with a GP to
obtain the PLD coefficients that best capture the instrumental noise signal.

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
      
def PLDBasis(fpix, time = None, breakpoints = None, pld_order = 3, cross_terms = True, 
             do_pca = True, max_components = 300):
  '''
  Returns the basis vectors :math:`\mathbf{X}` we're going to use for the PLD model.
  :math:`\mathbf{X}` has shape (`npts`, `n_components`). The first column (`[:,0]`)
  is a vector of ones. The following columns are the principal components,
  in decreasing order of explained variance. This is essentially Equation (6) in the paper.
  
  ::
      
                  [x00  x01  ...  x0N    1]
       X  =       [x10  x11  ...  x1N    1] 
                  [...  ...  ...  ...  ...]
                  [xM0  xM1  ...  xMN    1]
  
  :param ndarray fpix: The pixel flux array, shape (`npts`, `npixels`)
  :param ndarray time: The time array. Only necessary if `breakpoints` are specified. Default `None`
  :param list breakpoints: The time(s) at which to split the light curve. Default `None`
  :param int pld_order: The order of PLD to use. Default `3`
  :param bool cross_terms: If `True`, includes the cross terms :math:`\prod_{k != j} p_{ij}p_{ik}`. Default `True`
  :param int max_components: Maximum number of PCA components to produce. Default `300`
  
  :returns (X, n_components): The design matrix :math:`\mathbf{X}` and the number of \
                              regressors in each submatrix (if no breakpoints are \
                              specified, this is just the number of columns in :math:`\mathbf{X}`)
  
  '''
  
  # Our basis vector matrix: one set of basis vectors per chunk
  X = None
  
  # Fractional pixel fluxes
  npts, npix = fpix.shape
  frac = fpix / np.sum(fpix, axis = 1).reshape(-1, 1)
  
  # The number of signals
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
    if do_pca:
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
  
  :param ndarray C: The coefficients returned by :py:func:`PLDCoeffs`
  :param ndarray X: The design matrix returned by :py:func:`PLDBasis`
  
  :returns: :math:`\mathbf{X\cdot C}`
  '''

  return np.dot(C, X.T)

def PLDCoeffs(X, Y, time, errors, gp, mask = []):
  '''
  Get the PLD coefficients by GLS regression with a GP.
  
  :param ndarray X: The design matrix returned by :py:func:`PLDBasis`
  :param ndarray Y: The dependent variable array (the SAP flux)
  :param ndarray time: The independent variable array (the timestamps)
  :param ndarray errors: The standard errors on `Y`
  :param george.GP gp: The pre-initialized :py:mod:`george` gaussian process object
  :param list mask: The indices of points in `Y` to mask when computing the coefficients
  
  :returns C: The PLD coefficient array
  
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
  Reduce the dimensionality of :math:`\mathbf{X}` from `npc` components to `n` components.
  Trivial in the case where there's no breakpoints, slightly tricky when we have
  submatrices in :math:`\mathbf{X}`.
  
  :param ndarray X: The design matrix returned by :py:func:`PLDBasis`
  :param int n: The desired number of components (per submatrix)
  :param int npc: The original number of components (per submatrix)
  
  '''
  
  # Number of chunks in X
  nchunks = X.shape[1] // npc
  inds = np.concatenate([np.arange(i * npc + n, (i + 1) * npc) for i in range(nchunks)])
  return np.delete(X, inds, 1)

def ComputeScatter(X, Y, time, errors, gp, mask = [], niter = 30, nmasks = 10):
  '''
  Compute the median CDPP in the de-trended light curve and the 
  median CDPP in small masked portions of the light curve to
  gauge overfitting. This is the cross-validation step.
  
  :param ndarray X: The design matrix returned by :py:func:`PLDBasis`
  :param ndarray Y: The dependent variable array (the SAP flux)
  :param ndarray time: The independent variable array (the timestamps)
  :param ndarray errors: The standard errors on `Y`
  :param george.GP gp: The pre-initialized :py:mod:`george` gaussian process object
  :param list mask: The indices of points in `Y` to mask when computing the coefficients
  :param int niter: The number of iterations per principal component
  :param int nmasks: The number of masks to apply per iteration. The set of all these masks \
                     is our "validation" set
  
  :returns (masked_scatter, unmasked_scatter): A tuple containing the masked (validation) and unmasked (training) CDPP
  
  '''
  
  # Setup some variables
  masked_scatter = []
  mask_orig = np.array(mask, dtype = int)
  mask_new = np.array(mask_orig)
  
  # Mask transits and outliers
  M = Mask(mask_orig)
  mT = M(time)
  mY = M(Y)
  mX = M(X)
  mE = M(errors)
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
  
  '''
  EVEREST 0.1 BUG
  ---------------
  
  I think there's a bug in the for loop below. Above, we identify the indices ``inds``
  in the **masked** time array, but in the lines below we select those indices in the
  **unmasked** arrays. So what ends up happening is our validation sets aren't
  generally 13-cadence contiguous, especially if the light curve has many outliers.
  However, this doesn't make a terribly big difference. If anything, it
  leads to underfitting, since the CDPP in the validation set will be higher (since
  the set spans a larger time window), forcing the code to select fewer principal
  components. In the next version, the code below will be replaced with the following:
  
  for n in range(niter):
    masks = [np.arange(s, s + 13) for s in [np.random.choice(i) for i in inds if len(i)]]
    mask_new = np.concatenate(masks)
    C = PLDCoeffs(mX, mY, mT, mE, gp, mask_new)
    M = [PLDModel(C, mX[m]) for m in masks]
    masked_scatter.append(np.median([1.e6 * np.std((mY[m] - M[i] + med) / med) / np.sqrt(13) for i, m in enumerate(masks)]))  
  
  '''
  
  # Compute the precision several times and take the median
  for n in range(niter):
    
    # Get all our masks
    masks = [np.arange(s, s + 13) for s in [np.random.choice(i) for i in inds if len(i)]]

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
  the PLD coeffs. This performs iterative sigma-clipping.
  
  :param ndarray time: The independent variable array (the timestamps)
  :param ndarray flux: The dependent variable array (the SAP flux)
  :param ndarray fpix: The pixel fluxes, shape `(npts, npix)`
  :param ndarray ferr: The errors on the dependent variable array
  :param list mask: The indices of points in `Y` to mask when computing the coefficients
  :param float sigma: The outlier tolerance in standard deviations
  
  
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