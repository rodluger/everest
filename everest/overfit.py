#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`overfit.py` - Overfitting metrics
------------------------------------------

EXPERIMENTAL!

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .transit import TransitShape
from .utils import prange
import numpy as np
import matplotlib.pyplot as pl
try:
  from choldate import cholupdate, choldowndate
except:
  cholupdate = None
  choldowndate = None
from scipy.linalg import cholesky, cho_solve
import logging
log = logging.getLogger(__name__)

class Overfit(object):
  '''
  Generic overfitting metric class.
  
  '''
  
  def __init__(self, t, y, X, XL, XLX, C, Kinv, mask = [], tau = None, w = 9, rank2 = True, **kwargs):
    '''
    
    '''
    
    if cholupdate is None or choldowndate is None:
      raise Exception("Please install the `choldate` package.")
    self._unmasked(t, y, C, XLX, Kinv, mask = mask, tau = tau, **kwargs)
    self._masked(t, y, X, XL, C, Kinv, mask = mask, w = w, tau = tau, rank2 = rank2, **kwargs)
    self.time = t
    self.mask = mask
    
  def plot_unmasked(self, depth = None):
    '''
  
    '''
    
    if depth:
      metric = 1 - (self._O2 + self._O3 / depth) / self._O1
    else:
      metric = 1 - self._O2 / self._O1
    
    # Set up
    fig = pl.figure(figsize = (10, 4))
    fig.subplots_adjust(bottom = 0.2)
    ax = pl.subplot2grid((1, 5), (0, 0), colspan = 4, rowspan = 1)
    axh = pl.subplot2grid((1, 5), (0, 4), colspan = 1, rowspan = 1)
    ax.set_xlabel('Time (days)', fontsize = 14)
    ax.set_ylabel('Unmasked overfitting', fontsize = 14)
    
    # Plot the metric as a function of time
    ax.plot(np.delete(self.time, self.mask), np.delete(metric, self.mask), 'k.', alpha = 0.5, ms = 2)
    ax.plot(np.delete(self.time, self.mask), np.delete(metric, self.mask), 'k-', alpha = 0.1, lw = 0.5)
    ax.plot(self.time[self.mask], metric[self.mask], 'r.', alpha = 0.25, ms = 2)
    
    # Bound 99% of data
    y = np.delete(metric, self.mask)
    N = int(0.99 * len(y))
    hi, lo = y[np.argsort(y)][[N,-N]]
    fsort = y[np.argsort(y)]
    pad = (hi - lo) * 0.1
    ylim = (lo - pad, hi + pad)
    ax.set_ylim(ylim)
    
    # Plot the histogram
    m1 = np.delete(metric, self.mask)
    m2 = metric
    axh.hist(m1, bins = 30, range = ylim, orientation = "horizontal", histtype = "step", fill = False, color = 'k')
    axh.hist(m2, bins = 30, range = ylim, orientation = "horizontal", histtype = "step", fill = False, color = 'r', alpha = 0.3)
    axh.yaxis.tick_right()
    axh.set_ylim(*ax.get_ylim())
    axh.set_xticklabels([])
    
    # Return
    return fig, ax, axh
      
  def plot_masked(self, depth = 0.001):
    '''
    
    '''
    
    # Overfitting metric
    metric = self._O5 / depth
    
    # Set up
    fig = pl.figure(figsize = (10, 4))
    fig.subplots_adjust(bottom = 0.2)
    ax = pl.subplot2grid((1, 5), (0, 0), colspan = 4, rowspan = 1)
    axh = pl.subplot2grid((1, 5), (0, 4), colspan = 1, rowspan = 1)
    ax.set_xlabel('Time (days)', fontsize = 14)
    ax.set_ylabel('Masked overfitting', fontsize = 14)
    
    # Plot the metric as a function of time
    ax.plot(np.delete(self.time, self.mask), np.delete(metric, self.mask), 'k.', alpha = 0.5, ms = 2)
    ax.plot(np.delete(self.time, self.mask), np.delete(metric, self.mask), 'k-', alpha = 0.1, lw = 0.5)
    ax.plot(self.time[self.mask], metric[self.mask], 'r.', alpha = 0.25, ms = 2)
    
    # Bound 99% of data
    y = np.delete(metric, self.mask)
    N = int(0.99 * len(y))
    hi, lo = y[np.argsort(y)][[N,-N]]
    fsort = y[np.argsort(y)]
    pad = (hi - lo) * 0.1
    ylim = (lo - pad, hi + pad)
    ax.set_ylim(ylim)
    
    # Plot the histogram
    m1 = np.delete(metric, self.mask)
    m2 = metric
    axh.hist(m1, bins = 30, range = ylim, orientation = "horizontal", histtype = "step", fill = False, color = 'k')
    axh.hist(m2, bins = 30, range = ylim, orientation = "horizontal", histtype = "step", fill = False, color = 'r', alpha = 0.3)
    axh.yaxis.tick_right()
    axh.set_ylim(*ax.get_ylim())
    axh.set_xticklabels([])
    
    # Return
    return fig, ax, axh
  
  def plot_joint(self, depth = 0.001):
    '''
  
    '''
  
    # TODO
    pass
  
  def plot_marginal(self, depth = 0.001, std = 0.001, mu = 0.001):
    '''
    
    '''
    
    # TODO
    pass

  def _unmasked(self, t, y, C, XLX, Kinv, mask = [], tau = None, **kwargs):
    '''
    Computes the unmasked overfitting metric.
    
    '''
  
    # Some masked arrays/matrices
    mt = np.delete(t, mask)
    my = np.delete(y, mask)
    MC = np.delete(np.delete(C, mask, axis = 0), mask, axis = 1)
  
    # Compute the full `model` and the de-trended flux `f`
    model = np.dot(XLX, np.linalg.solve(MC, my))
    model -= np.nanmedian(model)
    f = y - model
  
    # The transit model
    if tau is None:
      tau = TransitShape(dur = 0.1)
  
    # The regression matrix
    log.info("Unmasked light curve: pre-computing some stuff...")
    R = np.linalg.solve(MC, XLX.T).T
  
    # Arrays we'll use for the overfitting calculations
    self._O1 = np.zeros_like(t)
    self._O2 = np.zeros_like(t)
    self._O3 = np.zeros_like(t)
    self._O4 = np.zeros_like(t)
  
    # Loop over all cadences
    log.info("Unmasked light curve: looping through the dataset...")
    for k in prange(len(t)):
    
      # Evaluate the transit model...
      TAU = tau(t, t0 = t[k])
      i = np.where(TAU < 0)[0]
      TAU = TAU.reshape(-1,1)
    
      # ...and the masked transit model
      MTAU = tau(mt, t0 = t[k])
      mi = np.where(MTAU < 0)[0]
      MTAU = MTAU.reshape(-1,1)
    
      # Sparse algebra
      A = np.dot(np.dot(TAU[i].T, Kinv[i,:][:,i]), TAU[i])
      B = np.dot(TAU[i].T, Kinv[i,:])
      C = TAU - np.dot(R[:,mi], MTAU[mi])
    
      # Compute the overfitting metric
      self._O1[k] = A
      self._O2[k] = np.dot(B, C)
      self._O3[k] = np.dot(B, f)
      self._O4[k] = np.dot(B, y)

    return
  
  def _default_mask(self, w, N, mask):
    '''
    Returns the mask for the arrays and matrices 
    in the zeroth iteration of the Cholesky update process.
    Also returns `r`, the number of outliers in the
    first `w` cadences.
  
    '''
  
    # The mask applied to the default timeseries
    mask0 = list(np.arange(w)) + list(mask)
  
    # For each outlier in the first `w` cadences, add
    # the last unmasked cadence to the default mask
    r = 0
    j = N - 1
    for i in mask:
      if i in range(w):
        while j in mask0: 
          j -= 1
        mask0.append(j)
        r += 1
  
    # Sort the mask
    mask0 = np.array(sorted(list(set(mask0))), dtype = int)
  
    return mask0, r
  
  def _masked(self, t, y, X, XL, C, Kinv, mask = [], w = 9, tau = None, rank2 = True, **kwargs):
    '''
    Computes the masked overfitting metric.
  
    '''
  
    # Get the data for an Everest target
    # `t` is the time array
    # `y` is the normalized raw flux
    # `X` is the PLD design matrix
    # `XL` is the regularized, unmasked design matrix
    # `C` is the normalized full covariance
    # `Kinv` is the inverse astrophysical covariance
    N = len(y)

    # Get the default mask and apply it
    log.info("Masked light curve: pre-computing some stuff...")
    mask0, r = self._default_mask(w, N, mask)
    # Construct a version of the data where
    # the first `w` indices are masked. This is the 
    # "0" dataset. `C0` is the full covariance 
    # matrix for this dataset.
    C0 = np.delete(np.delete(C, mask0, axis = 0), mask0, axis = 1)
    # Compute the upper Cholesky factorization for the "0"
    # dataset. We will update this iteratively.
    U_n = cholesky(C0)
    # A term in the model equation and the data vector, which
    # we will also update iteratively.
    P_n = np.dot(XL, np.delete(X, mask0, axis = 0).T)
    y_n = np.delete(y, mask0)
    # Indices (for bookkeeping)
    inds = np.arange(N)
    inds_n = np.delete(inds, mask0)

    # Define our replacement rules `iold` and `inew`.
    rep = np.delete(inds, mask)[::-1]
    rev = False
    iold = np.array(inds)[w:]
    inew = np.array(inds)[:N - w]

    # First, flag the indices where we will
    # need to make replacements due to outliers.
    fold = np.zeros_like(iold, dtype = bool)
    fnew = np.zeros_like(inew, dtype = bool)
    for m in mask:
      fold[np.where(iold == m)] = True
      fnew[np.where(inew == m)] = True

    # Now, go through the timeseries and amend
    # the replacement rules when there are outliers.
    for n in range(N - w):

      # Old index
      if fold[n]:
        iold[n] = rep[r]
        r += 1
    
      # New index
      if fnew[n]:
        r -= 1
        inew[n] = rep[r]
  
      # If we're past the halfway point, we will
      # borrow data from the beginning rather than
      # from the end.
      if (not rev) and (n > (N // 2)) and (r == 0):
        rep = rep[::-1]
        rev = True

    # The indices in `data_n` at which we are storing each
    # element of the original timeseries `data` make up the
    # index map `imap`
    imap = np.zeros(N, dtype = int) - 1
    for n, i in enumerate(inds_n):
      imap[i] = n
    
    # The overfitting metric(s)
    self._O5 = np.zeros(N) * np.nan
  
    # The transit model
    if tau is None:
      tau = TransitShape(dur = 0.1)
  
    # The base case. Uncommented, but it's what you'd
    # get if n = -1 in the loop below, so check out
    # the comments down there.
    j = (w + 1) // 2 - 1
    M = cho_solve((U_n, False), y_n)  
    m = np.dot(P_n, M)
    TAU = tau(t, t0 = t[j])
    i = np.where(TAU < 0)[0]
    TAU = TAU[i].reshape(-1,1)
    den = np.dot(np.dot(TAU.T, Kinv[i,:][:,i]), TAU)
    num = np.dot(TAU.T, Kinv[i,:])
    self._O5[:j+1] = -np.dot(num, y - m) / den
  
    # Iterate!  
    log.info("Masked light curve: looping through the dataset...")
    for n in prange(N - 2 * w):
            
      # The mask for this iteration goes from `n + 1` to `n + w`
      # (endpoints inclusive).
      # The center point of the mask is therefore
      j = n + (w + 1) // 2
        
      if iold[n] == inew[n]:
      
        # No change relative to the last iteration!
        self._O5[j] = self._O5[j - 1]
    
      else:
        
        # Update one element in our data and our index list
        inds_n[imap[iold[n]]] = inds[inew[n]]
        # Note that
        # y_n = y[inds_n] and
        # P_n = np.dot(XL, X[inds_n,:].T)
        y_n[imap[iold[n]]] = y[inew[n]]
        P_n[:,imap[iold[n]]] = np.dot(X[inew[n],:], XL.T)  

        # Now we update the covariance.
        # If no outliers are masked, this is equivalent to
        # R = np.delete(C[inew[n]] - C[iold[n]], range(n, n + w))
        # R[n] = C[inew[n], inew[n]] - C[iold[n], iold[n]]
        R = C[inew[n], inds_n] - C[iold[n], inds_n]
        R[imap[iold[n]]] = C[inew[n], inew[n]] - C[iold[n], iold[n]]
                
        # Rank 2 Cholesky update, R[n] positive
        if (rank2) and (R[imap[iold[n]]] > 0):
      
          # 1. Update
          v = np.array(R) / np.sqrt(R[imap[iold[n]]])
          cholupdate(U_n, v)
      
          # 2. Downdate
          v = np.array(R) / np.sqrt(R[imap[iold[n]]])
          v[imap[iold[n]]] = 0
          choldowndate(U_n, v)
    
        # Rank 2 Cholesky update, R[n] negative
        elif (rank2) and (R[imap[iold[n]]] < 0):
      
          # 1. Update
          v = np.array(R) / np.sqrt(-R[imap[iold[n]]])
          v[imap[iold[n]]] = 0
          cholupdate(U_n, v)
      
          # 2. Downdate
          v = np.array(R) / np.sqrt(-R[imap[iold[n]]])
          choldowndate(U_n, v)
    
        # Rank 3 Cholesky update. Slower but may be more stable
        else:
      
          # 1. Update
          v = np.array(R)
          v[imap[iold[n]]] = 1
          cholupdate(U_n, v)

          # 2. Downdate
          v = np.array(R)
          v[imap[iold[n]]] = 0
          choldowndate(U_n, v)

          # 3. Final up/downdate
          v = np.zeros_like(R)
          if R[imap[iold[n]]] >= 1:
            v[imap[iold[n]]] = np.sqrt(R[imap[iold[n]]] - 1)
            cholupdate(U_n, v)
          else:
            v[imap[iold[n]]] = np.sqrt(1 - R[imap[iold[n]]])
            choldowndate(U_n, v)
          
        # Now we can solve our linear equation
        try:
          M = cho_solve((U_n, False), y_n)
        except:
          import pdb; pdb.set_trace()
        m = np.dot(P_n, M)
        
        # Evaluate the sparse transit model
        TAU = tau(t, t0 = t[j])
        i = np.where(TAU < 0)[0]
        TAU = TAU[i].reshape(-1,1)
  
        # Dot the transit model in
        den = np.dot(np.dot(TAU.T, Kinv[i,:][:,i]), TAU)
        num = np.dot(TAU.T, Kinv[i,:])

        # Compute the overfitting metric
        # Divide this number by a depth
        # to get the overfitting for that
        # particular depth.
        self._O5[j] = -np.dot(num, y - m) / den
    
        # Finally, update the index map
        imap[inew[n]] = imap[iold[n]]
        imap[iold[n]] = -1
  
    # Fill in the endpoints
    self._O5[j:] = self._O5[j]
  
class SavedOverfitObject(Overfit):
  '''
  
  '''
  
  def __init__(self, time, mask, O1, O2, O3, O4, O5):
    '''
    
    '''
    
    self._O1 = O1
    self._O2 = O2
    self._O3 = O3
    self._O4 = O4
    self._O5 = O5
    self.time = time
    self.mask = mask