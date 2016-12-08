#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`cbv.py`
----------------

**EXPERIMENTAL**

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from ...config import EVEREST_DAT
from .aux import GetK2Campaign
import numpy as np
try:
  from wpca import WPCA
except:
  WPCA = None
from scipy.signal import savgol_filter, medfilt
import george
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import os
import logging
log = logging.getLogger(__name__)

def GetChunk(time, breakpoints, b):
  '''
  Returns the indices corresponding to a given light curve chunk.
  
  :param int b: The index of the chunk to return

  '''

  M = np.arange(len(time))
  if b > 0:
    res = M[(M > breakpoints[b - 1]) & (M <= breakpoints[b])]
  else:
    res = M[M <= breakpoints[b]]
  return res

def fObj(res, **kwargs):
  '''
  SOM likelihood function. Nothing fancy.
  
  '''
  
  return -0.5 * np.sum(res ** 2)

def Channels(module):
  '''
  Returns the channels contained in the given module.
  
  '''
  
  nums = {2:1, 3:5, 4:9, 6:13, 7:17, 8:21, 9:25, 
          10:29, 11:33, 12:37, 13:41, 14:45, 15:49, 
          16:53, 17:57, 18:61, 19:65, 20:69, 22:73, 
          23:77, 24:81}
  
  if module in nums:
    return [nums[module], nums[module] + 1, 
            nums[module] + 2, nums[module] + 3]
  else:
    return None

def GetLightCurves(campaign, model = 'nPLD', module = range(99), **kwargs):
  '''
  
  '''
  
  # Get the channel numbers
  channels = []
  for m in np.atleast_1d(module):
    c = Channels(m)
    if c is not None:
      channels.extend(c)
  
  # Get the EPIC numbers
  all = GetK2Campaign(campaign)
  stars = np.array([s[0] for s in all if s[2] in channels and 
          os.path.exists(
          os.path.join(EVEREST_DAT, 'k2', 'c%02d' % int(campaign),
          ('%09d' % s[0])[:4] + '00000', 
          ('%09d' % s[0])[4:], model + '.npz'))], dtype = int)
  N = len(stars)
  assert N > 0, "No light curves found."
  
  # Loop over all stars and store the fluxes in a list
  fluxes = []
  masks = []
  for n in range(N):
    
    print("Processing light curve %d/%d..." % (n + 1, N))
    
    # De-trended light curve file name
    nf = os.path.join(EVEREST_DAT, 'k2', 'c%02d' % int(campaign),
                   ('%09d' % stars[n])[:4] + '00000', 
                   ('%09d' % stars[n])[4:], model + '.npz')
    
    # Get the data
    data = np.load(nf)
    t = data['time']
    if n == 0:
      time = t
      breakpoints = data['breakpoints']
    
    # Get de-trended light curve
    y = data['fraw'] - data['model'] 
    
    # Interpolate over outliers 
    m = np.array(list(set(np.concatenate([data['outmask'], data['badmask'], 
                                          data['nanmask'], data['transitmask']]))), 
                                          dtype = int)
    y = np.interp(t, np.delete(t, m), np.delete(y, m))
    
    # Throw out light curves with more than 20% variation 
    # above/below the median, as they could mess up PCA
    my = np.nanmedian(y)
    if len(np.where((y < 0.8 * my) | (y > 1.2 * my))[0]):
      continue
    
    # Append fluxes to our running list
    fluxes.append(y)
    
    # Append the outlier mask to our running list
    mask = np.zeros_like(y)
    mask[m] = 1
    masks.append(mask)
  
  # Remove data points that are missing in all light curves
  masks = np.array(masks)
  bad = np.where(np.sum(masks, axis = 0) == masks.shape[0])[0]
  
  # Re-compute breakpoints.
  segs = np.zeros_like(time)
  for i in range(1,len(time)):
    segs[i] = np.argmax(breakpoints > i - 1)
  segs[bad] = -1
  segs = np.delete(segs, bad)  
  breakpoints = np.where(np.diff(segs))[0]
  
  # Now delete the bad points from the arrays
  time = np.delete(time, bad)
  fluxes = np.delete(np.array(fluxes), bad, axis = 1)
  
  # Add back the last breakpoint (signals end of light curve)
  breakpoints = np.append(breakpoints, [len(time) - 1])
  
  # -*- HACK -*-
  # We need to shift the breakpoint a little for C00 so that it 
  # lines up with the data gap. It's a long story, but this
  # hack fixes an issue with discontinuities.
  if campaign == 0:
    breakpoints[0] += 2
  
  return time, breakpoints, fluxes

def SOM(time, fluxes, nflx, alpha_0 = 0.1, ki = 3, kj = 1, 
        niter = 100, periodic = False, **kwargs):
  '''
  
  ''' 
   
  # Initialize the SOM
  N = len(fluxes)
  t = 0
  sig_0 = max(ki, kj)
  alpha = alpha_0
  sig = sig_0
  som = np.random.random((ki, kj, len(time))) 
    
  # Normalize the fluxes
  if nflx is None:
    maxf = np.max(fluxes, axis = 0)
    minf = np.min(fluxes, axis = 0)
    nflx = (fluxes - minf) / (maxf - minf)
  
  # Do `niter` iterations of the SOM
  for t in range(niter):
  
    print("Iteration %d/%d..." % (t + 1, niter))
  
    # Loop over all light curves
    for n in range(N):
  
      # Find best matching pixel for this light curve
      bi, bj = np.unravel_index(np.argmax([fObj(nflx[n] - som[i][j]) 
                                           for i in range(ki) for j in range(kj)]), 
                                           (ki, kj))
    
      # Trick: flip light curve upside-down and see if we get a better match
      bi_, bj_ = np.unravel_index(np.argmax([fObj(1 - nflx[n] - som[i][j]) 
                                             for i in range(ki) for j in range(kj)]), 
                                             (ki, kj))
      if fObj(1 - nflx[n] - som[bi_][bj_]) > fObj(nflx[n] - som[bi][bj]):
        nflx[n] = 1 - nflx[n]
        bi, bj = bi_, bj_
    
      # Update the Kohonen layer
      for i in range(ki):
        for j in range(kj):
        
          # Compute distance
          if periodic:
            dx2 = min((bi - i) ** 2, (ki - (bi - i)) ** 2)
            dy2 = min((bj - j) ** 2, (kj - (bj - j)) ** 2)
          else:
            dx2 = (bi - i) ** 2
            dy2 = (bj - j) ** 2
          d2 = dx2 + dy2
        
          # Update the pixel
          som[i][j] += alpha * np.exp(-d2 / (2 * sig ** 2)) * (nflx[n] - som[i][j])
  
    # Update the params
    sig = sig_0 * np.exp(-t / niter * np.log(max(ki, kj)))
    alpha = alpha_0 * (1 - t / niter)
  
  return som

def GetCBVs(fluxes, nflx, som, nreg = 3, pca = False, plot = False, **kwargs):
  '''
  
  '''
  
  # Find the best matching neuron for each light curve.
  # Store the indices in the `lcinds` map.
  ki, kj, sz = som.shape
  N = len(fluxes)
  lcinds = [[[] for j in range(kj)] for i in range(ki)]
  evratio = [[[] for j in range(kj)] for i in range(ki)]
  for n in range(N):
    bi, bj = np.unravel_index(np.argmin([np.sum((nflx[n] - som[i][j]) ** 2) 
                                         for i in range(ki) for j in range(kj)]), 
                                         (ki, kj))
    lcinds[bi][bj].append(n)
  
  # Get our CBVs
  if pca:
  
    # Get the top `nreg` principal components in each pixel
    pc = [[[] for j in range(kj)] for i in range(ki)]
    for i in range(ki):
      for j in range(kj):
        n = lcinds[i][j]
        pca = WPCA(n_components = nreg)
        pc[i][j] = pca.fit_transform((fluxes[n] / np.nanmedian(fluxes[n], 
                                      axis = 1).reshape(-1, 1)).T, 
                                      weights = np.array([np.nanmedian(np.sqrt(f)) * 
                                      np.ones_like(f) for f in fluxes[n]]).T).T
        evratio[i][j] = pca.explained_variance_ratio_
  
    # We're going to guess that the SOM pixel corresponding
    # to the instrumental component is the one with the largest
    # value of `n * evratio`, where `n` is the number of light curves
    # in the pixel and `evratio` is the explained variance ratio
    # of the first principal component. This isn't rigorous, but
    # seems to work in practice. We then create a simple design
    # matrix out of the SOM array in that pixel for fitting.
    ij = np.argmax([len(lcinds[i][j]) * evratio[i][j][0] 
                    for i in range(ki) for j in range(kj)])
    i, j = np.unravel_index(ij, (ki, kj))
    X = np.hstack([pc[i][j].T, np.ones(sz).reshape(-1, 1)])
    pixels = [(i,j),]
    
  else:
  
    # We're going to use the SOM arrays themselves as CBVs, ordered by
    # the number of light curves in each
    pc = None
    lflat = [lcinds[i][j] for i in range(ki) for j in range(kj)]
    sflat = np.array([som[i][j] for i in range(ki) for j in range(kj)])
    pixels = np.argsort([-len(n) for n in lflat])[:nreg]
    X = np.hstack([sflat[pixels].T, np.ones((sz, 1))])
    pixels = [np.unravel_index(p, (ki, kj)) for p in pixels]
    
  if plot:
    PlotSOM(nflx, som, lcinds, pc, pixels, **kwargs)
  
  return X

def PlotSOM(nflx, som, lcinds, pc, pixels, **kwargs):
  '''
  
  '''
  
  # Setup
  ki, kj, _ = som.shape
  N, sz = nflx.shape
  fig, ax = pl.subplots(1, figsize = (ki * 5, kj * 5))
  
  # Plot the light curves, and SOMs / PCs in each pixel
  for i in range(ki):
    for j in range(kj):
      
      # The light curves
      for n in lcinds[i][j]:
        ax.plot(np.linspace(i + 0.5 - 0.4, i + 0.5 + 0.4, sz), 
                j + 0.5 - 0.4 + 0.8 * nflx[n], 
                'k-', alpha = min(1., 100. / N))
      
      if pc is not None:
      
        # The principal components
        for p, c in zip(pc[i][j], ['r', 'b', 'g', 'y', 'purple']):
          ax.plot(np.linspace(i + 0.5 - 0.4, i + 0.5 + 0.4, sz), 
                  j + 0.5 - 0.4 + 0.8 * (p - np.min(p)) / (np.max(p) - np.min(p)), 
                  ls = '-', color = c)
                  
      else:
      
        # The Kohonen pixel
        ax.plot(np.linspace(i + 0.5 - 0.4, i + 0.5 + 0.4, sz), 
                j + 0.5 - 0.4 + 0.8 * som[i][j], 'r-')
      
      # Indicate the number of light curves in this pixel
      ax.annotate('%d' % (len(lcinds[i][j])), 
                  xy = (i + 0.5, j + 0.95), 
                  xycoords = 'data', va = 'center', ha = 'center', 
                  zorder = 99, fontsize = 18, color = 'b' if (i,j) in pixels else 'k')
      
      
  ax.get_xaxis().set_major_locator(MaxNLocator(integer=True))
  ax.get_yaxis().set_major_locator(MaxNLocator(integer=True))
  ax.grid(True, which = 'major', axis = 'both', linestyle = '-')
  ax.set_xlim(0, ki)
  ax.set_ylim(0, kj)
  pl.show()

def Run(campaign, clobber = False, smooth_in = True, smooth_out = True, **kwargs):
  '''
  
  '''
  
  # Get the light curves
  if clobber or not os.path.exists('lcs.npz'):
    time, breakpoints, fluxes = GetLightCurves(campaign, **kwargs)
    np.savez('lcs.npz', time = time, breakpoints = breakpoints, fluxes = fluxes)
  else:
    lcs = np.load('lcs.npz')
    time = lcs['time']
    breakpoints = lcs['breakpoints']
    fluxes = lcs['fluxes']
  
  # DEBUG
  quit()
    
  X = [None for b in range(len(breakpoints))]
  for b in range(len(breakpoints)):
    
    # Get the indices for this light curve segment
    inds = GetChunk(time, breakpoints, b)
    t = time[inds]
    f = fluxes[:, inds]
    
    # Smooth the fluxes to compute the SOM?
    if smooth_in:
      for n in range(len(f)):
        f[n] = savgol_filter(f[n], 99, 2)
    
    # Normalize the fluxes
    maxf = np.max(f, axis = 1).reshape(-1, 1)
    minf = np.min(f, axis = 1).reshape(-1, 1)
    nflx = (f - minf) / (maxf - minf)

    # Compute the SOM
    som = SOM(t, f, nflx, **kwargs)
    
    # Now get the CBVs
    X[b] = GetCBVs(f, nflx, som, **kwargs)
    
    # Smooth the regressors?
    if smooth_out:
      for n in range(X[b].shape[1] - 1):
        X[b][:,n] = savgol_filter(X[b][:,n], 99, 2)
    
  for n in range(len(fluxes)):
    
    model = [None for b in range(len(breakpoints))]
    for b in range(len(breakpoints)):
      
      # Get the indices for this light curve segment
      inds = GetChunk(time, breakpoints, b)
      
      # GP regression
      gp = george.GP((0.005 * np.nanmedian(fluxes[n,inds])) ** 2 * george.kernels.Matern32Kernel(1000. ** 2))
      gp.compute(time[inds], np.ones_like(time[inds]) * np.nanmedian(np.sqrt(fluxes[n,inds])))
      A = np.dot(X[b].T, gp.solver.apply_inverse(X[b]))
      B = np.dot(X[b].T, gp.solver.apply_inverse(fluxes[n,inds]))
      C = np.linalg.solve(A, B)
      model[b] = np.dot(X[b], C)
    model = np.concatenate(model)
    
    fig, ax = pl.subplots(2, figsize = (12, 6))
    ax[0].plot(time, fluxes[n], 'k.', markersize = 3)
    ax[0].plot(time, model, 'r-')
    ax[1].plot(time, fluxes[n] - model + np.nanmedian(fluxes[n]), 'k.', markersize = 3)
    pl.show()