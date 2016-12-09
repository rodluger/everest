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
from scipy.signal import savgol_filter, medfilt
import george
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import os, sys
import logging
log = logging.getLogger(__name__)

def GetChunk(time, breakpoints, b, mask = []):
  '''
  Returns the indices corresponding to a given light curve chunk.
  
  :param int b: The index of the chunk to return

  '''

  M = np.delete(np.arange(len(time)), mask, axis = 0)
  if b > 0:
    res = M[(M > breakpoints[b - 1]) & (M <= breakpoints[b])]
  else:
    res = M[M <= breakpoints[b]]
  return res

def fObj(res, **kwargs):
  '''
  SOM likelihood function. Nothing fancy.
  
  '''
  
  return -0.5 * np.nansum(res ** 2)

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

def GetLightCurves(campaign, module, model = 'nPLD', **kwargs):
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
  
  # Re-compute breakpoints
  segs = np.zeros_like(time)
  for i in range(1,len(time)):
    segs[i] = np.argmax(breakpoints > i - 1)
  segs[bad] = -1
  segs = np.delete(segs, bad)  
  breakpoints = np.where(np.diff(segs))[0]
  
  # Now delete the bad points from the arrays
  time_orig = np.array(time)
  time = np.delete(time, bad)
  fluxes = np.delete(np.array(fluxes), bad, axis = 1)
  
  # Add back the last breakpoint (signals end of light curve)
  breakpoints = np.append(breakpoints, [len(time) - 1])
  
  # -*- HACKS -*-
  if campaign == 0:
    breakpoints[0] = 501
  elif campaign == 1:
    breakpoints[0] = 1818
  
  return time_orig, time, breakpoints, fluxes

def SOM(time, fluxes, nflx, alpha_0 = 0.1, ki = 3, kj = 1, 
        niter = 30, periodic = False, **kwargs):
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
    maxf = np.nanmax(fluxes, axis = 0)
    minf = np.nanmin(fluxes, axis = 0)
    nflx = (fluxes - minf) / (maxf - minf)
  
  # Do `niter` iterations of the SOM
  for t in range(niter):
  
    # Loop over all light curves
    for n in range(N):
  
      # Find best matching pixel for this light curve
      bi, bj = np.unravel_index(np.nanargmax([fObj(nflx[n] - som[i][j]) 
                                           for i in range(ki) for j in range(kj)]), 
                                           (ki, kj))
    
      # Trick: flip light curve upside-down and see if we get a better match
      bi_, bj_ = np.unravel_index(np.nanargmax([fObj(1 - nflx[n] - som[i][j]) 
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

def GetRegressor(fluxes, nflx, som, **kwargs):
  '''
  
  '''
  
  # Find the best matching neuron for each light curve.
  # Store the indices in the `lcinds` map.
  ki, kj, sz = som.shape
  N = len(fluxes)
  lcinds = [[[] for j in range(kj)] for i in range(ki)]
  evratio = [[[] for j in range(kj)] for i in range(ki)]
  for n in range(N):
    bi, bj = np.unravel_index(np.argmin([np.nansum((nflx[n] - som[i][j]) ** 2) 
                                         for i in range(ki) for j in range(kj)]), 
                                         (ki, kj))
    lcinds[bi][bj].append(n)
  
  # We're going to use the SOM pixel with the largest number
  # of light curves as our CBV regressor
  lflat = [lcinds[i][j] for i in range(ki) for j in range(kj)]
  sflat = np.array([som[i][j] for i in range(ki) for j in range(kj)])
  pixel = np.nanargmax([len(n) for n in lflat])
  ij = np.unravel_index(pixel, (ki, kj))
    
  # Plot the results
  PlotSOM(nflx, som, lcinds, ij, **kwargs)
  
  return sflat[pixel]

def PlotSOM(nflx, som, lcinds, ij, figname = 'som.pdf', **kwargs):
  '''
  
  '''
  
  # Setup
  ki, kj, _ = som.shape
  N, sz = nflx.shape
  fig, ax = pl.subplots(1, figsize = (ki * 5, kj * 5))
  fig.subplots_adjust(left = 0.01, right = 0.99, bottom = 0.02, top = 0.98)
  
  # Plot the light curves, and SOMs / PCs in each pixel
  for i in range(ki):
    for j in range(kj):
      
      # The light curves
      for n in lcinds[i][j]:
        ax.plot(np.linspace(i + 0.5 - 0.4, i + 0.5 + 0.4, sz), 
                j + 0.5 - 0.4 + 0.8 * nflx[n], 
                'k-', alpha = min(1., 100. / N))
      
      # The Kohonen pixel
      ax.plot(np.linspace(i + 0.5 - 0.4, i + 0.5 + 0.4, sz), 
              j + 0.5 - 0.4 + 0.8 * som[i][j], 'r-', lw = 2)
      
      # The mean light curve, just for reference
      ax.plot(np.linspace(i + 0.5 - 0.4, i + 0.5 + 0.4, sz), 
              j + 0.5 - 0.4 + 0.8 * np.mean([nflx[n] for n in lcinds[i][j]], axis = 0), 'b-')
      
      # Indicate the number of light curves in this pixel
      ax.annotate('%d' % (len(lcinds[i][j])), 
                  xy = (i + 0.5, j + 0.95), 
                  xycoords = 'data', va = 'center', ha = 'center', 
                  zorder = 99, fontsize = 18, color = 'b' if ij == (i,j) else 'k')
      
  ax.get_xaxis().set_major_locator(MaxNLocator(integer=True))
  ax.get_yaxis().set_major_locator(MaxNLocator(integer=True))
  ax.grid(True, which = 'major', axis = 'both', linestyle = '-')
  ax.set_xlim(0, ki)
  ax.set_ylim(0, kj)
  ax.set_xticklabels([])
  ax.set_yticklabels([])
  fig.savefig(figname)
  pl.close()

def FitLightCurve(time, flux, breakpoints, X, mask = [], align = 'model'):
  '''
  
  '''
  
  # Loop over all the light curve segments
  model = [None for b in range(len(breakpoints))]
  weights = [None for b in range(len(breakpoints))]
  for b in range(len(breakpoints)):
    
    # Get the indices for this light curve segment
    inds = GetChunk(time, breakpoints, b, mask = mask)
        
    # Ordinary least squares
    mX = np.delete(X[b], mask, axis = 0)
    A = np.dot(mX.T, mX)
    B = np.dot(mX.T, flux[inds])
    weights[b] = np.linalg.solve(A, B)
    model[b] = np.dot(X[b], weights[b])
    
    # Vertical alignment
    if b == 0:
      model[b] -= np.nanmedian(model[b])
    else:
      if align == 'model':
        model[b] += (model[b - 1][-1] - model[b][0])
      elif align == 'flux':
        model[b] += (model[b - 1][-1] - model[b][0]) - (flux[inds[0] - 1] - flux[inds[0]])
      else:
        raise ValueError("Invalid value for `align`.")
  
  # Join model and normalize  
  model = np.concatenate(model)
  model -= np.nanmedian(model)
  
  fig, ax = pl.subplots(2, figsize = (12, 6), sharex = True)
  ax[0].plot(time, flux, 'k.', markersize = 3)
  ax[0].plot(time, model + np.nanmedian(flux), 'r-')
  ax[1].plot(time, flux - model, 'k.', markersize = 3)
  ax[0].set_ylabel('Original Flux', fontsize = 16)
  ax[1].set_ylabel('Corrected Flux', fontsize = 16)
  ax[1].set_xlabel('Time (BJD - 2454833)', fontsize = 16)
  ax[0].margins(0.01, None)
  
  pl.show()

def Compute(time, fluxes, breakpoints, smooth_in = True, smooth_out = True, 
            nreg = 2, path = '', **kwargs):
  '''
  
  '''
  
  X = [None for b in range(len(breakpoints))]
  for b in range(len(breakpoints)):
    
    # Get the indices for this light curve segment
    inds = GetChunk(time, breakpoints, b)
    t = time[inds]
    f = np.array(fluxes[:, inds])
    
    # Smooth the fluxes to compute the SOM?
    if smooth_in:
      for i in range(len(f)):
        f[i] = savgol_filter(f[i], 99, 2)
    
    # The design matrix for the current chunk
    X[b] = np.ones_like(t).reshape(-1, 1)
    
    # Run `nreg` iterations
    for n in range(nreg):
      
      log.info('Computing: segment %d/%d, iteration %d/%d...' %
              (b + 1, len(breakpoints), n + 1, nreg))
      
      # Normalize the fluxes
      maxf = np.nanmax(f, axis = 1).reshape(-1, 1)
      minf = np.nanmin(f, axis = 1).reshape(-1, 1)
      nflx = (f - minf) / (maxf - minf)

      # Compute the SOM
      som = SOM(t, f, nflx, **kwargs)
    
      # Now get the regressor
      cbv = GetRegressor(f, nflx, som, figname = 
                         os.path.join(path, 'Seg%02d_Iter%02d.pdf' % 
                         (b + 1, n + 1)), **kwargs)
      if smooth_out:
        cbv = savgol_filter(cbv, 99, 2)
      
      # Append to our design matrix
      X[b] = np.hstack([X[b], cbv.reshape(-1, 1)])
      
      # Recursion. We de-trend each of the light curves with
      # our design matrix so far, then repeat.
      if n < nreg - 1:
        A = np.dot(X[b].T, X[b])
        for i in range(len(f)):
          f[i] = fluxes[i, inds] - np.dot(X[b], np.linalg.solve(A, np.dot(X[b].T, fluxes[i, inds])))

  return X

def Run(campaign, module, model = 'nPLD', clobber = False, quiet = False, **kwargs):
  '''
  
  '''
  
  # Output path
  path = os.path.join(EVEREST_DAT, 'k2', 'cbv', 'c%02d' % campaign, '%2d' % module, model)
  if not os.path.exists(path):
    os.makedirs(path)
  
  # Set up the log file
  root = logging.getLogger()
  root.handlers = []
  root.setLevel(logging.DEBUG)
  sh = logging.StreamHandler(sys.stdout)
  if quiet:
    sh.setLevel(logging.CRITICAL)
  else:
    sh.setLevel(logging.DEBUG)
  sh_formatter = logging.Formatter("[%(funcName)s()]: %(message)s")
  sh.setFormatter(sh_formatter)
  root.addHandler(sh)
  
  # Get the design matrix
  xfile = os.path.join(path, 'X.npz')
  if clobber or not os.path.exists(xfile):
    
    # Get the light curves
    log.info('Obtaining light curves...')
    lcfile = os.path.join(path, 'lcs.npz')
    if clobber or not os.path.exists(lcfile):
      time_orig, time, breakpoints, fluxes = GetLightCurves(campaign, module, model = model, **kwargs)
      np.savez(lcfile, time_orig = time_orig, time = time, breakpoints = breakpoints, fluxes = fluxes)
    else:
      lcs = np.load(lcfile)
      time_orig = lcs['time_orig']
      time = lcs['time']
      breakpoints = lcs['breakpoints']
      fluxes = lcs['fluxes']
    
    # Get the design matrix  
    X = Compute(time, fluxes, breakpoints, path = path, **kwargs)
    np.savez(xfile, X = X)
  
  else:
      
    X = np.load(xfile)['X']
  
  
  # DEBUG
  lcfile = os.path.join(path, 'lcs.npz')
  lcs = np.load(lcfile)
  time = lcs['time']
  breakpoints = lcs['breakpoints']
  fluxes = lcs['fluxes']
  for f in fluxes:
    FitLightCurve(time, f, breakpoints, X)
  