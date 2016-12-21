#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
sysrem.py
---------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from ...config import EVEREST_DAT
from .aux import GetK2Campaign, Campaign, Channels
import os
import numpy as np
import matplotlib.pyplot as pl
import george
from george.kernels import Matern32Kernel, WhiteKernel
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

def GetStars(campaign, module, model = 'nPLD', **kwargs):
  '''
  
  '''
  
  # Get the channel numbers
  if module == 'all':
    channels = range(99)
  else:
    channels = Channels(module)
  
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
  errors = []
  kpars = []
  
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
    err = data['fraw_err']
    
    # De-weight outliers and bad timestamps 
    m = np.array(list(set(np.concatenate([data['outmask'], data['badmask'], 
                                          data['nanmask'], data['transitmask']]))), 
                                          dtype = int)
    
    # Interpolate over the outliers
    y = np.interp(t, np.delete(t, m), np.delete(y, m))
    err = np.interp(t, np.delete(t, m), np.delete(err, m))
       
    # Append to our running lists
    fluxes.append(y)
    errors.append(err)
    kpars.append(data['kernel_params'])
    
  return time, breakpoints, np.array(fluxes), np.array(errors), np.array(kpars)

def SysRem(time, flux, err, nrec = 5, niter = 50, **kwargs):
  '''
  
  '''
  
  nflx, tlen = flux.shape
  
  # Get normalized fluxes
  med = np.nanmedian(flux, axis = 1).reshape(-1, 1)
  y = flux - med
  
  # Compute the inverse of the variances
  invvar = 1. / err ** 2
  
  # The CBVs for this set of fluxes
  cbvs = np.zeros((nrec, tlen))

  # Recover `nrec` components
  for n in range(nrec):
    
    # Initialize the weights and regressors
    c = np.zeros(nflx)
    a = np.ones(tlen)
    f = y * invvar
    
    # Perform `niter` iterations
    for i in range(niter):
    
      # Compute the `c` vector (the weights)
      c = np.dot(f, a) / np.dot(invvar, a ** 2)

      # Compute the `a` vector (the regressors)
      a = np.dot(c, f) / np.dot(c ** 2, invvar)
  
    # Remove this component from all light curves
    y -= np.outer(c, a)
    
    # Save this regressor
    cbvs[n] = a - np.nanmedian(a)
    
  return cbvs

def GetCBVs(campaign, module, model = 'nPLD', clobber = False, **kwargs):
  '''
  
  '''
  
  # Output path
  path = os.path.join(EVEREST_DAT, 'k2', 'cbv', 'c%02d' % campaign, str(module), model)
  if not os.path.exists(path):
    os.makedirs(path)
  
  # Get the design matrix
  xfile = os.path.join(path, 'X.npz')
  if clobber or not os.path.exists(xfile):
    
    # Get the light curves
    log.info('Obtaining light curves...')
    lcfile = os.path.join(path, 'lcs.npz')
    if clobber or not os.path.exists(lcfile):
      time, breakpoints, fluxes, errors, kpars = GetStars(campaign, module, model = model, **kwargs)
      np.savez(lcfile, time = time, breakpoints = breakpoints, fluxes = fluxes, errors = errors, kpars = kpars)
    else:
      lcs = np.load(lcfile)
      time = lcs['time']
      breakpoints = lcs['breakpoints']
      fluxes = lcs['fluxes']
      errors = lcs['errors']
      kpars = lcs['kpars']
    
    # Compute the design matrix  
    X = np.ones((len(time), 1 + kwargs.get('nrec', 5)))
    
    # Loop over the segments
    new_fluxes = np.zeros_like(fluxes)
    for i, b in enumerate(range(len(breakpoints))):
      
      # Get the current segment's indices
      inds = GetChunk(time, breakpoints, b)
    
      # Update the error arrays with the white GP component
      for j in range(len(errors)):
        errors[j] = np.sqrt(errors[j] ** 2 + kpars[j][0] ** 2)
    
      # Get de-trended fluxes
      X[inds,1:] = SysRem(time[inds], fluxes[:,inds], errors[:,inds])
    
    # Save
    np.savez(xfile, X = X)
  
  else:
    
    # Load from disk
    X = np.load(xfile)['X']
  
  return X
 
def GetLightCurve(EPIC, campaign, model = 'nPLD', **kwargs):
  '''
  
  '''
  
  # De-trended light curve file name
  nf = os.path.join(EVEREST_DAT, 'k2', 'c%02d' % int(campaign),
                    ('%09d' % EPIC)[:4] + '00000', 
                    ('%09d' % EPIC)[4:], model + '.npz')
  
  # Get the data
  data = np.load(nf)
  time = data['time']
  breakpoints = data['breakpoints']
  flux = data['fraw'] - data['model'] 
  mask = np.array(list(set(np.concatenate([data['outmask'], data['badmask'], 
                                          data['nanmask'], data['transitmask']]))), 
                                          dtype = int)
  flux[data['nanmask']] = np.nan
  return time, flux, breakpoints, mask
  
def Fit(EPIC, campaign = None, module = None, model = 'nPLD', **kwargs):
  '''
  Fit a given K2 target with the CBV design matrix and plot the results.
  
  '''
  
  # Get the data
  if campaign is None:
    campaign = Campaign(EPIC)
  if module is None:
    module = Module(EPIC)
  time, flux, breakpoints, mask = GetLightCurve(EPIC, campaign, **kwargs)
  path = os.path.join(EVEREST_DAT, 'k2', 'cbv', 'c%02d' % campaign, str(module), model)
  
  # Get the design matrix  
  X = GetX(campaign, module, model = model, **kwargs)
  
  # Loop over all the light curve segments
  model = [None for b in range(len(breakpoints))]
  weights = [None for b in range(len(breakpoints))]
  for b in range(len(breakpoints)):
    
    # Get the indices for this light curve segment
    inds = GetChunk(time, breakpoints, b)
    masked_inds = GetChunk(time, breakpoints, b, mask = mask)

    # Ordinary least squares
    mX = X[masked_inds]
    A = np.dot(mX.T, mX)
    B = np.dot(mX.T, flux[masked_inds])
    weights[b] = np.linalg.solve(A, B)
    model[b] = np.dot(X[inds], weights[b])
    
    # Vertical alignment
    if b == 0:
      model[b] -= np.nanmedian(model[b])
    else:
      # Match the first finite model point on either side of the break
      # We could consider something more elaborate in the future
      i0 = -1 - np.argmax([np.isfinite(model[b - 1][-i]) for i in range(1, len(model[b - 1]) - 1)])
      i1 = np.argmax([np.isfinite(model[b][i]) for i in range(len(model[b]))])
      model[b] += (model[b - 1][i0] - model[b][i1])
  
  # Join model and normalize  
  model = np.concatenate(model)
  model -= np.nanmedian(model)
  
  # Plot the light curve, model, and corrected light curve
  pl.switch_backed('Agg')
  fig, ax = pl.subplots(2, figsize = (12, 6), sharex = True)
  ax[0].plot(np.delete(time, mask), np.delete(flux, mask), 'k.', markersize = 3)
  ax[0].plot(np.delete(time, mask), np.delete(model, mask) + np.nanmedian(flux), 'r-')
  ax[1].plot(np.delete(time, mask), np.delete(flux - model, mask), 'k.', markersize = 3)
  ax[0].set_ylabel('Original Flux', fontsize = 16)
  ax[1].set_ylabel('Corrected Flux', fontsize = 16)
  ax[1].set_xlabel('Time (BJD - 2454833)', fontsize = 16)
  
  # Force the axis range based on the non-outliers
  ax[0].margins(0.01, None)
  ax[0].set_xlim(ax[0].get_xlim())
  ax[0].set_ylim(ax[0].get_ylim())
  ax[1].set_xlim(ax[1].get_xlim())
  ax[1].set_ylim(ax[1].get_ylim())
  
  # Now plot the outliers
  ax[0].plot(time[mask], flux[mask], 'r.', markersize = 3, alpha = 0.5)
  ax[1].plot(time[mask], flux[mask] - model[mask], 'r.', markersize = 3, alpha = 0.5)
  
  pl.suptitle('EPIC %d' % EPIC, fontsize = 20)
  fig.savefig(os.path.join(path, '%d.pdf' % EPIC))
  pl.close()