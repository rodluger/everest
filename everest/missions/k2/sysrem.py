#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
sysrem.py
---------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from ...config import EVEREST_DAT
from ...utils import InitLog
from .aux import GetK2Campaign, Campaign, Channels
import os
import numpy as np
import matplotlib.pyplot as pl
import george
from george.kernels import Matern32Kernel, WhiteKernel
from scipy.signal import savgol_filter
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
  channels = Channels(module)
  
  # Get the EPIC numbers
  all = GetK2Campaign(campaign)
  stars = np.array([s[0] for s in all if s[2] in channels and 
          os.path.exists(
          os.path.join(EVEREST_DAT, 'k2', 'c%02d' % int(campaign),
          ('%09d' % s[0])[:4] + '00000', 
          ('%09d' % s[0])[4:], model + '.npz'))], dtype = int)
  N = len(stars)
  assert N > 0, "No light curves found for campaign %d, module %d." % (campaign, module)

  # Loop over all stars and store the fluxes in a list
  fluxes = []
  errors = []
  kpars = []
  
  for n in range(N):

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

def SysRem(time, flux, err, nrec = 5, niter = 50, sv_win = 99, sv_order = 2, **kwargs):
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
    
    # Save this regressor after smoothing it a bit
    cbvs[n] = savgol_filter(a - np.nanmedian(a), sv_win, sv_order)
    
  return cbvs

def GetCBVs(campaign, module = None, model = 'nPLD', clobber = False, **kwargs):
  '''
  
  '''
  
  # Initialize logging?
  if len(logging.getLogger().handlers) == 0:
    InitLog(file_name = None, screen_level = logging.DEBUG)
  
  # All modules?
  if module is None:
  
    # We're going to plot the CBVs on the CCD
    fig = [None for n in range(1, kwargs.get('nrec', 5) + 1)]
    ax = [None for n in range(1, kwargs.get('nrec', 5) + 1)]
    for n in range(1, kwargs.get('nrec', 5) + 1):
      fig[n], ax[n] = pl.subplots(5, 5, figsize = (9, 9))
      fig[n].subplots_adjust(wspace = 0.025, hspace = 0.025)
      ax[n] = [None] + list(ax[n].flatten())
      for axis in [ax[n][1], ax[n][5], ax[n][21], ax[n][25]]:
        axis.set_visible(False)
      for i in range(1, 25):
        ax[n][i].set_xticks([])
        ax[n][i].set_yticks([])
        ax[n][i].annotate('%02d' % i, (0.5, 0.5), 
                      va = 'center', ha = 'center',
                      color = 'k', fontsize = 60, alpha = 0.05)
    
    # Get the CBVs
    for module in range(2, 25):
      X = GetCBVs(campaign, module = module, model = model, clobber = clobber, **kwargs)
      if X is not None:
        for n in range(1, min(5, X.shape[1])):
          ax[n][module].plot(X[:,n])
    
    for n in range(1, kwargs.get('nrec', 5) + 1):
      figname = os.path.join(EVEREST_DAT, 'k2', 'cbv', 'c%02d' % campaign, model + '_%02d.pdf' % n)
      fig[n].savefig(figname)
      pl.close(fig[n])
    
    return
  
  log.info('Computing CBVs for campaign %d, module %d...' % (campaign, module))
    
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
      try:
        time, breakpoints, fluxes, errors, kpars = GetStars(campaign, module, model = model, **kwargs)
      except:
        np.savez(lcfile, time = None, breakpoints = None, fluxes = None, errors = None, kpars = None)
        np.savez(xfile, X = None)
        return None
      np.savez(lcfile, time = time, breakpoints = breakpoints, fluxes = fluxes, errors = errors, kpars = kpars)
    else:
      lcs = np.load(lcfile)
      time = lcs['time']
      breakpoints = lcs['breakpoints']
      fluxes = lcs['fluxes']
      errors = lcs['errors']
      kpars = lcs['kpars']
    
    # Compute the design matrix  
    log.info('Running SysRem...')
    X = np.ones((len(time), 1 + kwargs.get('nrec', 5)))
    
    # Loop over the segments
    new_fluxes = np.zeros_like(fluxes)
    for b in range(len(breakpoints)):
      
      # Get the current segment's indices
      inds = GetChunk(time, breakpoints, b)
    
      # Update the error arrays with the white GP component
      for j in range(len(errors)):
        errors[j] = np.sqrt(errors[j] ** 2 + kpars[j][0] ** 2)
    
      # Get de-trended fluxes
      X[inds,1:] = SysRem(time[inds], fluxes[:,inds], errors[:,inds], **kwargs).T
      
    # Save
    np.savez(xfile, X = X)
  
  else:
    
    # Load from disk
    X = np.load(xfile)['X'][()]
  
  if X is not None:    
    return X[:,:kwargs.get('nrec', 5) + 1]
  else:
    return None