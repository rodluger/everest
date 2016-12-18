#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
sysrem.py
---------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from ...config import EVEREST_DAT
from .aux import GetK2Campaign, Campaign, Module
import os
import numpy as np
import matplotlib.pyplot as pl
import george
from george.kernels import Matern32Kernel, WhiteKernel

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

def Channels(module):
  '''
  Returns the channels contained in the given K2 module.
  
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
    # TODO: Keep the original copy as well in the future
    y = np.interp(t, np.delete(t, m), np.delete(y, m))
    err = np.interp(t, np.delete(t, m), np.delete(err, m))
       
    # Append to our running lists
    fluxes.append(y)
    errors.append(err)
    kpars.append(data['kernel_params'])
    
  return time, breakpoints, np.array(fluxes), np.array(errors), kpars

def SysRem(time, flux, err, nrec = 3, niter = 10, kernels = None):
  '''
  
  '''
  
  nflx, tlen = flux.shape
  
  # Get normalized fluxes
  med = np.nanmedian(flux, axis = 1).reshape(-1, 1)
  y = flux - med
  
  # Compute the inverse of the variances
  invvar = 1. / err ** 2
  
  # Our linear model for each light curve
  model = np.zeros((nflx, tlen))

  # Recover `nrec` components
  for n in range(nrec):
    
    # Initialize the weights and regressors
    c = np.zeros(nflx)
    a = np.ones(tlen)
    ri = y * invvar
  
    # Perform `niter` iterations
    for i in range(niter):
    
      print("Running iteration %d/%d..." % (i + 1, niter))
    
      # Compute the `c` vector (the weights)
      if kernels is None:
        
        # This is the solution assuming a diagonal covariance 
        # matrix for each light curve
        c = np.dot(ri, a) / np.dot(invvar, a ** 2)
    
      else:
        
        # This is the solution allowing red noise in each light curve (slower)
        X = a.reshape(-1,1)
        for j in range(nflx):
        
          print("Fitting light curve %d/%d..." % (j + 1, nflx))
        
          gp = george.GP(kernels[j])
          gp.compute(time, err[j])
          A = np.dot(X.T, gp.solver.apply_inverse(X))
          B = np.dot(X.T, gp.solver.apply_inverse(y[j]))
          c[j] = np.linalg.solve(A, B)
    
      # Compute the `a` vector (the regressors)
      a = np.dot(c, ri) / np.dot(c ** 2, invvar)
  
    # The linear model for this step
    m = np.outer(c, a)
  
    # Add to running model
    model += m
  
    # Remove this component
    y -= m
  
  # Add the median back in
  y += med
  
  return y

def Test():
  '''
  
  '''
  
  # Input
  lcfile = os.path.join(EVEREST_DAT, 'k2', 'cbv', 'test.npz')
  if not os.path.exists(lcfile):
    time, breakpoints, fluxes, errors, kpars = GetStars(2, 18, model = 'nPLD')
    np.savez(lcfile, time = time, breakpoints = breakpoints, fluxes = fluxes, 
             errors = errors, kpars = kpars)
  else:
    data = np.load(lcfile)
    time = data['time']
    breakpoints = data['breakpoints']
    fluxes = data['fluxes']
    errors = data['errors']
    kpars = data['kpars']
  
  # Output
  outfile = os.path.join(EVEREST_DAT, 'k2', 'cbv', 'test_out.npz')
  if not os.path.exists(outfile):
  
    # Get the kernels
    kernels = [WhiteKernel(k[0] ** 2) + k[1] ** 2 * Matern32Kernel(k[2] ** 2) for k in kpars]
  
    # Let's just do the first segment for now
    inds = GetChunk(time, breakpoints, 0)
  
    # Get the new fluxes
    new_fluxes = SysRem(time[inds], fluxes[:,inds], errors[:,inds], kernels = kernels, nrec = 2, niter = 5)
  
    # Save
    np.savez(outfile, new_fluxes = new_fluxes)
  
  else:
    
    inds = GetChunk(time, breakpoints, 0)
    data = np.load(outfile)
    new_fluxes = data['new_fluxes']
    
  # Set up the plot
  fig, axes = pl.subplots(5, 5, figsize = (12, 8))
  fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.95, bottom = 0.05, wspace = 0.025, hspace = 0.025)
  for i, ax in enumerate(axes.flatten()):
    ax.set_xticklabels([])
    ax.set_yticklabels([])
  
    f = fluxes[i][inds] - np.nanmedian(fluxes[i][inds])
    ax.plot(time[inds], f, 'b-', alpha = 0.5)
    model = f - new_fluxes[i][inds]
    model -= np.nanmedian(model)
    ax.plot(time[inds], model, 'r-')
    
    fsort = f[np.argsort(f)]
    lo = fsort[int(0.05 * len(fsort))]
    hi = fsort[int(0.95 * len(fsort))]
    pad = (hi - lo) * 0.05
    ax.set_ylim(lo-pad,hi+pad)
    
  fig.savefig(os.path.join(EVEREST_DAT, 'k2', 'cbv', 'test_out.pdf'))