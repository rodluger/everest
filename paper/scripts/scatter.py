#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
scatter.py
----------

Compute the scatter plot (3rd order PLD, 2 chunks) for EPIC 201497682.
This may take up to half an hour to run, as I'm computing it on a very
fine grid and over many iterations.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np
import os, sys
EVEREST_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IMG_PATH = os.path.join(EVEREST_ROOT, 'paper', 'tex', 'images')
sys.path.insert(1, EVEREST_ROOT)
import everest
from everest.detrend import PLDBasis, ComputeScatter, SliceX
from everest.gp import GetGP
from everest.utils import InitLog, GetMasks, Breakpoints
from everest.data import GetK2Data
from everest.kernels import KernelModels
import matplotlib.pyplot as pl
import george
import logging
log = logging.getLogger(__name__)
InitLog(screen_level = logging.DEBUG)

def GetData(EPIC = 201497682):

  # Grab the data
  log.info('Pre-processing the data...')
  k2star = GetK2Data(EPIC, apnum = 15)
  time = k2star.time
  bkg = k2star.bkg
  bkgerr = k2star.bkgerr

  # Calculate the total background-corrected flux in the aperture
  apidx = np.where(k2star.apertures[15] & 1 & ~np.isnan(k2star.fpix[0]))
  fpix = np.array([f[apidx] for f in k2star.fpix], dtype='float64') - \
                   bkg.reshape(bkg.shape[0], 1)
  perr = np.array([f[apidx] for f in k2star.perr], dtype='float64')
  perr = np.sqrt(perr ** 2 + bkgerr.reshape(bkgerr.shape[0], 1) ** 2)
  flux = np.sum(fpix, axis = 1)
  ferr = np.sqrt(np.sum(perr ** 2, axis = 1))
  npix_total = fpix.shape[1]

  # Sort the pixel fluxes in decreasing order of total photons
  sidx = np.argsort(np.nansum(fpix, axis = 0))[::-1]
  fpix = fpix[:, sidx]
  perr = perr[:, sidx]

  # Obtain transit and outlier masks
  log.info('Obtaining masks...')
  mask, _, _, _ = GetMasks(time, flux, fpix, ferr, 5, planets = k2star.planets)

  # Get the GP  
  log.info('Computing the GP...')                                                
  knum, kpars, white, _, kchisq, _ = GetGP(EPIC, time, fpix, ferr, mask = mask, niter = 2)['data']
  
  # Get the basis vectors
  breakpoints = Breakpoints(k2star.campaign, time, mask)
  X = PLDBasis(fpix, time = time, pld_order = 3, cross_terms = True, max_components = 500,
               breakpoints = breakpoints)
  
  # Save
  np.savez(os.path.join('npz', 'scatter_in.npz'), time = time, X = X, flux = flux, 
           ferr = ferr, mask = mask, knum = knum, kpars = kpars)

def GetScatter(EPIC = 201497682):

  # Grab the data
  try:
    data = np.load(os.path.join('npz', 'scatter_in.npz'))
  except:
    GetData(EPIC = EPIC)
    data = np.load(os.path.join('npz', 'scatter_in.npz'))
  time = data['time']
  X, npctot = data['X']
  flux = data['flux']
  ferr = data['ferr']
  mask = data['mask']
  knum = data['knum']
  kpars = data['kpars']
    
  # Get the scatter
  log.info('Computing the scatter...')
  kernel = KernelModels[knum]
  kernel[:] = kpars
  gp = george.GP(kernel.george_kernel())
  npc = np.arange(1, X.shape[1], 10)
  pred = np.zeros_like(npc)
  real = np.zeros_like(npc)
  for i, n in enumerate(npc):
    log.info('%d/%d...' % (n, X.shape[1]))
    pred[i], real[i] = ComputeScatter(SliceX(X, n, npctot), flux, time, ferr, gp, niter = 30, nmasks = 10)

  # Save
  np.savez(os.path.join('npz', 'scatter_out.npz'), npc = npc, pred = pred, real = real)

def PlotScatter(EPIC = 201497682):

  # Grab the data
  try:
    data = np.load(os.path.join('npz', 'scatter_out.npz'))
  except:
    GetScatter(EPIC = EPIC)
    data = np.load(os.path.join('npz', 'scatter_out.npz'))
  npc = data['npc']
  pred = data['pred']
  real = data['real']
  npc_pred = np.arange(1, npc[-1])
  sffrms = 70.5
  phtrms = 32.7
  
  # Plot it
  fig = pl.figure(figsize = (12, 8))
  pl.plot(npc, pred, 'r.', alpha = 0.25)
  pl.plot(npc, real, 'b.', alpha = 0.25)
  
  # Mark RMS
  pl.axhline(sffrms, color = 'k', ls = ':', alpha = 0.5, lw = 2)
  pl.annotate('K2SFF', xy = (400, sffrms + 1), xycoords = 'data', ha = 'center', alpha = 0.5)
  pl.axhline(phtrms, color = 'k', ls = ':', alpha = 0.5, lw = 2)
  pl.annotate('Photon limit', xy = (400, phtrms + 1), xycoords = 'data', ha = 'center', alpha = 0.5)
  
  # Fit the scatter with a GP. GP params are hard-coded for now.
  sig = 1.4826 * np.nanmedian(np.abs(pred - np.nanmedian(pred)))
  amp = 10 * sig
  tau = 100
  gp_scatter = george.GP(amp ** 2 * george.kernels.Matern32Kernel(tau ** 2))
  gp_scatter.compute(npc, np.ones_like(npc) * sig)
  
  # Predicted scatter
  fpred, _ = gp_scatter.predict(pred - np.nanmedian(pred), npc_pred)
  fpred += np.nanmedian(pred)
  pl.plot(npc_pred, fpred, 'r-', label = 'Masked')

  # Computed scatter
  freal, _ = gp_scatter.predict(real - np.nanmedian(real), npc_pred)
  freal += np.nanmedian(real)
  pl.plot(npc_pred, freal, 'b-', label = 'Unmasked') 

  # The minimum predicted scatter. Here we minimize the sum of the predicted
  # scatter and the difference between the predicted and computed scatter.
  # This isn't super rigorous, but we want both quantities to be small, so
  # this should further prevent us from choosing an outlier in fP1 as the minimum.
  met = fpred + np.abs(fpred - freal)
  k = np.nanargmin(met)
  mps = [met[k], npc_pred[k]]
  
  # Mark the best value
  pl.axvline(mps[1], color = 'k', lw = 2, alpha = 0.5, ls = '--', label = 'Best Model')
  pl.axhline(mps[0], color = 'k', lw = 2, alpha = 0.5, ls = '--')
  
  # Add labels       
  pl.ylabel('Scatter (ppm)', fontsize = 20)
  pl.xlabel('Number of Principal Components', fontsize = 20)
  pl.legend(loc = 'upper right', fontsize = 14)
  for label in pl.gca().get_xticklabels() + pl.gca().get_yticklabels():
    label.set_fontsize(14)
  
  fig.savefig(os.path.join(IMG_PATH, 'scatter.png'), bbox_inches = 'tight')
  fig.savefig(os.path.join(IMG_PATH, 'scatter.pdf'), bbox_inches = 'tight')

PlotScatter()