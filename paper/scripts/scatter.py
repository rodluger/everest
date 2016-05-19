#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Figure 5
--------

This script reproduces Figure 5 in the paper.
It computes the cross-validation plot (3rd order PLD, 1 breakpoint) for EPIC 201497682.
This may take up to half an hour to run, as I'm computing it on a very
fine grid and over many iterations.
The source code is available 
`here <https://github.com/rodluger/everest/blob/master/paper/scripts/scatter.py>`_.

.. figure:: ../paper/tex/images/scatter.png
    :width: 500px
    :align: center
    :height: 100px
    :alt: alternate text
    :figclass: align-center

    **Figure 5** De-trended light curve precision as a function of the number of principal 
    components for EPIC 201497682. The blue dots are the median 6-hr precision (in ppm) 
    of the unmasked sec- tions of the light curve (the training set); the red dots are 
    the median precision in 6-hr chunks that were masked during the de-trending step 
    (the validation set). Solid curves indicate our GP fit to the data points. Initially, 
    the scatter decreases in both cases as the number of components is increased. 
    However, above about 50 components, while the scatter in the training set continues 
    to decrease, the scatter in the validation set (where the model is extrapolated) 
    begins to grow. This is the signature of overfitting. We therefore choose 50 principal 
    components for the de-trending, yielding a precision of 55 ppm (versus about 70 ppm for 
    the K2SFF de-trended flux).

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np
import os, sys
EVEREST_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IMG_PATH = os.path.join(EVEREST_ROOT, 'paper', 'tex', 'images')
sys.path.insert(1, EVEREST_ROOT)
import everest
from everest.detrend import PLDBasis, PLDCoeffs, PLDModel, ComputeScatter, SliceX
from everest.gp import GetGP
from everest.utils import InitLog, Breakpoints, RMS, Mask
from everest.compute import GetMasks
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
  npc = np.arange(1, npctot, 10)
  pred = np.zeros_like(npc)
  real = np.zeros_like(npc)
  
  for i, n in enumerate(npc):
    log.info('%d/%d...' % (n, npctot))
    pred[i], real[i] = ComputeScatter(SliceX(X, n, npctot), flux, time, ferr, gp, 
                                      mask = mask, niter = 30, nmasks = 10)

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
  fig = pl.figure(figsize = (8, 6))
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
  pl.plot(npc_pred, fpred, 'r-', label = 'Validation set', lw = 2)

  # Computed scatter
  freal, _ = gp_scatter.predict(real - np.nanmedian(real), npc_pred)
  freal += np.nanmedian(real)
  pl.plot(npc_pred, freal, 'b-', label = 'Training set', lw = 2) 

  # The minimum predicted scatter.
  k = np.nanargmin(fpred)
  mps = [fpred[k], npc_pred[k]]
  
  # Mark the best value
  pl.axvline(mps[1], color = 'k', lw = 2, alpha = 0.5, ls = '--', label = 'Best model')
  pl.axhline(mps[0], color = 'k', lw = 2, alpha = 0.5, ls = '--')
  
  # Add labels       
  pl.ylabel('Scatter (ppm)', fontsize = 20)
  pl.xlabel('Number of Principal Components', fontsize = 20)
  pl.legend(loc = 'upper right', fontsize = 14)
  for label in pl.gca().get_xticklabels() + pl.gca().get_yticklabels():
    label.set_fontsize(14)
  
  fig.savefig(os.path.join(IMG_PATH, 'scatter.png'), bbox_inches = 'tight')
  fig.savefig(os.path.join(IMG_PATH, 'scatter.pdf'), bbox_inches = 'tight')

if __name__ == '__main__':
  PlotScatter()