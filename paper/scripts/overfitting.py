#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Figure 4
--------

This script reproduces Figure 4 in the paper. It de-trends EPIC 201367065 with
third order PLD, keeping all 8435 components, and leading to gross overfitting.
The source code is available 
`here <https://github.com/rodluger/everest/blob/master/paper/scripts/overfitting.py>`_.

.. figure:: ../paper/tex/images/overfitting.png
    :width: 500px
    :align: center
    :height: 100px
    :alt: alternate text
    :figclass: align-center

    **Figure 4** *Top:* Third order PLD applied to EPIC 201367065, but this time 
    keeping all basis vectors. Compare to Figure 1. While the median scatter improved 
    by a factor of about 4, the scatter in the transits (which were masked during the 
    de-trending) increased by a factor of several thousand. *Bottom:* The same figure, 
    but zoomed out to show the in-transit scatter.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np
import os, sys
EVEREST_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IMG_PATH = os.path.join(EVEREST_ROOT, 'paper', 'tex', 'images')
sys.path.insert(1, EVEREST_ROOT)
import everest
from everest.data import GetK2Data
from everest.utils import RMS
from everest.compute import GetMasks
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import george
from itertools import combinations_with_replacement as multichoose

def PLDBasis(fpix, pld_order = 1, cross_terms = True):
  frac = fpix / np.sum(fpix, axis = 1).reshape(-1, 1)
  if not cross_terms:
    x = np.hstack([frac ** (n + 1) for n in range(pld_order)])
  else:
    x = np.empty(shape = (frac.shape[0], 0), dtype = 'float64')
    for n in range(1, pld_order + 1):  
      x = np.hstack([x, np.product(list(multichoose(frac.T, n)), axis = 1).T])
  return x

def PLDModel(C, X):
  return np.dot(C, X.T)

def PLDCoeffs(X, Y, time, errors, gp):
  gp.compute(time, errors)
  A = np.dot(X.T, gp.solver.apply_inverse(X))
  B = np.dot(X.T, gp.solver.apply_inverse(Y))
  C = np.linalg.solve(A, B)
  return C

def GetData(EPIC = 201367065):
  k2star = GetK2Data(EPIC, apnum = 15)
  time = k2star.time
  bkg = k2star.bkg
  bkgerr = k2star.bkgerr
  apidx = np.where(k2star.apertures[15] & 1 & ~np.isnan(k2star.fpix[0]))
  fpix = np.array([f[apidx] for f in k2star.fpix], dtype='float64') - \
                   bkg.reshape(bkg.shape[0], 1)
  perr = np.array([f[apidx] for f in k2star.perr], dtype='float64')
  perr = np.sqrt(perr ** 2 + bkgerr.reshape(bkgerr.shape[0], 1) ** 2)
  flux = np.sum(fpix, axis = 1)
  ferr = np.sqrt(np.sum(perr ** 2, axis = 1))
  sidx = np.argsort(np.nansum(fpix, axis = 0))[::-1]
  fpix = fpix[:, sidx]
  mask, _, _, _ = GetMasks(time, flux, fpix, ferr, 5., planets = k2star.planets)  
  M = lambda x: np.delete(x, mask, axis = 0)
  U = lambda x: x[mask]
  gp = george.GP(np.std(flux) ** 2 * george.kernels.Matern32Kernel(30. ** 2))
  return time, fpix, flux, ferr, M, U, gp

def GetDetrended(EPIC = 201367065):

  try:
    data = np.load(os.path.join('npz', 'overfitting.npz'))
    time = data['time']
    flux = data['flux']
    Y = data['Y']
    rms = data['rms']
    rms_in_transit = data['rms_in_transit']
    
  except:
  
    time, fpix, flux, ferr, mask, unmask, gp = GetData(EPIC)
    X = PLDBasis(fpix, pld_order = 3, cross_terms = True)
    C = PLDCoeffs(mask(X), mask(flux), mask(time), mask(ferr), gp)
    M = PLDModel(C, X)
    Y = flux - M
    f = (flux - M + np.median(flux)) / np.median(flux)
    rms = RMS(mask(f))
    rms_in_transit = RMS(unmask(f))
    
    time = time[:1767]
    flux = flux[:1767]
    Y = Y[:1767]
    
    Y -= np.median(Y)
    
    np.savez(os.path.join('npz', 'overfitting.npz'), time = time, flux = flux, 
             Y = Y, rms = rms, rms_in_transit = rms_in_transit)  

  return time, flux, Y, rms, rms_in_transit

if __name__ == '__main__':
  fig, ax = pl.subplots(2, figsize = (7, 4), sharex = True)
  time, flux, Y, rms, rms_in_transit = GetDetrended()

  ax[0].plot(time, (Y + np.median(flux)) / np.median(flux), 'k.', alpha = 0.3, label = '%.1f ppm' % rms)
  ax[0].set_ylabel('De-trended Flux', fontsize = 12)
  ax[0].text(2014.4, 1.00425, "6-hr scatter: %.1f ppm" % rms, ha="right", va="top", fontsize=12,
             bbox = dict(fc = 'w'))
  ax[0].yaxis.set_major_locator(MaxNLocator(5))
  ax[0].margins(0.01, None)
  ax[0].set_ylim(0.996, 1.0050)
  ax[0].ticklabel_format(useOffset=False)

  ax[1].plot(time, (Y + np.median(flux)) / np.median(flux), 'k.', alpha = 0.3, label = '%.1f ppm' % rms)
  ax[1].set_ylabel('De-trended Flux', fontsize = 12, labelpad = 20.5)
  ax[1].yaxis.set_major_locator(MaxNLocator(5))
  ax[1].margins(0.01, None)
  ax[1].set_ylim(0.86, 1.14)
  ax[1].ticklabel_format(useOffset=False)
  ax[1].set_xlabel('Time (days)', fontsize = 16)

  fig.savefig(os.path.join(IMG_PATH, 'overfitting.png'), bbox_inches = 'tight')
  fig.savefig(os.path.join(IMG_PATH, 'overfitting.pdf'), bbox_inches = 'tight')