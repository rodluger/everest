#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
third_order.py
--------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np
import os, sys
EVEREST_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IMG_PATH = os.path.join(EVEREST_ROOT, 'paper', 'tex', 'images')
sys.path.insert(1, EVEREST_ROOT)
import everest
from everest.data import GetK2Data
from everest.utils import GetMasks, RMS
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
  # Keep only the 10 brightest pixels
  fpix = fpix[:, :10]
  mask, _, _, _ = GetMasks(time, flux, fpix, ferr, 5., planets = k2star.planets)  
  M = lambda x: np.delete(x, mask, axis = 0)
  gp = george.GP(np.std(flux) ** 2 * george.kernels.Matern32Kernel(2. ** 2))
  return time, fpix, flux, ferr, M, gp

def GetDetrended(EPIC = 201367065):

  try:
    data = np.load(os.path.join('npz', 'third_order.npz'))
    time = data['time']
    flux = data['flux']
    Y1 = data['Y1']
    Y2 = data['Y2']
    Y3 = data['Y3']
    RMS0 = data['RMS0']
    RMS1 = data['RMS1']
    RMS2 = data['RMS2']
    RMS3 = data['RMS3']
    
  except:
  
    time, fpix, flux, ferr, mask, gp = GetData(EPIC)
    RMS0 = RMS(flux / np.median(flux))
  
    X = PLDBasis(fpix, pld_order = 1, cross_terms = True)
    C = PLDCoeffs(mask(X), mask(flux), mask(time), mask(ferr), gp)
    M = PLDModel(C, X)
    Y1 = flux - M
    f = (flux - M + np.median(flux)) / np.median(flux)
    RMS1 = RMS(f)

    X = PLDBasis(fpix, pld_order = 2, cross_terms = True)
    C = PLDCoeffs(mask(X), mask(flux), mask(time), mask(ferr), gp)
    M = PLDModel(C, X)
    Y2 = flux - M
    f = (flux - M + np.median(flux)) / np.median(flux)
    RMS2 = RMS(f)

    X = PLDBasis(fpix, pld_order = 3, cross_terms = True)
    C = PLDCoeffs(mask(X), mask(flux), mask(time), mask(ferr), gp)
    M = PLDModel(C, X)
    Y3 = flux - M
    f = (flux - M + np.median(flux)) / np.median(flux)
    RMS3 = RMS(f)

    time = time[:1767]
    flux = flux[:1767]
    Y1 = Y1[:1767]
    Y2 = Y2[:1767]
    Y3 = Y3[:1767]
    
    Y1 -= np.median(Y1)
    Y2 -= np.median(Y2)
    Y3 -= np.median(Y3)
    
    np.savez(os.path.join('npz', 'third_order.npz'), time = time, flux = flux, 
             Y1 = Y1, Y2 = Y2, Y3 = Y3, RMS0 = RMS0, RMS1 = RMS1, 
             RMS2 = RMS2, RMS3 = RMS3)  

  return time, flux, Y1, Y2, Y3, RMS0, RMS1, RMS2, RMS3
    
fig, ax = pl.subplots(4, figsize = (6, 8))
time, flux, Y1, Y2, Y3, RMS0, RMS1, RMS2, RMS3 = GetDetrended()

ax[0].plot(time, flux / np.median(flux), 'k.', alpha = 0.3, label = '%.1f ppm' % RMS0)
ax[1].plot(time, (Y1 + np.median(flux)) / np.median(flux), 'k.', alpha = 0.3, label = '%.1f ppm' % RMS1)
ax[2].plot(time, (Y2 + np.median(flux)) / np.median(flux), 'k.', alpha = 0.3, label = '%.1f ppm' % RMS2)
ax[3].plot(time, (Y3 + np.median(flux)) / np.median(flux), 'k.', alpha = 0.3, label = '%.1f ppm' % RMS3)

ax[0].set_ylabel('Raw Flux', fontsize = 14)
ax[1].set_ylabel('1st Order PLD', fontsize = 14)
ax[2].set_ylabel('2nd Order PLD', fontsize = 14)
ax[3].set_ylabel('3rd Order PLD', fontsize = 14)

for axis in ax:
  axis.legend(loc = 'upper right', numpoints = 3, fontsize = 10)
  axis.yaxis.set_major_locator(MaxNLocator(5))
  axis.margins(0.01, None)
  axis.set_ylim(0.996, 1.0050)
  axis.ticklabel_format(useOffset=False)
for axis in ax[:-1]:
  axis.set_xticklabels([])

ax[-1].set_xlabel('Time (days)', fontsize = 18)

fig.savefig(os.path.join(IMG_PATH, 'third_order.png'), bbox_inches = 'tight')
fig.savefig(os.path.join(IMG_PATH, 'third_order.pdf'), bbox_inches = 'tight')