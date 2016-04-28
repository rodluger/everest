#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
sinusoid.py
-----------
This example shows how PLD can fail for variable stars. Even though the variable
signal is astrophysical -- i.e., present in all pixels -- PLD will try to fit out
the variability by exploiting the white noise in each of the pixels. In other words,
by inflating the white noise of the fit, PLD is able to find linear combinations of
the fractional pixel fluxes that approximate the astrophysical variability. This
generally results in a poor fit with much higher noise than the original data.
In order to fix this, we use a GP to fit the variability along the time dimension. 
It works beautifully!

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
EVEREST_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IMG_PATH = os.path.join(EVEREST_ROOT, 'paper', 'tex', 'images')
sys.path.insert(1, EVEREST_ROOT)
import everest
from everest.detrend import PLDBasis, PLDCoeffs, PLDModel
from everest.utils import RMS, MedianFilter
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import kplr
import george
from george.kernels import Matern32Kernel, WhiteKernel

def GetKeplerData(koi, fnum):
  '''
  Retrieves data from a single Kepler TPF file.
  
  '''
  
  client = kplr.API()
  tpf = client.koi(koi).get_target_pixel_files(short_cadence = False)
  with tpf[fnum].open() as f:  
    ap = f[2].data
    idx = np.where(ap & 2)
    qdata = f[1].data
  time = np.array(qdata.field('TIME'), dtype='float64')
  nan_inds = list(np.where(np.isnan(time))[0])
  time = np.delete(time, nan_inds)
  flux = np.array(qdata.field('FLUX'), dtype='float64')
  flux = np.delete(flux, nan_inds, 0)
  fpix = np.array([f[idx] for f in flux], dtype='float64')
  perr = np.array([f[idx] for f in qdata.field('FLUX_ERR')], dtype='float64')
  perr = np.delete(perr, nan_inds, 0)
  fsum = np.sum(fpix, axis = 1)
  ferr = np.sum(perr**2, axis = 1)**0.5
  quality = qdata.field('QUALITY')
  quality = np.delete(quality, nan_inds)
  qual_inds = []
  for b in [1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17]:
    qual_inds += list(np.where(quality & 2**(b-1))[0])
  nan_inds = list(np.where(np.isnan(fsum))[0])
  bad_inds = np.array(sorted(list(set(qual_inds + nan_inds))))

  time = np.delete(time, bad_inds)
  fpix = np.delete(fpix, bad_inds, 0)
  perr = np.delete(perr, bad_inds, 0)

  return time, fpix, perr

def GetDetrendedData():
  '''
  
  '''
  
  if os.path.exists(os.path.join('npz', 'sinusoid.npz')):
    data = np.load(os.path.join('npz', 'sinusoid.npz'))
    F1, F2, F3, F4 = data['F']
    P1, P2, P3, P4 = data['P']
    G1, G2, G3, G4 = data['G']
    time = data['time']
    
  else:
    # Get the data for TPF #3 (Quarter 4) of KOI 1275.01
    time, fpix, perr = GetKeplerData(1275.01, 3)
    npts = time.shape[0]
    flux = np.sum(fpix, axis = 1)
    ferr = np.sqrt(np.sum(perr ** 2, axis = 1))

    # Now add a sinusoidal signal (period 25 days, amplitude of 0.5%)
    sin_signal = (1. + 0.005 * np.sin(2 * np.pi / 25. * time).reshape(npts, 1))
    fpix_sin = fpix * sin_signal
    flux_sin = np.sum(fpix_sin, axis = 1)

    # Get the transit masks
    transit1 = np.where((time >= 368.3) & (time <= 369.3))[0]
    transit2 = np.where((time >= 418.6) & (time <= 419.6))[0]
    mask = np.append(transit1, transit2)
    MASK = lambda x: np.delete(x, mask, axis = 0)

    # 1. Get our PLD model for the real data, no GP
    gp = george.GP(WhiteKernel(0.))
    X, _ = PLDBasis(fpix, pld_order = 1, 
                 cross_terms = False, max_components = 100)
    C = PLDCoeffs(X, flux, time, ferr, gp, mask)    
    F1 = flux
    P1 = PLDModel(C, X)

    # 2. Detrend the sinusoid, no GP
    gp = george.GP(WhiteKernel(0.))
    X, _ = PLDBasis(fpix_sin, pld_order = 1, 
                 cross_terms = False, max_components = 100)
    C = PLDCoeffs(X, flux_sin, time, ferr, gp, mask)    
    F2 = flux_sin
    P2 = PLDModel(C, X)

    # 3. Detrend the sinusoid, polynomial
    poly_order = 10
    gp = george.GP(WhiteKernel(0.))
    X, _ = PLDBasis(fpix_sin, pld_order = 1, 
                 cross_terms = False, max_components = 100)
    t = (time - time[0]) / (time[-1] - time[0])
    poly = np.vstack([t ** n for n in range(1, poly_order + 1)]).T
    X = np.hstack([X, poly])
    C = PLDCoeffs(X, flux_sin, time, ferr, gp, mask)    
    F3 = flux_sin
    P3 = PLDModel(C, X)

    # 4. Detrend the sinusoid with a GP
    gp = george.GP(100 ** 2 * Matern32Kernel(20. ** 2))
    X, _ = PLDBasis(fpix_sin, pld_order = 1, 
                 cross_terms = False, max_components = 100)
    C = PLDCoeffs(X, flux_sin, time, ferr, gp, mask)    
    F4 = flux_sin
    P4 = PLDModel(C, X)

    # Now we model out the rest of the variability in each of the
    # lightcurves with a GP
    gp = george.GP(100 ** 2 * Matern32Kernel(10.))
    gp.compute(MASK(time), MASK(ferr))
    G = []
    for i, y in enumerate([F1 - P1, F2 - P2, F3 - P3, F4 - P4]):
      med = np.median(MASK(y))
      mu, _ = gp.predict(MASK(y) - med, time)
      G.append(mu + med)
    G1, G2, G3, G4 = G

    np.savez(os.path.join('npz', 'sinusoid.npz'), 
             time = time,
             F = [F1, F2, F3, F4], 
             P = [P1, P2, P3, P4], 
             G = [G1, G2, G3, G4])

  return time, F1, F2, F3, F4, P1, P2, P3, P4, G1, G2, G3, G4

# Grab the detrended data
time, F1, F2, F3, F4, P1, P2, P3, P4, G1, G2, G3, G4 = GetDetrendedData()

# The folding function
per = 50.28
t0 = 368.77
fold = lambda t: (t - t0 - per / 2.) % per - per / 2.
fold2 = lambda t, y: np.array(sorted(zip(fold(t), y))).T

# Plot
fig = pl.figure(figsize = (20,11))
fig.subplots_adjust(wspace = 0.05, hspace = 0.125, left = 0.075, right = 0.95, top = 0.85)
ax = np.array([[pl.subplot2grid((40, 4), (0, j), colspan=1, rowspan=9) for j in range(4)],
               [pl.subplot2grid((40, 4), (10, j), colspan=1, rowspan=9) for j in range(4)],
               [pl.subplot2grid((40, 4), (20, j), colspan=1, rowspan=9) for j in range(4)],
               [pl.subplot2grid((40, 4), (32, j), colspan=1, rowspan=7) for j in range(4)]])

# 1a. Raw data
ax[0,0].plot(time, F1, 'k.', alpha = 0.3, label = 'SAP Flux')
ax[0,0].plot(time, P1, 'r-', label = 'PLD Model')

# 1b. PLD-detrended data
ax[1,0].plot(time, F1 - P1, 'k.', alpha = 0.3)
ax[2,0].plot(time, F1 - P1 - G1, 'k.', alpha = 0.3)
ax[3,0].plot(fold(time), F1 - P1 - G1, 'k.', alpha = 1)
t, y = fold2(time, MedianFilter(F1 - P1 - G1, 10))
ax[3,0].plot(t, y, 'r-', lw = 2)

# 2a. Raw data times sinusoid
ax[0,1].plot(time, F2, 'k.', alpha = 0.3, label = 'SAP Flux')
ax[0,1].plot(time, P2, 'r-', label = 'PLD Model')

# 2b. PLD-detrended sinusoidal data
ax[1,1].plot(time, F2 - P2, 'k.', alpha = 0.3)
ax[2,1].plot(time, F2 - P2 - G2, 'k.', alpha = 0.3)
ax[3,1].plot(fold(time), F2 - P2 - G2, 'k.', alpha = 1)
t, y = fold2(time, MedianFilter(F2 - P2 - G2, 10))
ax[3,1].plot(t, y, 'r-', lw = 2)

# 3a. Raw data times sinusoid
ax[0,2].plot(time, F3, 'k.', alpha = 0.3, label = 'SAP Flux')
ax[0,2].plot(time, P3, 'r-', label = 'PLD Model')

# 3b. PLD + GP-detrended sinusoidal data
ax[1,2].plot(time, F3 - P3, 'k.', alpha = 0.3)
ax[2,2].plot(time, F3 - P3 - G3, 'k.', alpha = 0.3)
ax[3,2].plot(fold(time), F3 - P3 - G3, 'k.', alpha = 1)
t, y = fold2(time, MedianFilter(F3 - P3 - G3, 10))
ax[3,2].plot(t, y, 'r-', lw = 2)

# 4a. Raw data times sinusoid
ax[0,3].plot(time, F4, 'k.', alpha = 0.3, label = 'SAP Flux')
ax[0,3].plot(time, P4, 'r-', label = 'PLD Model')

# 4b. PLD + GP-detrended sinusoidal data
ax[1,3].plot(time, F4 - P4, 'k.', alpha = 0.3)
ax[2,3].plot(time, F4 - P4 - G4, 'k.', alpha = 0.3)
ax[3,3].plot(fold(time), F4 - P4 - G4, 'k.', alpha = 1)
t, y = fold2(time, MedianFilter(F4 - P4 - G4, 10))
ax[3,3].plot(t, y, 'r-', lw = 2)

# Fix plot ranges
for n in [0,1,2]:
  ylim = ax[n,0].get_ylim() + ax[n,1].get_ylim() + ax[n,2].get_ylim() + ax[n,3].get_ylim()
  ax[n,0].set_ylim(min(ylim), max(ylim))
  ax[n,1].set_ylim(min(ylim), max(ylim))
  ax[n,2].set_ylim(min(ylim), max(ylim))
  ax[n,3].set_ylim(min(ylim), max(ylim))
ax[0,0].margins(0.01, None)
xlim = ax[0,0].get_xlim()
for axis in ax.flatten():
  axis.set_xlim(*xlim)
for axis in [ax[3,0], ax[3,1], ax[3,2], ax[3,3]]:
  axis.set_xlim(-0.75, 0.75)
  axis.set_ylim(-100, 100)

# Ticks and tick labels
for i in [0,1,2,3]:
  for j in [1,2,3]:
    ax[i,j].set_yticklabels([])
for i in [0,1]:
  for j in [0,1,2,3]:
    ax[i,j].set_xticklabels([])

# Print the RMS
ax[2,0].annotate("%.1f ppm" % RMS((F1 - P1 - G1 + np.median(F1)) / np.median(F1)), xy = (0.05,0.9), xycoords = 'axes fraction', horizontalalignment = 'left', verticalalignment = 'top', bbox = dict(fc = "w"))
ax[2,1].annotate("%.1f ppm" % RMS((F2 - P2 - G2 + np.median(F2)) / np.median(F2)), xy = (0.05,0.9), xycoords = 'axes fraction', horizontalalignment = 'left', verticalalignment = 'top', bbox = dict(fc = "w"))
ax[2,2].annotate("%.1f ppm" % RMS((F3 - P3 - G3 + np.median(F3)) / np.median(F3)), xy = (0.05,0.9), xycoords = 'axes fraction', horizontalalignment = 'left', verticalalignment = 'top', bbox = dict(fc = "w"))
ax[2,3].annotate("%.1f ppm" % RMS((F4 - P4 - G4 + np.median(F4)) / np.median(F4)), xy = (0.05,0.9), xycoords = 'axes fraction', horizontalalignment = 'left', verticalalignment = 'top', bbox = dict(fc = "w"))

# Labels
#pl.suptitle('KOI 1275.01 Q4', fontsize = 25)
ax[0,0].set_title(r'a. Original', fontsize = 24, y = 1.03)
ax[0,1].set_title(r'b. PLD', fontsize = 24, y = 1.03)
ax[0,2].set_title(r'c. PLD + Polynomial', fontsize = 24, y = 1.03)
ax[0,3].set_title(r'd. PLD + GP', fontsize = 24, y = 1.03)
ax[0,0].set_ylabel('SAP Flux', fontsize = 14)
ax[1,0].set_ylabel('PLD Residuals', fontsize = 14)
ax[2,0].set_ylabel('Final Residuals', fontsize = 14)
ax[3,0].set_ylabel('Folded Residuals', fontsize = 12)
ax[2,0].set_xlabel('Time (days)', fontsize = 14)
ax[2,1].set_xlabel('Time (days)', fontsize = 14)
ax[2,2].set_xlabel('Time (days)', fontsize = 14)
ax[2,3].set_xlabel('Time (days)', fontsize = 14)
ax[3,0].set_xlabel('Time from transit center (days)', fontsize = 12)
ax[3,1].set_xlabel('Time from transit center (days)', fontsize = 12)
ax[3,2].set_xlabel('Time from transit center (days)', fontsize = 12)
ax[3,3].set_xlabel('Time from transit center (days)', fontsize = 12)

# Legends
ax[0,0].legend(loc = 'upper left', fontsize = 10, numpoints = 3)
ax[0,1].legend(loc = 'upper left', fontsize = 10, numpoints = 3)
ax[0,2].legend(loc = 'upper left', fontsize = 10, numpoints = 3)
ax[0,3].legend(loc = 'upper left', fontsize = 10, numpoints = 3)

for n in range(4):
  ax[0,n].yaxis.set_major_locator(MaxNLocator(5))
  ax[1,n].yaxis.set_major_locator(MaxNLocator(5))
  ax[2,n].yaxis.set_major_locator(MaxNLocator(5))

# Save
fig.savefig(os.path.join(IMG_PATH, 'sinusoid.png'), bbox_inches = 'tight')
fig.savefig(os.path.join(IMG_PATH, 'sinusoid.pdf'), bbox_inches = 'tight')