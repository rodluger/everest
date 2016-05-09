#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
precision.py
------------
       
      
'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
import numpy as np
import matplotlib.pyplot as pl
import re

# User
campaigns = [4] # [1,2,3,4,5,6]
save = False
fig_precision1 = False
fig_precision2 = False
fig_comp_k2sff = False
fig_comp_k2sc = False
fig_comp_k2varcat = False
fig_comp_k2sff_k2sc = True

# Get the raw **Kepler** rms
kep_star, kep_kepmag, _, kep_raw = np.loadtxt(os.path.join('CDPP', 'kepler.tsv'), unpack = True)

# Get the raw K2 and EVEREST rms for each campaign
k2_star = [[] for i in range(7)]
k2_kepmag = [[] for i in range(7)]
k2_raw = [[] for i in range(7)]
k2_ever = [[] for i in range(7)]
k2_phot = [[] for i in range(7)]
for c in range(7):
  try:
    k2_star[c], k2_kepmag[c], _, k2_raw[c], k2_phot[c] = np.loadtxt(os.path.join('CDPP', 'k2raw_C%02d.tsv' % c), unpack = True)
    s, _, k2_ever[c] = np.loadtxt(os.path.join('CDPP', 'everest_C%02d.tsv' % c), unpack = True)
  except:
    continue
  if not np.allclose(k2_star[c], s):
    raise Exception("Input tables misaligned!")
  k2_star[c] = [int(e) for e in k2_star[c]]

# Now the figures
if fig_precision1:
  '''
  
  '''
  
  x = np.concatenate([k2_kepmag[c] for c in campaigns])
  y = np.concatenate([k2_ever[c] for c in campaigns])
  
  # Set up the figure
  fig, ax = pl.subplots(1, figsize = (8,6))
  ax.plot(kep_kepmag, kep_raw, 'y.', alpha = 0.025)
  ax.plot(x, y, 'b.', alpha = 0.15)

  # Dummy points for legend
  ax.plot(-1, -1, 'yo', label = 'Raw Kepler')
  ax.plot(-1, -1, 'bo', label = 'K2 detrended with EVEREST')
  leg = ax.legend(loc = 'upper left', numpoints = 1, handletextpad = 0)
  
  # Tweaks
  ax.set_ylabel('6-hr CDPP (ppm)', fontsize = 18)
  ax.set_xlabel('Kepler Magnitude', fontsize = 18)
  ax.set_xlim(10,16)
  ax.set_ylim(0, 500)

  # Bin for quantiles
  bins = np.arange(10.5,16,0.5)
  b_kep = np.zeros_like(bins) * np.nan
  b_k2 = np.zeros_like(bins) * np.nan
  for b, bin in enumerate(bins):
    i = np.where((kep_kepmag >= bin - 0.5) & (kep_kepmag < bin + 0.5))[0]
    b_kep[b] = np.nanmedian(kep_raw[i])
    i = np.where((x >= bin - 0.5) & (x < bin + 0.5))[0]
    b_k2[b] = np.nanmedian(y[i])
    
  pl.errorbar(bins - 0.025, b_kep, 1.4826 * np.median(np.abs(b_kep - np.median(b_kep))), color = 'y', fmt = 'o', lw = 1, ecolor = 'k', capthick = 1, zorder = 99)
  pl.errorbar(bins + 0.025, b_k2, 1.4826 * np.median(np.abs(b_k2 - np.median(b_k2))), color = 'b', fmt = 'o', lw = 1, ecolor = 'k', capthick = 1, zorder = 99)
  
  if save:
    fig.savefig('../tex/images/precision1.png' % name, bbox_inches = 'tight')
    fig.savefig('../tex/images/precision1.pdf' % name, bbox_inches = 'tight')
    pl.close()
  
if fig_precision2:
  '''
  
  '''

  x = np.concatenate([k2_kepmag[c] for c in campaigns])
  r = np.concatenate([k2_raw[c] for c in campaigns])
  e = np.concatenate([k2_ever[c] for c in campaigns])
  p = np.concatenate([k2_phot[c] for c in campaigns])

  # Set up the figure
  fig, ax = pl.subplots(1, figsize = (8,6))
  ax.plot(x, r, 'r.', alpha = 0.025)
  ax.plot(x, e, 'b.', alpha = 0.025)
  

  # Plot the median values
  bins = np.arange(10.5,19.,0.5)
  b_raw = np.zeros_like(bins) * np.nan
  b_k2 = np.zeros_like(bins) * np.nan
  b_pht = np.zeros_like(bins) * np.nan
  for b, bin in enumerate(bins):
    i = np.where((x >= bin - 0.5) & (x < bin + 0.5))[0]
    b_raw[b] = np.nanmedian(r[i])
    b_pht[b] = np.nanmedian(p[i])
    b_k2[b] = np.nanmedian(e[i])
  
  # Let's fit a straight line to the photon limit curve
  phot_lin = 10 ** np.poly1d(np.polyfit(bins, np.log10(b_pht), 1))(bins)
  
  ax.plot(bins, b_raw, 'r--', label = 'Raw', lw = 2)
  ax.plot(bins, b_k2, ls = '-', color = 'b', label = 'EVEREST', lw = 2)
  ax.plot(bins, phot_lin, 'y--', label = 'Photon limit', lw = 2)
  ax.legend(loc = 'upper left')

  # Tweaks
  ax.set_ylabel('6-hr CDPP (ppm)', fontsize = 18)
  ax.set_xlabel('Kepler Magnitude', fontsize = 18)
  ax.set_yscale('log')
  ax.set_xlim(10,19)
  ax.set_ylim(10 ** 0.75, 1e4)
  
  if save:
    fig.savefig('../tex/images/precision2.png' % name, bbox_inches = 'tight')
    fig.savefig('../tex/images/precision2.pdf' % name, bbox_inches = 'tight')
    pl.close()

if fig_comp_k2sff:
  '''
  
  '''

  # Get the EVEREST rms for each campaign
  evr_epic = np.concatenate([k2_star[c] for c in campaigns])
  evr_cdpp = np.concatenate([k2_ever[c] for c in campaigns])
  k = np.concatenate([k2_kepmag[c] for c in campaigns])

  # Get the K2SFF rms for each campaign
  sff_epic = []
  sff_cdpp = []
  for c in campaigns:
    try:
      a, _, b = np.loadtxt(os.path.join('CDPP', 'k2sff_C%02d.tsv' % c), unpack = True)
      sff_epic.extend(a)
      sff_cdpp.extend(b)
    except:
      continue
  
  # Now align with EVEREST for the one-to-one comparison
  stars = list(set(evr_epic) & set(sff_epic))
  x = []
  y = []
  for star in stars:
    i = np.argmax(evr_epic == star)
    j = np.argmax(sff_epic == star)
    x.append(k[i])
    y.append(((evr_cdpp[i] - sff_cdpp[j]) / sff_cdpp[j]))
  x = np.array(x)
  y = np.array(y)
  
  # Plot the equivalent of Fig. 10 in Aigrain+16
  fig, ax = pl.subplots(1, figsize = (9,6))
  ax.plot(x, y, 'b.', alpha = 0.2)
  ax.set_ylim(-1,1)
  ax.set_xlim(11,19)
  ax.axhline(0, color = 'gray', lw = 2, zorder = -99, alpha = 0.5)
  ax.axhline(0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)
  ax.axhline(-0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)

  bins = np.arange(10,19.5,0.5)
  by = np.zeros_like(bins) * np.nan
  for b, bin in enumerate(bins):
    i = np.where((x >= bin - 0.5) & (x < bin + 0.5))[0]
    by[b] = np.median(y[i])
  pl.plot(bins, by, 'k-', lw = 2)

  # Tweaks
  ax.set_ylabel(r'$\frac{\mathrm{CDPP}_{\mathrm{EVEREST}} - \mathrm{CDPP}_{\mathrm{K2SFF}}}{\mathrm{CDPP}_{\mathrm{K2SFF}}}$', fontsize = 22)
  ax.set_xlabel('Kepler Magnitude', fontsize = 18)
  
  if save:
    fig.savefig('../tex/images/comparison_k2sff.png' % name, bbox_inches = 'tight')
    fig.savefig('../tex/images/comparison_k2sff.pdf' % name, bbox_inches = 'tight')
    pl.close()

if fig_comp_k2sc:
  '''
  
  '''

  # Get the EVEREST rms for each campaign
  evr_epic = np.concatenate([k2_star[c] for c in campaigns])
  evr_cdpp = np.concatenate([k2_ever[c] for c in campaigns])
  k = np.concatenate([k2_kepmag[c] for c in campaigns])

  # Get the K2SC rms for each campaign
  sc_epic = []
  sc_cdpp = []
  for c in campaigns:
    try:
      a, _, b = np.loadtxt(os.path.join('CDPP', 'k2sc_C%02d.tsv' % c), unpack = True)
      sc_epic.extend(a)
      sc_cdpp.extend(b)
    except:
      continue
  
  # Now align with EVEREST for the one-to-one comparison
  stars = list(set(evr_epic) & set(sc_epic))
  x = []
  y = []
  for star in stars:
    i = np.argmax(evr_epic == star)
    j = np.argmax(sc_epic == star)
    x.append(k[i])
    y.append(((evr_cdpp[i] - sc_cdpp[j]) / sc_cdpp[j]))
  x = np.array(x)
  y = np.array(y)
  
  # Plot the equivalent of Fig. 10 in Aigrain+16
  fig, ax = pl.subplots(1, figsize = (9,6))
  ax.plot(x, y, 'b.', alpha = 0.2)
  ax.set_ylim(-1,1)
  ax.set_xlim(11,19)
  ax.axhline(0, color = 'gray', lw = 2, zorder = -99, alpha = 0.5)
  ax.axhline(0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)
  ax.axhline(-0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)

  bins = np.arange(10,19.5,0.5)
  by = np.zeros_like(bins) * np.nan
  for b, bin in enumerate(bins):
    i = np.where((x >= bin - 0.5) & (x < bin + 0.5))[0]
    by[b] = np.median(y[i])
  pl.plot(bins, by, 'k-', lw = 2)

  # Tweaks
  ax.set_ylabel(r'$\frac{\mathrm{CDPP}_{\mathrm{EVEREST}} - \mathrm{CDPP}_{\mathrm{K2SC}}}{\mathrm{CDPP}_{\mathrm{K2SC}}}$', fontsize = 22)
  ax.set_xlabel('Kepler Magnitude', fontsize = 18)
  
  if save:
    fig.savefig('../tex/images/comparison_k2sc.png' % name, bbox_inches = 'tight')
    fig.savefig('../tex/images/comparison_k2sc.pdf' % name, bbox_inches = 'tight')
    pl.close()

if fig_comp_k2varcat:
  '''
  
  '''

  # Get the EVEREST rms for each campaign
  evr_epic = np.concatenate([k2_star[c] for c in campaigns])
  evr_cdpp = np.concatenate([k2_ever[c] for c in campaigns])
  k = np.concatenate([k2_kepmag[c] for c in campaigns])

  # Get the K2varcat rms for each campaign
  varcat_epic = []
  varcat_cdpp = []
  for c in campaigns:
    try:
      a, _, b = np.loadtxt(os.path.join('CDPP', 'k2varcat_C%02d.tsv' % c), unpack = True)
      varcat_epic.extend(a)
      varcat_cdpp.extend(b)
    except:
      continue
  
  # Now align with EVEREST for the one-to-one comparison
  stars = list(set(evr_epic) & set(varcat_epic))
  x = []
  y = []
  for star in stars:
    i = np.argmax(evr_epic == star)
    j = np.argmax(varcat_epic == star)
    x.append(k[i])
    y.append(((evr_cdpp[i] - varcat_cdpp[j]) / varcat_cdpp[j]))
  x = np.array(x)
  y = np.array(y)
  
  # Plot the equivalent of Fig. 10 in Aigrain+16
  fig, ax = pl.subplots(1, figsize = (9,6))
  ax.plot(x, y, 'b.', alpha = 0.2)
  ax.set_ylim(-1,1)
  ax.set_xlim(11,19)
  ax.axhline(0, color = 'gray', lw = 2, zorder = -99, alpha = 0.5)
  ax.axhline(0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)
  ax.axhline(-0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)

  bins = np.arange(10,19.5,0.5)
  by = np.zeros_like(bins) * np.nan
  for b, bin in enumerate(bins):
    i = np.where((x >= bin - 0.5) & (x < bin + 0.5))[0]
    by[b] = np.median(y[i])
  pl.plot(bins, by, 'k-', lw = 2)

  # Tweaks
  ax.set_ylabel(r'$\frac{\mathrm{CDPP}_{\mathrm{EVEREST}} - \mathrm{CDPP}_{\mathrm{K2VARCAT}}}{\mathrm{CDPP}_{\mathrm{K2VARCAT}}}$', fontsize = 22)
  ax.set_xlabel('Kepler Magnitude', fontsize = 18)
  
  if save:
    fig.savefig('../tex/images/comparison_k2varcat.png' % name, bbox_inches = 'tight')
    fig.savefig('../tex/images/comparison_k2varcat.pdf' % name, bbox_inches = 'tight')
    pl.close()

if fig_comp_k2sff_k2sc:
  '''
  
  '''

  # Get the kepmags for each campaign
  evr_epic = np.concatenate([k2_star[c] for c in campaigns])
  k = np.concatenate([k2_kepmag[c] for c in campaigns])

  # Get the K2SC rms for each campaign
  sc_epic = []
  sc_cdpp = []
  for c in campaigns:
    try:
      a, _, b = np.loadtxt(os.path.join('CDPP', 'k2sc_C%02d.tsv' % c), unpack = True)
      sc_epic.extend(a)
      sc_cdpp.extend(b)
    except:
      continue

  # Get the K2SFF rms for each campaign
  sff_epic = []
  sff_cdpp = []
  for c in campaigns:
    try:
      a, _, b = np.loadtxt(os.path.join('CDPP', 'k2sff_C%02d.tsv' % c), unpack = True)
      sff_epic.extend(a)
      sff_cdpp.extend(b)
    except:
      continue
  
  # Now align with EVEREST for the one-to-one comparison
  stars = list(set(evr_epic) & set(sff_epic) & set(sc_epic))
  x = []
  y = []
  for star in stars:
    h = np.argmax(sc_epic == star)
    i = np.argmax(evr_epic == star)
    j = np.argmax(sff_epic == star)
    x.append(k[i])
    y.append(((sc_cdpp[h] - sff_cdpp[j]) / sff_cdpp[j]))
  x = np.array(x)
  y = np.array(y)
  
  # Plot the equivalent of Fig. 10 in Aigrain+16
  fig, ax = pl.subplots(1, figsize = (9,6))
  ax.plot(x, y, 'b.', alpha = 0.2)
  ax.set_ylim(-1,1)
  ax.set_xlim(11,19)
  ax.axhline(0, color = 'gray', lw = 2, zorder = -99, alpha = 0.5)
  ax.axhline(0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)
  ax.axhline(-0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)

  bins = np.arange(10,19.5,0.5)
  by = np.zeros_like(bins) * np.nan
  for b, bin in enumerate(bins):
    i = np.where((x >= bin - 0.5) & (x < bin + 0.5))[0]
    by[b] = np.median(y[i])
  pl.plot(bins, by, 'k-', lw = 2)

  # Tweaks
  ax.set_ylabel(r'$\frac{\mathrm{CDPP}_{\mathrm{K2SC}} - \mathrm{CDPP}_{\mathrm{K2SFF}}}{\mathrm{CDPP}_{\mathrm{K2SFF}}}$', fontsize = 22)
  ax.set_xlabel('Kepler Magnitude', fontsize = 18)
  
  if save:
    fig.savefig('../tex/images/comparison_k2sff_k2sc.png' % name, bbox_inches = 'tight')
    fig.savefig('../tex/images/comparison_k2sff_k2sc.pdf' % name, bbox_inches = 'tight')
    pl.close()

if not save:
  pl.show()
  pl.close()