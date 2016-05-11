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
from matplotlib.ticker import MaxNLocator
import re

# The campaigns we'll use for the synthesis plots
campaigns = [0,1,2,3,4,5]
figures = [1,2,3,4,5,6]

# Get the raw **Kepler** rms
kep_star, kep_kepmag, _, kep_raw = np.loadtxt(os.path.join('CDPP', 'kepler.tsv'), unpack = True)

# Get the raw K2 and EVEREST rms for each campaign
k2_star = [[] for i in range(8)]
k2_kepmag = [[] for i in range(8)]
k2_raw = [[] for i in range(8)]
k2_ever = [[] for i in range(8)]
k2_phot = [[] for i in range(8)]
for c in range(8):
  try:
    k2_star[c], k2_kepmag[c], _, k2_raw[c], k2_phot[c] = np.loadtxt(os.path.join('CDPP', 'k2raw_C%02d.tsv' % c), unpack = True)
    s, _, k2_ever[c] = np.loadtxt(os.path.join('CDPP', 'everest_C%02d.tsv' % c), unpack = True)
  except:
    continue
  if not np.allclose(k2_star[c], s):
    raise Exception("Input tables misaligned!")
  k2_star[c] = [int(e) for e in k2_star[c]]

# HACK: Campaign 0 magnitudes are reported to the nearest tenth
# We'll add some artificial scatter to make it plot nicer
k2_kepmag[0] = np.array(k2_kepmag[0]) + 0.1 * np.random.rand(len(k2_star[0])) - 0.05

# Define the figure functions
def fig_comp_kepler(fig, ax, campaigns, errorbars = True, labels = True):
  '''
  
  '''
  
  x = np.concatenate([k2_kepmag[c] for c in campaigns])
  y = np.concatenate([k2_ever[c] for c in campaigns])
  
  ax.plot(kep_kepmag, kep_raw, 'y.', alpha = 0.025)
  ax.plot(x, y, 'b.', alpha = 0.1)

  if labels:
    # Dummy points for legend
    ax.plot(-1, -1, 'yo', label = 'Raw Kepler')
    ax.plot(-1, -1, 'bo', label = 'K2 detrended with EVEREST')
    leg = ax.legend(loc = 'upper left', numpoints = 1, handletextpad = 0)
  
    # Labels
    ax.set_ylabel('6-hr CDPP (ppm)', fontsize = 18)
    ax.set_xlabel('Kepler Magnitude', fontsize = 18)
  
  # Tweaks
  ax.set_xlim(10,16)
  ax.set_ylim(0, 500)

  # Bin for quantiles
  if errorbars:
    bins = np.arange(10.5,16,0.5)
    b_kep = np.zeros_like(bins) * np.nan
    b_k2 = np.zeros_like(bins) * np.nan
    for b, bin in enumerate(bins):
      i = np.where((kep_kepmag >= bin - 0.5) & (kep_kepmag < bin + 0.5))[0]
      b_kep[b] = np.nanmedian(kep_raw[i])
      i = np.where((x >= bin - 0.5) & (x < bin + 0.5))[0]
      b_k2[b] = np.nanmedian(y[i])
    
    ax.errorbar(bins - 0.025, b_kep, 1.4826 * np.median(np.abs(b_kep - np.median(b_kep))), color = 'y', fmt = 'o', lw = 1, ecolor = 'k', capthick = 1, zorder = 99)
    ax.errorbar(bins + 0.025, b_k2, 1.4826 * np.median(np.abs(b_k2 - np.median(b_k2))), color = 'b', fmt = 'o', lw = 1, ecolor = 'k', capthick = 1, zorder = 99)

def fig_precision(fig, ax, campaigns, labels = True):
  '''
  
  '''

  x = np.concatenate([k2_kepmag[c] for c in campaigns])
  r = np.concatenate([k2_raw[c] for c in campaigns])
  e = np.concatenate([k2_ever[c] for c in campaigns])
  p = np.concatenate([k2_phot[c] for c in campaigns])

  # Set up the figure
  ax.plot(x, r, 'r.', alpha = 0.025)
  ax.plot(x, e, 'b.', alpha = 0.025)
  
  # Plot the median values
  bins = np.arange(10.5,19.,0.5)
  b_raw = np.zeros_like(bins) * np.nan
  b_k2 = np.zeros_like(bins) * np.nan
  mode = np.zeros_like(bins) * np.nan
  b_pht = np.zeros_like(bins) * np.nan
  for b, bin in enumerate(bins):
    i = np.where((x >= bin - 0.5) & (x < bin + 0.5))[0]
    b_raw[b] = np.nanmedian(r[i])
    b_pht[b] = np.nanmedian(p[i])
    b_k2[b] = np.nanmedian(e[i])
    
    # Compute the mode?
    h, edg = np.histogram(np.log10(e[i]), range = (np.log10(b_pht[b]), np.log10(b_raw[b])), bins = 100)
    j = np.argmax(h)
    mode[b] = 10 ** edg[j]
    
  # Let's fit a straight line to the photon limit curve
  phot_lin = 10 ** np.poly1d(np.polyfit(bins, np.log10(b_pht), 1))(bins)
  
  ax.plot(bins, b_raw, 'r--', label = 'Raw', lw = 2)
  ax.plot(bins, b_k2, ls = '-', color = 'b', label = 'EVEREST', lw = 2)
  ax.plot(bins, phot_lin, 'y--', label = 'Photon limit', lw = 2)
  
  # The mode... Perhaps someday.
  #ax.plot(bins, mode, ls = '--', color = 'b', label = 'EVEREST (mode)', lw = 2)
  
  if labels:
    ax.legend(loc = 'upper left')
    ax.set_ylabel('6-hr CDPP (ppm)', fontsize = 18)
    ax.set_xlabel('Kepler Magnitude', fontsize = 18)
  
  try:
    ax.set_yscale('log')
  except:
    pass
  ax.set_xlim(10,19)
  ax.set_ylim(10 ** 0.75, 1e4)

def fig_comp_k2sff(fig, ax, campaigns, labels = True):
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
  
  ax.plot(x, y, 'b.', alpha = 0.1)
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
  ax.plot(bins, by, 'k-', lw = 2)

  if labels:
    ax.set_ylabel(r'$\frac{\mathrm{CDPP}_{\mathrm{EVEREST}} - \mathrm{CDPP}_{\mathrm{K2SFF}}}{\mathrm{CDPP}_{\mathrm{K2SFF}}}$', fontsize = 22)
    ax.set_xlabel('Kepler Magnitude', fontsize = 18)

def fig_comp_k2sc(fig, ax, campaigns, labels = True):
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
  ax.plot(x, y, 'b.', alpha = 0.1)
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
  ax.plot(bins, by, 'k-', lw = 2)

  if labels:
    ax.set_ylabel(r'$\frac{\mathrm{CDPP}_{\mathrm{EVEREST}} - \mathrm{CDPP}_{\mathrm{K2SC}}}{\mathrm{CDPP}_{\mathrm{K2SC}}}$', fontsize = 22)
    ax.set_xlabel('Kepler Magnitude', fontsize = 18)

def fig_comp_k2varcat(fig, ax, campaigns, labels = True):
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
  ax.plot(x, y, 'b.', alpha = 0.1)
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
  ax.plot(bins, by, 'k-', lw = 2)

  if labels:
    ax.set_ylabel(r'$\frac{\mathrm{CDPP}_{\mathrm{EVEREST}} - \mathrm{CDPP}_{\mathrm{K2VARCAT}}}{\mathrm{CDPP}_{\mathrm{K2VARCAT}}}$', fontsize = 22)
    ax.set_xlabel('Kepler Magnitude', fontsize = 18)

def fig_comp_k2sff_k2sc(fig, ax, campaigns, labels = True):
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
  
  # Plot Fig. 10 in Aigrain+16
  ax.plot(x, y, 'b.', alpha = 0.1)
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
  ax.plot(bins, by, 'k-', lw = 2)

  if labels:
    ax.set_ylabel(r'$\frac{\mathrm{CDPP}_{\mathrm{K2SC}} - \mathrm{CDPP}_{\mathrm{K2SFF}}}{\mathrm{CDPP}_{\mathrm{K2SFF}}}$', fontsize = 22)
    ax.set_xlabel('Kepler Magnitude', fontsize = 18)

# 1. Comparison to raw Kepler
# ---------------------------
if 1 in figures:
  fig, ax = pl.subplots(1, figsize = (8,6))
  fig_comp_kepler(fig, ax, campaigns = campaigns)
  fig.savefig('../tex/images/comparison_kepler.png', bbox_inches = 'tight')
  pl.close()
  fig, ax = pl.subplots(2, 3, figsize = (14, 6))
  fig.subplots_adjust(wspace = 0.05, hspace = 0.075)
  ax = ax.flatten()
  for n in range(6):
    fig_comp_kepler(fig, ax[n], campaigns = [n], errorbars = False, labels = False)
    if n not in [3,4,5]: 
      ax[n].set_xticklabels([])
    else:
      ax[n].xaxis.set_major_locator(MaxNLocator(prune='upper', integer=True))
    if n not in [0,3]: 
      ax[n].set_yticklabels([])
    ax[n].annotate('C%02d' % n, xy = (0.02, 0.96), xycoords = 'axes fraction', 
                   ha = 'left', va = 'top', fontsize = 14)
  fig.text(0.5, 0.015, 'Kepler Magnitude', ha='center', va='center', fontsize = 22)
  fig.text(0.08, 0.5, '6-hr CDPP (ppm)', fontsize = 22, 
           ha='center', va='center', rotation='vertical')
  fig.savefig('../tex/images/comparison_kepler_by_campaign.png', bbox_inches = 'tight')
  pl.close()

# 2. Overall CDPP
# ---------------
if 2 in figures:
  fig, ax = pl.subplots(1, figsize = (8,6))
  fig_precision(fig, ax, campaigns)
  fig.savefig('../tex/images/precision.png', bbox_inches = 'tight')
  pl.close()
  fig, ax = pl.subplots(2, 3, figsize = (14, 6))
  fig.subplots_adjust(wspace = 0.05, hspace = 0.075)
  ax = ax.flatten()
  for n in range(6):
    fig_precision(fig, ax[n], campaigns = [n], labels = False)
    if n not in [3,4,5]: 
      ax[n].set_xticklabels([])
    else:
      ax[n].xaxis.set_major_locator(MaxNLocator(prune='upper', integer=True))
    if n not in [0,3]: 
      ax[n].set_yticklabels([])
    ax[n].annotate('C%02d' % n, xy = (0.02, 0.96), xycoords = 'axes fraction', 
                   ha = 'left', va = 'top', fontsize = 14)
  fig.text(0.5, 0.015, 'Kepler Magnitude', ha='center', va='center', fontsize = 22)
  fig.text(0.08, 0.5, '6-hr CDPP (ppm)', fontsize = 22, 
           ha='center', va='center', rotation='vertical')
  fig.savefig('../tex/images/precision_by_campaign.png', bbox_inches = 'tight')
  pl.close()

# 3. Comparison to K2SFF
# ----------------------
if 3 in figures:
  fig, ax = pl.subplots(1, figsize = (9,6))
  fig_comp_k2sff(fig, ax, campaigns)
  fig.savefig('../tex/images/comparison_k2sff.png', bbox_inches = 'tight')
  pl.close()
  
  fig, ax = pl.subplots(2, 3, figsize = (12, 8))
  fig.subplots_adjust(wspace = 0.05, hspace = 0.075, left = 0.15)
  ax = ax.flatten()
  for n in range(6):
    fig_comp_k2sff(fig, ax[n], campaigns = [n], labels = False)
    if n not in [3,4,5]: 
      ax[n].set_xticklabels([])
    else:
      ax[n].xaxis.set_major_locator(MaxNLocator(prune='upper', integer=True))
    if n not in [0,3]: 
      ax[n].set_yticklabels([])
    ax[n].annotate('C%02d' % n, xy = (0.02, 0.96), xycoords = 'axes fraction', 
                   ha = 'left', va = 'top', fontsize = 14)
  fig.text(0.5, 0.015, 'Kepler Magnitude', ha='center', va='center', fontsize = 24)
  fig.text(0.08, 0.5, r'$\frac{\mathrm{CDPP}_{\mathrm{EVEREST}} - \mathrm{CDPP}_{\mathrm{K2SFF}}}{\mathrm{CDPP}_{\mathrm{K2SFF}}}$', fontsize = 28, 
           ha='center', va='center', rotation='vertical')
  fig.savefig('../tex/images/comparison_k2sff_by_campaign.png', bbox_inches = 'tight')
  pl.close()

# 4. Comparison to K2SC
# ---------------------
if 4 in figures:
  fig, ax = pl.subplots(1, figsize = (9,6))
  fig_comp_k2sc(fig, ax, campaigns)
  fig.savefig('../tex/images/comparison_k2sc.png', bbox_inches = 'tight')
  pl.close()
  
  fig, ax = pl.subplots(1, 2, figsize = (12, 3.5))
  fig.subplots_adjust(wspace = 0.05, hspace = 0.075, bottom = 0.15)
  ax = ax.flatten()
  fig_comp_k2sff(fig, ax[0], campaigns = [4], labels = False)
  fig_comp_k2sff(fig, ax[1], campaigns = [5], labels = False)
  ax[0].xaxis.set_major_locator(MaxNLocator(prune='upper', integer=True))
  ax[1].xaxis.set_major_locator(MaxNLocator(prune='upper', integer=True))
  ax[1].set_yticklabels([])
  ax[0].annotate('C04', xy = (0.02, 0.96), xycoords = 'axes fraction', ha = 'left', va = 'top', fontsize = 14)
  ax[1].annotate('C05', xy = (0.02, 0.96), xycoords = 'axes fraction', ha = 'left', va = 'top', fontsize = 14)
  fig.text(0.5, 0.015, 'Kepler Magnitude', ha='center', va='center', fontsize = 18)
  ax[0].set_ylabel(r'$\frac{\mathrm{CDPP}_{\mathrm{EVEREST}} - \mathrm{CDPP}_{\mathrm{K2SC}}}{\mathrm{CDPP}_{\mathrm{K2SC}}}$', fontsize = 18)
  fig.savefig('../tex/images/comparison_k2sc_by_campaign.png', bbox_inches = 'tight')
  pl.close()

# 5. Comparison to K2VARCAT
# -------------------------
if 5 in figures:
  fig, ax = pl.subplots(1, figsize = (9,6))
  fig_comp_k2varcat(fig, ax, campaigns)
  fig.savefig('../tex/images/comparison_k2varcat.png', bbox_inches = 'tight')
  pl.close()
  
  fig, ax = pl.subplots(2, 3, figsize = (12, 8))
  fig.subplots_adjust(wspace = 0.05, hspace = 0.075, left = 0.15)
  ax = ax.flatten()
  for n in range(6):
    fig_comp_k2varcat(fig, ax[n], campaigns = [n], labels = False)
    if n not in [3,4,5]: 
      ax[n].set_xticklabels([])
    else:
      ax[n].xaxis.set_major_locator(MaxNLocator(prune='upper', integer=True))
    if n not in [0,3]: 
      ax[n].set_yticklabels([])
    ax[n].annotate('C%02d' % n, xy = (0.02, 0.96), xycoords = 'axes fraction', 
                   ha = 'left', va = 'top', fontsize = 14)
  fig.text(0.5, 0.015, 'Kepler Magnitude', ha='center', va='center', fontsize = 24)
  fig.text(0.08, 0.5, r'$\frac{\mathrm{CDPP}_{\mathrm{EVEREST}} - \mathrm{CDPP}_{\mathrm{K2VARCAT}}}{\mathrm{CDPP}_{\mathrm{K2VARCAT}}}$', fontsize = 28, 
           ha='center', va='center', rotation='vertical')
  fig.savefig('../tex/images/comparison_k2varcat_by_campaign.png', bbox_inches = 'tight')
  pl.close()

# 6. Compare K2SC to K2SFF
# ------------------------
if 6 in figures:
  fig, ax = pl.subplots(1, figsize = (9,6))
  fig_comp_k2sff_k2sc(fig, ax, campaigns)
  fig.savefig('../tex/images/comparison_k2sff_k2sc.png', bbox_inches = 'tight')
  pl.close()
  
  fig, ax = pl.subplots(1, 2, figsize = (12, 3.5))
  fig.subplots_adjust(wspace = 0.05, hspace = 0.075, bottom = 0.15)
  ax = ax.flatten()
  fig_comp_k2sff_k2sc(fig, ax[0], campaigns = [4], labels = False)
  fig_comp_k2sff_k2sc(fig, ax[1], campaigns = [5], labels = False)
  ax[0].xaxis.set_major_locator(MaxNLocator(prune='upper', integer=True))
  ax[1].xaxis.set_major_locator(MaxNLocator(prune='upper', integer=True))
  ax[1].set_yticklabels([])
  ax[0].annotate('C04', xy = (0.02, 0.96), xycoords = 'axes fraction', ha = 'left', va = 'top', fontsize = 14)
  ax[1].annotate('C05', xy = (0.02, 0.96), xycoords = 'axes fraction', ha = 'left', va = 'top', fontsize = 14)
  fig.text(0.5, 0.015, 'Kepler Magnitude', ha='center', va='center', fontsize = 18)
  ax[0].set_ylabel(r'$\frac{\mathrm{CDPP}_{\mathrm{K2SC}} - \mathrm{CDPP}_{\mathrm{K2SFF}}}{\mathrm{CDPP}_{\mathrm{K2SFF}}}$', fontsize = 18)
  fig.savefig('../tex/images/comparison_k2sff_k2sc_by_campaign.png', bbox_inches = 'tight')
  pl.close()