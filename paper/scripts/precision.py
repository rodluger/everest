#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
precision.py
------------
       
TODO: Work on this!     
        
'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
import numpy as np
import matplotlib.pyplot as pl
import re

# Get the raw K2 rms
epic, kepmag, raw, raw2, phot = np.loadtxt(os.path.join('CDPP', '_k2raw.tsv'), unpack = True)
epic = [int(e) for e in epic]

# Get the raw **Kepler** rms
kic, kep_kepmag, kep_raw, kep_raw2 = np.loadtxt(os.path.join('CDPP', 'kepler.tsv'), unpack = True)

# Get the EVEREST detrended data
anc_epic, anc, _, anc2, _ = np.loadtxt('/Users/rodrigo/src/anchor/anchorrms.csv', unpack = True)
anc_epic = np.array(anc_epic, dtype = int)
anc_kepmag = np.zeros(len(anc_epic)) * np.nan
for i, e in enumerate(anc_epic):
  where = np.argmax(e == epic)
  if where:
    anc_kepmag[i] = kepmag[where]

fig_precision1 = True


if fig_precision1:
  
  # Set up the figure
  fig, ax = pl.subplots(1, figsize = (8,6))
  ax.plot(kep_kepmag, y1, 'y.', alpha = 0.025)
  ax.plot(anc_kepmag, y2, 'b.', alpha = 0.15)

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
  b_anc = np.zeros_like(bins) * np.nan
  for b, bin in enumerate(bins):
    i = np.where((kep_kepmag >= bin - 0.5) & (kep_kepmag < bin + 0.5))[0]
    b_kep[b] = np.nanmedian(kep_raw[i])
    i = np.where((anc_kepmag >= bin - 0.5) & (anc_kepmag < bin + 0.5))[0]
    b_anc[b] = np.nanmedian(anc2[i])
    
  pl.errorbar(bins - 0.025, b_kep, 1.4826 * np.median(np.abs(b_kep - np.median(b_kep))), color = 'y', fmt = 'o', lw = 1, ecolor = 'k', capthick = 1, zorder = 99)
  pl.errorbar(bins + 0.025, b_anc, 1.4826 * np.median(np.abs(b_anc - np.median(b_anc))), color = 'b', fmt = 'o', lw = 1, ecolor = 'k', capthick = 1, zorder = 99)
  
  pl.show(); quit()
  
  fig.savefig('../tex/images/precision1.png' % name, bbox_inches = 'tight')
  fig.savefig('../tex/images/precision1.pdf' % name, bbox_inches = 'tight')
  pl.close()