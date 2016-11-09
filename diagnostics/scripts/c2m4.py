#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
c2m4.py
-------

Run Everest 2.0 tests on Campaign 2, Module 4 stars.

'''

import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from everest2.data import GetK2Stars, GetK2Campaign
from everest2.config import EVEREST_SRC, EVEREST_DAT
from everest2.run import DownloadCampaign, RunCampaign
from everest2.model import rPLD
import numpy as np
import matplotlib.pyplot as pl
import os

download = False
run = False
gen = False
plot = True
methods = ['rPLD', 'nPLD', 'bPLD']

# Get stars on the same module as EPIC 205071984
if not os.path.exists(os.path.join(EVEREST_SRC, 'tables', 'C99.csv')):
  epics, kepmags, channels = np.array(GetK2Stars()[2]).T
  epics = np.array(epics, dtype = int)
  inds = np.where(((channels == 9) | (channels == 10) | (channels == 11) | (channels == 12)))
  epics = epics[inds]
  kepmags = kepmags[inds]
  channels = channels[inds]
  with open(os.path.join(EVEREST_SRC, 'tables', 'C99.csv'), 'w') as f:
    for i in range(len(epics)):
      print("%d,%.3f,%d" % (epics[i], kepmags[i], channels[i]), file = f)

# Download
if download:
  DownloadCampaign(99)

# Run
if run:
  RunCampaign(99)

# Generate cdpp file
if gen:
  for method in methods:
    stars = np.array([s[0] for s in GetK2Campaign(99)], dtype = int)
    with open('c2m4_everest2%s.tsv' % method, 'w') as f:
      for i, star in enumerate(stars):
        print('Processing target %d/%d...' % (i + 1, len(stars)))
        nf = os.path.join(EVEREST_DAT, 'c02', 
                         ('%09d' % star)[:4] + '00000', 
                         ('%09d' % star)[4:], method + '.npz')
        try:
          data = np.load(nf)
          print("{:>09d} {:>15.3f} {:>15.3f}".format(star, data['cdpp6'][()], data['cdpp6v'][()]), file = f)
        except:
          print("{:>09d} {:>15.3f} {:>15.3f}".format(star, np.nan, np.nan), file = f)
  
# Plot
if plot:
  for method in methods:
    
    # Everest 1.0
    epic, kp, cdpp6 = np.loadtxt('c2m4_everest.tsv', unpack = True)
    epic = np.array(epic, dtype = int)
    # Everest 2.0
    epic2, cdpp6_2, cdpp6v_2 = np.loadtxt('c2m4_everest2%s.tsv' % method, unpack = True)
    epic2 = np.array(epic2, dtype = int)
    assert np.array_equal(epic, epic2), "OOPS"
    
    '''
    # Figure 1
    # --------
    fig1 = pl.figure()
    pl.scatter(kp, cdpp6, color = 'b')
    pl.scatter(kp, cdpp6_2, color = 'y')
    pl.xlim(8, 17)
    pl.ylim(0, 750)
    pl.title(method)
    '''
    
    # Figure 2
    # --------
    fig2, ax = pl.subplots(1)
    y = (cdpp6_2 - cdpp6) / cdpp6
    ax.plot(kp, y, 'b.', alpha = 0.3)
    ax.set_ylim(-1,1)
    ax.set_xlim(11,17)
    ax.axhline(0, color = 'gray', lw = 2, zorder = -99, alpha = 0.5)
    ax.axhline(0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)
    ax.axhline(-0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)
    # Median
    bins = np.arange(10,16.5,0.5)
    by = np.zeros_like(bins) * np.nan
    for b, bin in enumerate(bins):
      i = np.where((y > -np.inf) & (y < np.inf) & (kp >= bin - 0.5) & (kp < bin + 0.5))[0]
      if len(i) > 10:
        by[b] = np.median(y[i])
    ax.plot(bins, by, 'k-', lw = 2)
    ax.set_ylabel(r'$\frac{\mathrm{CDPP}_{\mathrm{2.0}} - \mathrm{CDPP}_{1.0}}{\mathrm{CDPP}_{1.0}}$', fontsize = 22)
    ax.set_xlabel('Kepler Magnitude', fontsize = 18)  
    pl.title(method)
  
  pl.show()