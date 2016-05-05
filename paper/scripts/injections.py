#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
injections.py
-------------
       
TODO: Work on this!     
        
'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
EVEREST_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as pl
import re

# Some hard-coded stuff
folders = ['inject_0.0100m', 'inject_0.0100u', 'inject_0.0010m', 
           'inject_0.0010u', 'inject_0.0001m', 'inject_0.0001u']
ranges = [(0.5, 1.5), (0.5, 1.5), 
          (0.5, 1.5), (0.5, 1.5), 
          (0., 2.), (0., 2.)]
depths = [1e-2, 1e-2, 1e-3, 1e-3, 1e-4, 1e-4]
nbins = [30, 30, 30, 30, 30, 30]
ymax = [0.6, 0.6, 0.25, 0.25, 0.2, 0.2]
xticks = [[0.5, 0.75, 1., 1.25, 1.5],
          [0.5, 0.75, 1., 1.25, 1.5],
          [0.5, 0.75, 1., 1.25, 1.5],
          [0.5, 0.75, 1., 1.25, 1.5],
          [0., 0.5, 1., 1.5, 2.0],
          [0., 0.5, 1., 1.5, 2.0]]

def GetDepths(clobber = False):
  '''
  
  '''
    
  if not clobber and os.path.exists(os.path.join('npz', 'injections.npz')):
    data = np.load(os.path.join('npz', 'injections.npz'))
    D = data['D']
    E = data['E']
  else:
    D = []
    E = []
    for i, folder, depth, rng in zip(range(len(folders)), folders, depths, ranges):
      nstars = 0
      for campaign in range(7):
        stars = [s for s in os.listdir(os.path.join(EVEREST_ROOT, 'output', 'C%02d' % campaign)) if 
                 os.path.exists(os.path.join(EVEREST_ROOT, 'output', 
                 'C%02d' % campaign, s, folder, '%s.inj' % s))]
        nstars += len(stars)
        everest_depth = []
        everest_err = []
        for star in stars:
          with open(os.path.join(EVEREST_ROOT, 'output', 'C%02d' % campaign, star, 
                    folder, '%s.inj' % star), 'r') as f:
            lines = f.readlines()
            a, _, b = re.findall('([-0-9.]+)', lines[4])
            everest_depth.append(float(a))
            everest_err.append(float(b))
      print('%s: %d' % (folder, nstars))
      D.append(everest_depth)
      E.append(everest_err)
    np.savez(os.path.join('npz', 'injections.npz'), D = D, E = E)
  return D, E

fig, ax = pl.subplots(3,2, figsize = (10,12))
fig.subplots_adjust(hspace = 0.25)
ax = ax.flatten()
ax[0].set_title(r'Masked', fontsize = 18)
ax[1].set_title(r'Unmasked', fontsize = 18)
ax[0].set_ylabel(r'D$_0$ = 10$^{-2}$', rotation = 90, fontsize = 18, labelpad = 10)
ax[2].set_ylabel(r'D$_0$ = 10$^{-3}$', rotation = 90, fontsize = 18, labelpad = 10)
ax[4].set_ylabel(r'D$_0$ = 10$^{-4}$', rotation = 90, fontsize = 18, labelpad = 10)

D, E = GetDepths()

for i, axis in enumerate(ax): 

  # Plot the histograms
  everest_depth = np.array(D[i]) / depths[i]

  axis.hist(everest_depth, bins = nbins[i], range = ranges[i], color = 'b', 
            histtype = 'step', weights = np.ones_like(everest_depth)/len(everest_depth))
  axis.axvline(1., color = 'k', ls = '--')
  
  # Indicate the fraction above and below
  au = len(np.where(everest_depth < ranges[i][0])[0]) / len(everest_depth)
  al = len(np.where(everest_depth > ranges[i][1])[0]) / len(everest_depth)
  axis.annotate('%.2f' % al, xy = (0.01, 0.95), xycoords = 'axes fraction', 
                xytext = (0.1, 0.95), ha = 'left', va = 'center', color = 'b',
                arrowprops = dict(arrowstyle="->",color='b'))
  axis.annotate('%.2f' % au, xy = (0.99, 0.95), xycoords = 'axes fraction', 
                xytext = (0.9, 0.95), ha = 'right', va = 'center', color = 'b',
                arrowprops = dict(arrowstyle="->",color='b'))  
  
  # Indicate the median
  axis.annotate('M = %.2f' % np.median(everest_depth), xy = (0.3, 0.5), ha = 'right',
                xycoords = 'axes fraction', color = 'b', fontsize = 14)
  
  # Tweaks
  axis.set_xticks(xticks[i])
  axis.set_ylim(0, ymax[i])
  axis.set_xlabel(r'D/D$_0$', fontsize = 14)

pl.show()

#fig.savefig('../tex/images/injections.png', bbox_inches = 'tight')
#fig.savefig('../tex/images/injections.pdf', bbox_inches = 'tight')