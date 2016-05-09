#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
injections.py
-------------    
        
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
nums = [1, 2, 3, 4, 5, 6]
ranges = [(0.5, 1.5), (0.5, 1.5), 
          (0.5, 1.5), (0.5, 1.5), 
          (0., 2.), (0., 2.)]
depths = [1e-2, 1e-2, 1e-3, 1e-3, 1e-4, 1e-4]
nbins = [30, 30, 30, 30, 30, 30]
ymax = [0.5, 0.5, 0.3, 0.3, 0.1, 0.1]
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
    DC = data['DC']
    EC = data['EC']
  else:
    D = []; DC = []
    E = []; EC = []
    for i, folder, depth, rng, n in zip(range(len(folders)), folders, depths, ranges, nums):
      nstars = 0
      nstars_ctrl = 0
      everest_depth = []
      everest_err = []
      control_depth = []
      control_err = []
      
      # Loop over all campaigns
      for campaign in range(7):
      
        # Get the targets
        stars = [s for s in os.listdir(os.path.join(EVEREST_ROOT, 'output', 'C%02d' % campaign)) if 
                 os.path.exists(os.path.join(EVEREST_ROOT, 'output', 
                 'C%02d' % campaign, s, folder, '%s.ctrl1.inj' % s))]
        nstars += len(stars)

        # Now get the depths in the injection and the control
        for star in stars:
        
          # Injection
          with open(os.path.join(EVEREST_ROOT, 'output', 'C%02d' % campaign, star, 
                    folder, '%s.inj' % star), 'r') as f:
            lines = f.readlines()
            a, _, b = re.findall('([-0-9.]+)', lines[4])
            everest_depth.append(float(a))
            everest_err.append(float(b))
      
      # Loop over all campaigns
      for campaign in range(7):
      
        # Get the targets
        stars = [s for s in os.listdir(os.path.join(EVEREST_ROOT, 'output', 'C%02d' % campaign)) if 
                 os.path.exists(os.path.join(EVEREST_ROOT, 'output', 
                 'C%02d' % campaign, s, 'default', '%s.ctrl1.inj' % s))]
        nstars_ctrl += len(stars)
        
        # Now get the depths in the injection and the control
        # These must be first generated with ``inject_control.py``
        for star in stars:
        
          # Control injection
          if os.path.exists((os.path.join(EVEREST_ROOT, 'output', 'C%02d' % campaign, star, 
                             'default', '%s.ctrl%d.inj' % (star, n)))):
            with open(os.path.join(EVEREST_ROOT, 'output', 'C%02d' % campaign, star, 
                      'default', '%s.ctrl%d.inj' % (star, n)), 'r') as f:
              lines = f.readlines()
              a, _, b = re.findall('([-0-9.]+)', lines[4])
              control_depth.append(float(a))
              control_err.append(float(b))
      
      # Append to running list
      D.append(everest_depth)
      E.append(everest_err)
      DC.append(control_depth)
      EC.append(control_err)
        
      print('%s: %d/%d' % (folder, nstars, nstars_ctrl))
      
    np.savez(os.path.join('npz', 'injections.npz'), D = D, E = E, DC = DC, EC = EC)
  
  return D, E, DC, EC

# Set up the figure
fig, ax = pl.subplots(3,2, figsize = (10,12))
fig.subplots_adjust(hspace = 0.25)
ax = ax.flatten()
ax[0].set_title(r'Masked', fontsize = 18)
ax[1].set_title(r'Unmasked', fontsize = 18)
ax[0].set_ylabel(r'D$_0$ = 10$^{-2}$', rotation = 90, fontsize = 18, labelpad = 10)
ax[2].set_ylabel(r'D$_0$ = 10$^{-3}$', rotation = 90, fontsize = 18, labelpad = 10)
ax[4].set_ylabel(r'D$_0$ = 10$^{-4}$', rotation = 90, fontsize = 18, labelpad = 10)

# Get the depths
D, E, DC, EC = GetDepths(clobber = True)

# Plot
for i, axis in enumerate(ax): 

  # Plot the histograms
  everest_depth = np.array(D[i]) / depths[i]
  control_depth = np.array(DC[i]) / depths[i]
  
  try:
    axis.hist(everest_depth, bins = nbins[i], range = ranges[i], color = 'b', 
              histtype = 'step', weights = np.ones_like(everest_depth)/len(everest_depth))
  except:
    pass

  try:
    axis.hist(control_depth, bins = nbins[i], range = ranges[i], color = 'r', 
              histtype = 'step', weights = np.ones_like(control_depth)/len(control_depth))
  except:
    pass
    
  axis.axvline(1., color = 'k', ls = '--')
  
  # Indicate the fraction above and below
  if len(everest_depth):
    au = len(np.where(everest_depth > ranges[i][1])[0]) / len(everest_depth)
    al = len(np.where(everest_depth < ranges[i][0])[0]) / len(everest_depth)
    axis.annotate('%.2f' % al, xy = (0.01, 0.95), xycoords = 'axes fraction', 
                  xytext = (0.1, 0.95), ha = 'left', va = 'center', color = 'b',
                  arrowprops = dict(arrowstyle="->",color='b'))
    axis.annotate('%.2f' % au, xy = (0.99, 0.95), xycoords = 'axes fraction', 
                  xytext = (0.9, 0.95), ha = 'right', va = 'center', color = 'b',
                  arrowprops = dict(arrowstyle="->",color='b'))
  
  if len(control_depth):  
    cu = len(np.where(control_depth > ranges[i][1])[0]) / len(control_depth)
    cl = len(np.where(control_depth < ranges[i][0])[0]) / len(control_depth)
    axis.annotate('%.2f' % cl, xy = (0.01, 0.88), xycoords = 'axes fraction', 
                  xytext = (0.1, 0.88), ha = 'left', va = 'center', color = 'r',
                  arrowprops = dict(arrowstyle="->",color='r'))
    axis.annotate('%.2f' % cu, xy = (0.99, 0.88), xycoords = 'axes fraction', 
                  xytext = (0.9, 0.88), ha = 'right', va = 'center', color = 'r',
                  arrowprops = dict(arrowstyle="->",color='r'))
                
  # Indicate the median
  axis.annotate('M = %.2f' % np.median(everest_depth), xy = (0.3, 0.5), ha = 'right',
                xycoords = 'axes fraction', color = 'b', fontsize = 14)
  axis.annotate('M = %.2f' % np.median(control_depth), xy = (0.7, 0.5), ha = 'left',
                xycoords = 'axes fraction', color = 'r', fontsize = 14)
  
  # Tweaks
  axis.set_xticks(xticks[i])
  axis.set_ylim(0, ymax[i])
  axis.set_xlabel(r'D/D$_0$', fontsize = 14)

pl.show()

#fig.savefig('../tex/images/injections.png', bbox_inches = 'tight')
#fig.savefig('../tex/images/injections.pdf', bbox_inches = 'tight')