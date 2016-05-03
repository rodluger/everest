#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
acor.py
-------
Illustrates the GP optimization procedure for EPIC 201497682.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
EVEREST_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IMG_PATH = os.path.join(EVEREST_ROOT, 'paper', 'tex', 'images')
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator

def PlotAcor():
  '''
  
  '''
  
  # This .npz file was generated a while back in an older
  # version of the EVEREST code...
  data = np.load(os.path.join('npz', 'acor.npz'))
  time = data['time']
  fpld = data['fpld']
  flux = data['flux']
  tfull = data['tfull']
  pfull = data['pfull']
  period = data['period']
  PS = data['PS']
  pers = data['pers']
  lags = data['lags']
  acor = data['acor']
  sigma = data['sigma']
  # Hack: Get the ``^2``s out of the ``Exp`` kernel for clarity
  #strkernel = data['strkernel']
  strkernel = '$15.6^2$ $\\times$ $\\mathbf{Exp}$($55.1$) $\\times$ ' + \
              '$\\mathbf{Cos}$($29.6$) $+$ $20.4^2$ $\\times$ ' + \
              '$\\mathbf{Exp}$($31.3$) $\\times$ $\\mathbf{Cos}$($11.8$)'
  kernel = data['kernel']
  EPIC = data['EPIC']
  
  # Plot
  fig, ax = pl.subplots(2, 2, figsize = (20, 8), squeeze = True)
  ax = ax.T.flatten()
  fig.subplots_adjust(left = 0.1, right = 0.95, hspace = 0.3, wspace = 0.15)
  
  # The raw data
  scale = 1.e4
  scale_str = '$10^4$'
  ax[0].plot(time, flux / scale, 'k.', alpha = 0.3, label = 'SAP flux')
  ax[0].plot(time, (flux - fpld) / scale, 'r-', alpha = 0.5, label = 'PLD model')
  ax[0].margins(0.01, None)
  ax[0].set_ylim(43000 / scale, 43600 / scale)
  
  # The detrended data
  ax[1].plot(time, fpld, 'k.', alpha = 0.3)
  ax[1].margins(0.01, None)
  
  # The periodogram
  ax[2].plot(period, PS, '-', c='black', lw=1, zorder=1)
  if len(pers):
    ax[2].axvline(pers[0], color = 'r', alpha = 0.5, label = 'Peak periods')
    [ax[2].axvline(p, color = 'r', alpha = 0.5) for p in pers[1:]]
  ax[2].legend(loc = 'upper left')
  ax[2].margins(0, None)
  
  # The autocorrelation function and fit
  ax[3].plot(lags, acor, 'k-', lw = 2, alpha = 1.)
  ax[3].fill_between(lags, acor - sigma, acor + sigma, alpha = 0.1, color = 'k')
  ax[3].plot(lags, kernel, 'r-', alpha = 1., label = strkernel)
  ax[3].legend(loc = 'upper right', fontsize = 14)
  ax[3].margins(0, None)

  # Labels
  ax[0].legend(loc = 'upper left', numpoints = 3)
  ax[0].set_ylabel('Flux / %s' % scale_str, fontsize = 18)
  ax[0].set_xlabel('Time (days)', fontsize = 18)
  ax[1].set_ylabel('Residuals', fontsize = 18, labelpad = 1.)
  ax[1].set_xlabel('Time (days)', fontsize = 18)
  ax[2].set_ylabel('Power', fontsize = 18, labelpad = 15.)
  ax[2].set_xlabel('Period (days)', fontsize = 18)
  ax[3].set_ylabel('Autocorrelation', fontsize = 18)
  ax[3].set_xlabel('Time Lag (days)', fontsize = 18)
  
  # Save
  fig.savefig(os.path.join(IMG_PATH, 'acor.png'), bbox_inches = 'tight')
  fig.savefig(os.path.join(IMG_PATH, 'acor.pdf'), bbox_inches = 'tight')

PlotAcor()