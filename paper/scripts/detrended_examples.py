#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
detrended_examples.py
---------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
import os, sys
EVEREST_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IMG_PATH = os.path.join(EVEREST_ROOT, 'paper', 'tex', 'images')
sys.path.insert(1, EVEREST_ROOT)
import everest
import argparse
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator

for i, EPIC, per, t0, ylim, yticks, ylimf, yticksf in \
                        zip([0,1], [201367065, 205071984], 
                            [10.05448, 8.9919156], 
                            [2456813.4177, 2456900.927],
                            [[0.997, 1.004], [0.988, 1.012]],
                            [[0.998, 0.999, 1., 1.001, 1.002, 1.003], 
                             [0.99, 0.995, 1., 1.005, 1.01]],
                            [[0.998, 1.0005], [0.9955, 1.0006]],
                            [[0.998, 0.999, 1.], [0.996, 0.997, 0.998, 0.999, 1.]]):

  data = np.load(os.path.join('npz', '%s.npz' % EPIC))['data'][()]
  lightcurve = data['lightcurve']
  detrended = data['detrended']
  data = {}
  data.update(lightcurve)
  data.update(detrended)

  time = data['time']
  time_k2sff = data['time_sff']
  flux = data['flux'] / np.median(data['flux'])
  fpld = data['fpld'] / np.median(data['fpld'])
  fwhite = data['fwhite'] / np.median(data['fwhite'])
  fwhite_sff = data['fwhite_sff'] / np.median(data['fwhite_sff'])
  tsff, fsff, _ = np.loadtxt(os.path.join('npz', '%s.sff' % EPIC), unpack = True, skiprows = 14)
  fsff /= np.median(fsff)
  F = lambda t: (t - (t0 - 2454833) - per / 2.) % per - per / 2.
  
  fig = pl.figure(figsize = (16, 6))
  fig.subplots_adjust(wspace = 0.1, hspace = 0.05)
  ax = [pl.subplot2grid((2, 10), (0, 0), colspan=7, rowspan=1),
        pl.subplot2grid((2, 10), (1, 0), colspan=7, rowspan=1),
        pl.subplot2grid((2, 10), (0, 7), colspan=3, rowspan=1),
        pl.subplot2grid((2, 10), (1, 7), colspan=3, rowspan=1)]

  ax[0].plot(tsff, fsff, 'r.', alpha = 0.3)
  ax[1].plot(time, fpld, 'k.', alpha = 0.3)
  ax[2].plot(F(time_k2sff), fwhite_sff, 'r.', alpha = 0.75)
  ax[3].plot(F(time), fwhite, 'k.', alpha = 0.75)
  
  for n in [0,1]:
    ax[n].ticklabel_format(useOffset=False)
    ax[n].margins(0.01, None)
    ax[n].set_ylim(ylim)
    ax[n].set_yticks(yticks)
  for n in [2,3]:
    ax[n].ticklabel_format(useOffset=False)
    ax[n].set_xlim(-0.35, 0.35)
    ax[n].yaxis.tick_right()
    ax[n].set_xticks([-0.25, 0., 0.25])
    ax[n].set_yticks(yticksf)
    ax[n].set_ylim(ylimf)
  ax[0].set_xticklabels([])
  ax[2].set_xticklabels([])
  
  ax[0].set_ylabel('K2SFF Flux', fontsize = 18)
  ax[1].set_ylabel('EVEREST Flux', fontsize = 18)
  ax[1].set_xlabel('Time (days)', fontsize = 18)
  ax[3].set_xlabel('Time (days)', fontsize = 18)

  fig.savefig(os.path.join(IMG_PATH, 'detrended%s.png' % i), bbox_inches = 'tight')
  fig.savefig(os.path.join(IMG_PATH, 'detrended%s.pdf' % i) , bbox_inches = 'tight')