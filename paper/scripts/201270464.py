#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Figure 9
--------

This script reproduces Figure 9 in the paper.
The source code is available 
`here <https://github.com/rodluger/everest/blob/master/paper/scripts/201270464.py>`_.

.. figure:: ../paper/tex/images/201270464.png
    :width: 600px
    :align: center
    :height: 100px
    :alt: alternate text
    :figclass: align-center

    **Figure 9** EPIC 201270464, a `Kp = 9.4` saturated eclipsing binary. Plotted 
    here is the raw flux (top), the K2SFF flux (center), and the EVEREST flux 
    (bottom). PLD washes out the stellar variability along with most of the 
    eclipses for some saturated stars.

This data was pre-downloaded and pre-detrended, so again this script doesn't
generate the figure *from scratch*. Also, note that the data here are actually 
from EPIC 201270176, but this star
shares the same aperture with 201270464 and 201270464 is much, much brighter,
so it doesn't make a difference, really.
We're not capturing all the flux from 201270464 in this aperture,
so not great for science, but great for the saturation example
in the paper.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
EVEREST_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IMG_PATH = os.path.join(EVEREST_ROOT, 'paper', 'tex', 'images')
sys.path.insert(1, EVEREST_ROOT)
import everest
import argparse
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator

if __name__ == '__main__':

  data = np.load(os.path.join('npz', '201270464.npz'))['data'][()]
  lightcurve = data['lightcurve']
  detrended = data['detrended']
  data = {}
  data.update(lightcurve)
  data.update(detrended)

  time = data['time']
  flux = data['flux'] / np.median(data['flux'])
  fpld = data['fpld'] / np.median(data['fpld'])
  tsff, fsff, _ = np.loadtxt(os.path.join('npz', '201270464.sff'), unpack = True, skiprows = 14)
  fsff /= np.median(fsff)

  fig, ax = pl.subplots(3, figsize = (10,6))
  fig.subplots_adjust(hspace = 0.1, wspace = 0.075)
  ax[0].plot(time, flux, 'k.', alpha = 0.3)
  ax[1].plot(tsff, fsff, 'k.', alpha = 0.3)
  ax[2].plot(time, fpld, 'k.', alpha = 0.3)

  ax[0].set_xticklabels([])
  ax[1].set_xticklabels([])
  for n in range(3):
    ax[n].yaxis.set_major_locator(MaxNLocator(5))
    ax[n].margins(0.01, None)
  ax[0].set_ylabel('Raw Flux', fontsize = 14)
  ax[1].set_ylabel('K2SFF Flux', fontsize = 14)
  ax[2].set_ylabel('EVEREST Flux', fontsize = 14)
  ax[2].set_xlabel('Time (days)', fontsize = 14)
  ax[0].set_ylim(0.75, 1.25)
  ax[1].set_ylim(0.98, 1.01)
  ax[2].set_ylim(0.98, 1.01)
  ax[0].set_yticks([0.75, 0.90, 1.05, 1.20])
  ax[0].set_yticklabels(['0.750', '0.900', '1.050', '1.200'])
  ax[1].set_yticks([0.985, 0.99, 0.995, 1., 1.005])
  ax[2].set_yticks([0.985, 0.99, 0.995, 1., 1.005])

  fig.savefig(os.path.join(IMG_PATH, '201270464.png'), bbox_inches = 'tight')
  fig.savefig(os.path.join(IMG_PATH, '201270464.pdf'), bbox_inches = 'tight')
  pl.close()