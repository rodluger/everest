#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Figure 7
--------

This script reproduces Figure 7 in the paper, which shows the eclipsing binary
EPIC 202072563 de-trended with both `everest` and `K2SFF`.
The source code is available 
`here <https://github.com/rodluger/everest/blob/master/paper/scripts/202072563.py>`_.

I pre-ran this target with :py:mod:`everest` and stored the raw and de-trended data
in the file 
`202072563.npz <https://github.com/rodluger/everest/blob/master/paper/scripts/202072563.npz>`_,
so this script isn't as transparent as the others. But it shouldn't be hard to write
your own script to plot this up from scratch.

.. figure:: ../paper/tex/images/202072563.png
    :width: 500px
    :align: center
    :height: 100px
    :alt: alternate text
    :figclass: align-center

    **Figure 7** Primary and secondary eclipses of EPIC 202072563, folded on a period 
    of 2.1237 days. From left to right, the columns show the raw SAP flux, the EVEREST 
    flux, the EVEREST flux obtained when explicitly masking the eclipses, and the K2SFF 
    flux. Masking is critical to preserving the transit depth during de-trending for 
    any pipeline that does not explicitly include a transit or eclipse model.
     
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
  data = np.load(os.path.join('npz', '202072563.npz'))
  time = data['time']
  flux = data['flux'] / np.median(data['flux']) + 0.005
  # Note that we adjust the baseline slightly so that ingress and egress are at 1.
  # The out-of-transit variability of the star throws off the median!
  fwhite = data['fwhite']
  fwhite_sff = data['fwhite_sff']
  fwhite_nomask = data['fwhite_nomask']

  # Eclipses
  per = 2.1237
  p0 = 56775.488
  s0 = 56775.488 + per / 2.
  P = lambda t: (t - (p0 - 54833) - per / 2.) % per - per / 2.
  S = lambda t: (t - (s0 - 54833) - per / 2.) % per - per / 2.

  fig, ax = pl.subplots(2,4, sharex = True, sharey = True, figsize = (12,6))
  fig.subplots_adjust(hspace = 0.1, wspace = 0.075)
  ax[0,0].plot(P(time), flux, 'k.', alpha = 0.3)
  ax[0,3].plot(P(time), fwhite_sff, 'k.', alpha = 0.3)
  ax[0,1].plot(P(time), fwhite_nomask, 'k.', alpha = 0.3)
  ax[0,2].plot(P(time), fwhite, 'k.', alpha = 0.3)

  ax[1,0].plot(S(time), flux, 'k.', alpha = 0.3)
  ax[1,3].plot(S(time), fwhite_sff, 'k.', alpha = 0.3)
  ax[1,1].plot(S(time), fwhite_nomask, 'k.', alpha = 0.3)
  ax[1,2].plot(S(time), fwhite, 'k.', alpha = 0.3)

  ax[0,0].set_xlim(-0.25, 0.25)
  ax[0,0].set_ylim(0.75, 1.05)
  ax[0,0].yaxis.set_ticks([0.75, 0.8, 0.85, 0.9, 0.95, 1.0])
  ax[0,0].set_title('Raw SAP', fontsize = 20)
  ax[0,3].set_title('K2SFF', fontsize = 20)
  ax[0,1].set_title('EVEREST', fontsize = 20)
  ax[0,2].set_title('EVEREST (Masked)', fontsize = 20)

  ax[0,0].set_ylabel('Primary', fontsize = 18)
  ax[1,0].set_ylabel('Secondary', fontsize = 18)
  ax[1,0].set_xlabel('Time (days)', fontsize = 18)
  ax[1,1].set_xlabel('Time (days)', fontsize = 18)
  ax[1,2].set_xlabel('Time (days)', fontsize = 18)
  ax[1,3].set_xlabel('Time (days)', fontsize = 18)

  fig.savefig(os.path.join(IMG_PATH, '202072563.png'), bbox_inches = 'tight')
  fig.savefig(os.path.join(IMG_PATH, '202072563.pdf'), bbox_inches = 'tight')
  pl.close()