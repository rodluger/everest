#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Figures 16-17
-------------

This script reproduces Figures 16-17 in the paper; the source code is available
`here <https://github.com/rodluger/everest/blob/master/paper/scripts/detrended_examples.py>`_.

The light curves are pre-downloaded and pre-detrended; I stored them in files in the `npz`
subdirectory.

.. figure:: ../paper/tex/images/detrended0.png
    :width: 600px
    :align: center
    :height: 100px
    :alt: alternate text
    :figclass: align-center

    **Figure 16** De-trended light curves for the campaign 1 star EPIC 201367065 (K2-3, 
    Crossfield et al. 2015). *Top:* The de-trended K2SFF flux (left) and the GP-smoothed 
    flux folded on the periods of the planets b, c, and d (right). *Bottom:* 
    The de-trended EVEREST flux. The 6-hr CDPP is 30.9 ppm for K2SFF and 16.6 ppm 
    for EVEREST, a factor of ∼ 2 improvement.

.. figure:: ../paper/tex/images/detrended1.png
    :width: 600px
    :align: center
    :height: 100px
    :alt: alternate text
    :figclass: align-center

    **Figure 17** De-trended light curves for EPIC 205071984, a campaign 2 star with 
    three known planet candidates (Sinukoff et al. 2015). As in Figure 16, the K2SFF 
    light curve and the folded transits of EPIC 205071984.01 (b), 205071984.02 (c), 
    and 205071984.03 (d) are shown at the top; the equivalent plots for EVEREST are 
    shown at the bottom. The 6-hr CDPP is 56.1 ppm for K2SFF and 24.0 ppm for EVEREST, 
    a factor of ~ 2 improvement.

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

if __name__ == '__main__':
  for i, EPIC, per, t0, ylim, yticks, ylimf, yticksf in \
                          zip([0,1], [201367065, 205071984], 
                              [(10.05448, 24.64745, 44.571), (8.9919156, 20.660987, 31.715868)], 
                              [(2456813.4177, 2456812.2759, 2456826.227), (2456900.927, 2456899.423, 2456903.792)],
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
  
    F1 = lambda t: (t - (t0[0] - 2454833) - per[0] / 2.) % per[0] - per[0] / 2.
    F2 = lambda t: (t - (t0[1] - 2454833) - per[1] / 2.) % per[1] - per[1] / 2.
    F3 = lambda t: (t - (t0[2] - 2454833) - per[2] / 2.) % per[2] - per[2] / 2.
  
    fig = pl.figure(figsize = (20, 6))
    fig.subplots_adjust(wspace = 1., hspace = 0.5)
    ax = [pl.subplot2grid((40, 130), (0, 0), colspan=90, rowspan=20),
          pl.subplot2grid((40, 130), (20, 0), colspan=90, rowspan=20),
          pl.subplot2grid((40, 130), (0, 94), colspan=20, rowspan=20),
          pl.subplot2grid((40, 130), (20, 94), colspan=20, rowspan=20),
          pl.subplot2grid((40, 130), (0, 114), colspan=20, rowspan=10),
          pl.subplot2grid((40, 130), (20, 114), colspan=20, rowspan=10),
          pl.subplot2grid((40, 130), (10, 114), colspan=20, rowspan=10),
          pl.subplot2grid((40, 130), (30, 114), colspan=20, rowspan=10)]

    ax[0].plot(tsff, fsff, 'r.', alpha = 0.3)
    ax[1].plot(time, fpld, 'k.', alpha = 0.3)
    ax[2].plot(F1(time_k2sff), fwhite_sff, 'r.', alpha = 0.75)
    ax[3].plot(F1(time), fwhite, 'k.', alpha = 0.75)
    ax[4].plot(F2(time_k2sff), fwhite_sff, 'r.', alpha = 0.75)
    ax[5].plot(F2(time), fwhite, 'k.', alpha = 0.75)
    ax[6].plot(F3(time_k2sff), fwhite_sff, 'r.', alpha = 0.75)
    ax[7].plot(F3(time), fwhite, 'k.', alpha = 0.75)
  
    ax[2].annotate('b', xy = (0.05, 0.05), ha = 'left', va = 'bottom', xycoords = 'axes fraction', fontweight = 'bold', color = 'r')
    ax[3].annotate('b', xy = (0.05, 0.05), ha = 'left', va = 'bottom', xycoords = 'axes fraction', fontweight = 'bold')
  
    ax[4].annotate('c', xy = (0.05, 0.05), ha = 'left', va = 'bottom', xycoords = 'axes fraction', fontweight = 'bold', color = 'r')
    ax[5].annotate('c', xy = (0.05, 0.05), ha = 'left', va = 'bottom', xycoords = 'axes fraction', fontweight = 'bold')
  
    ax[6].annotate('d', xy = (0.05, 0.05), ha = 'left', va = 'bottom', xycoords = 'axes fraction', fontweight = 'bold', color = 'r')
    ax[7].annotate('d', xy = (0.05, 0.05), ha = 'left', va = 'bottom', xycoords = 'axes fraction', fontweight = 'bold')
  
    for n in [0,1]:
      ax[n].ticklabel_format(useOffset=False)
      ax[n].margins(0.01, None)
      ax[n].set_ylim(ylim)
      ax[n].set_yticks(yticks)
    for n in [2,3,4,5,6,7]:
      ax[n].ticklabel_format(useOffset=False)
      ax[n].set_xlim(-0.35, 0.35)
      if n > 3:
        ax[n].yaxis.tick_right()
      ax[n].set_xticks([-0.25, 0., 0.25])
      for tick in ax[n].get_yticklabels() + ax[n].get_xticklabels():
        tick.set_fontsize(8)
      ax[n].set_yticks(yticksf)
      ax[n].set_ylim(ylimf)
    ax[0].set_xticklabels([])
    ax[2].set_xticklabels([])
    ax[4].set_xticklabels([])
    ax[6].set_xticklabels([])
    ax[0].set_ylabel('K2SFF Flux', fontsize = 18)
    ax[1].set_ylabel('EVEREST Flux', fontsize = 18)
    ax[1].set_xlabel('Time (days)', fontsize = 18)
  
    fig.text(0.795, 0.03, 'Time (days)', ha = 'center', va = 'center', fontsize = 18)
    fig.savefig(os.path.join(IMG_PATH, 'detrended%s.png' % i), bbox_inches = 'tight')
    fig.savefig(os.path.join(IMG_PATH, 'detrended%s.pdf' % i) , bbox_inches = 'tight')