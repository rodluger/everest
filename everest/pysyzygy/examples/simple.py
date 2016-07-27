#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`simple.py` - A simple light curve
------------------------------------------

Plots a simple light curve, the planet orbit
as seen from the top and from the observer's
viewpoint, and the orbital elements as a
function of time:

.. figure:: ../pysyzygy/img/simple.png
    :width: 600px
    :align: center
    :height: 100px
    :alt: alternate text
    :figclass: align-center

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import pysyzygy as ps
import numpy as np
import matplotlib.pyplot as pl

if __name__ == '__main__':

  # Transit model kwargs
  kwargs = dict(RpRs = 0.1, b = 0.5, 
                ecc = 0.5, w = 0.,
                MpMs = 0., rhos = 1.4,
                per = 5., t0 = 0.,
                q1 = 0.45, q2 = 0.3,
                maxpts = 20000,
                exptime = ps.KEPLONGEXP)

  # Set up the plot
  fig, ax = pl.subplots(2, 2, figsize = (10, 8))
  fig.subplots_adjust(hspace = 0.35)
  fig.suptitle('Pysyzygy Transit Model', fontsize = 18)
  ax = ax.flatten()

  # The transit lightcurve
  trn = ps.Transit(**kwargs)
  time = np.linspace(-0.15, 0.15, 1000)

  ax[0].plot(time, trn(time, 'binned'), 'b-', label = 'binned')
  ax[0].plot(time, trn(time, 'unbinned'), 'r-', label = 'unbinned')
  ax[0].margins(0., 0.25)
  ax[0].legend(loc = 'lower left', fontsize = 8)
  ax[0].set_xlabel('Time (days)', fontsize = 12)
  ax[0].set_ylabel('Normalized flux', fontsize = 12)
  ax[0].set_title('Transit lightcurve')

  # The orbit
  trn = ps.Transit(fullorbit = True, **kwargs)
  time = np.linspace(-5, 5, 1000)

  ax[1].plot(trn(time, 'x'), trn(time, 'y'), 'k-')
  ax[1].axvline(0., alpha = 0.5, ls = '--', color = 'k')
  ax[1].axhline(0., alpha = 0.5, ls = '--', color = 'k')
  ax[1].set_xlabel('X (stellar radii)', fontsize = 12)
  ax[1].set_ylabel('Y (stellar radii)', fontsize = 12)
  ax[1].set_aspect('equal', 'datalim')
  ax[1].set_title('Observer view')

  ax[2].plot(trn(time, 'x'), trn(time, 'z'), 'k-')
  ax[2].axvline(0., alpha = 0.5, ls = '--', color = 'k')
  ax[2].axhline(0., alpha = 0.5, ls = '--', color = 'k')
  ax[2].set_xlabel('X (stellar radii)', fontsize = 12)
  ax[2].set_ylabel('Z (stellar radii)', fontsize = 12)
  ax[2].set_title('Top view')

  ax[3].plot(time, trn(time, 'M'), label = 'Mean anomaly')
  ax[3].plot(time, trn(time, 'E'), label = 'Eccentric anomaly')
  ax[3].plot(time, trn(time, 'f'), label = 'True anomaly')
  ax[3].margins(0., 0.1)
  ax[3].set_xlabel('Time (days)', fontsize = 12)
  ax[3].set_ylabel('Anomaly (radians)', fontsize = 12)
  ax[3].set_title('Orbital elements')
  ax[3].legend(fontsize = 8)

  pl.show()
