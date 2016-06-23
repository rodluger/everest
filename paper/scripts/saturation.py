#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Figure 8
--------

This script reproduces Figure 8 in the paper.
The source code is available 
`here <https://github.com/rodluger/everest/blob/master/paper/scripts/saturation.py>`_.

.. figure:: ../paper/tex/images/saturation2.png
    :width: 400px
    :align: center
    :height: 100px
    :alt: alternate text
    :figclass: align-center

    **Figure 8** Fractional pixel fluxes :math:`p_{il} / \sum_k{p_{ik}}`
    for quarter 3 of Kepler-3, a `K_p = 9.2` hot-Jupiter host
    observed by the original `Kepler` mission. The panels
    are arranged according to the positions of the pixels on
    the detector, and the data is smoothed and folded on the orbital period
    of Kepler-3b. Saturated pixels are highlighted in red and are labeled
    with an **S**;
    overflow pixels are highlighted in blue and labeled with an **O**. PLD fails for this
    system because the transit signal is present in several of
    the basis vectors.

It also generates the **raw** pixel flux map (i.e., the unnormalized version of 
**Figure 8**), which I think is just as cool looking:

.. figure:: ../paper/tex/images/saturation1.png
    :width: 400px
    :align: center
    :height: 100px
    :alt: alternate text
    :figclass: align-center

    **Figure 8b** Same as above, but showing the unnormalized pixels.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
from everest.config import EVEREST_SRC
IMG_PATH = os.path.join(os.path.dirname(EVEREST_SRC), 'paper', 'tex', 'images')
import everest
from everest.utils import Smooth
import k2plr as kplr
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator

def GetData():
  try:
    data = np.load(os.path.join('npz', 'saturation.npz'))
    x = data['x']
    y = data['y']
    z = data['z']
  except:
    # Get the data
    client = kplr.API()
    kic = client.koi(3.01)
    tpf = kic.get_target_pixel_files(short_cadence = False)
    with tpf[3].open(clobber = False) as f:  
      aperture = f[2].data
      qdata = f[1].data
      quarter = f[0].header['QUARTER']
      crwd = f[1].header['CROWDSAP']
    # Get the arrays
    time = np.array(qdata.field('TIME'), dtype='float64')
    fpix = np.array(qdata.field('FLUX'), dtype='float64')
    fpix_opt = np.array([f[np.where(aperture & 2)] for f in fpix], dtype='float64')
    qual = np.array(qdata.field('QUALITY'), dtype=int)
    # Get bad indices
    t_nan_inds = list(np.where(np.isnan(time))[0]) 
    f_nan_inds = list(np.where(np.isnan(np.sum(fpix_opt, axis = 1)))[0])                           
    qual_inds = []
    for b in [1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17]:
      qual_inds += list(np.where(qual & 2 ** (b - 1))[0])
    bad_inds = np.array(sorted(list(set(qual_inds + t_nan_inds))))
    # Remove them
    time = np.delete(time, bad_inds)
    fpix_opt = np.delete(fpix_opt, bad_inds, 0)
    fpix = np.delete(fpix, bad_inds, 0)
    fsap = np.sum(fpix_opt, axis = 1)
    # Fold and get the pixel fluxes
    per = 4.8878
    t0 = 2.615
    fold = lambda t: (t - t0 - per / 2.) % per - per / 2.
    x = fold(time)
    y = np.zeros((15, 11, len(x)))
    z = np.zeros((15, 11, len(x)))
    for i in range(15):
      for j in range(11):
        f = fpix[:,i,j]
        med = np.nanmedian(f)
        y[i][j] = (f - Smooth(f) + med) / med
        f = fpix[:,i,j] / fsap
        med = np.nanmedian(f)
        z[i][j] = (f - Smooth(f) + med) / med
    np.savez(os.path.join('npz', 'saturation.npz'), x = x, y = y, z = z)
  return x, y, z

if __name__ == '__main__':
  
  # Grab the data
  x, y, z = GetData()

  # Fluxes
  fig, ax = pl.subplots(15,11, sharex = True, sharey = True, figsize = (11,15))
  fig.subplots_adjust(hspace = 0, wspace = 0)
  for i in range(15):
    for j in range(11):
      ax[i,j].plot(x, y[i][j], 'k.', markersize = 0.75)
      ax[i,j].set_xticks([])
      ax[i,j].set_yticks([])
      if np.all(np.isnan(y[i][j])):
        ax[i,j].set_visible(False)
  for j,i in [[5,3],[5,4],[5,5],[5,6],[5,7],[5,8],[5,9],[5,10],[5,11],[6,6],[6,7]]:
    ax[i,j].set_axis_bgcolor((1., 0.9, 0.9))
    ax[i,j].annotate('S', xy = (0.1, 0.9), xycoords = 'axes fraction',
                     ha = 'left', va = 'top', fontweight = 'bold', fontsize = 12)
  for j,i in [[5,2],[5,12],[6,5],[6,8]]:
    ax[i,j].set_axis_bgcolor((0.9, 0.9, 1.))
    ax[i,j].annotate('O', xy = (0.1, 0.9), xycoords = 'axes fraction',
                     ha = 'left', va = 'top', fontweight = 'bold', fontsize = 12)
  ax[0,0].set_xlim(-0.15, 0.15)
  ax[0,0].set_ylim(0.95,1.05)
  fig.savefig(os.path.join(IMG_PATH, 'saturation1.png'), bbox_inches = 'tight')
  pl.close()

  # Fractions
  fig, ax = pl.subplots(15,11, sharex = True, sharey = True, figsize = (11,15))
  fig.subplots_adjust(hspace = 0, wspace = 0)
  for i in range(15):
    for j in range(11):
      ax[i,j].plot(x, z[i][j], 'k.', markersize = 0.75)
      ax[i,j].set_xticks([])
      ax[i,j].set_yticks([])
      if np.all(np.isnan(y[i][j])):
        ax[i,j].set_visible(False)
  for j,i in [[5,3],[5,4],[5,5],[5,6],[5,7],[5,8],[5,9],[5,10],[5,11],[6,6],[6,7]]:
    ax[i,j].set_axis_bgcolor((1., 0.9, 0.9))
    ax[i,j].annotate('S', xy = (0.1, 0.9), xycoords = 'axes fraction',
                     ha = 'left', va = 'top', fontweight = 'bold', fontsize = 12)
  for j,i in [[5,2],[5,12],[6,5],[6,8]]:
    ax[i,j].set_axis_bgcolor((0.9, 0.9, 1.))
    ax[i,j].annotate('O', xy = (0.1, 0.9), xycoords = 'axes fraction',
                     ha = 'left', va = 'top', fontweight = 'bold', fontsize = 12)
  ax[0,0].set_xlim(-0.15, 0.15)
  ax[0,0].set_ylim(0.95,1.05)
  fig.savefig(os.path.join(IMG_PATH, 'saturation2.png'), bbox_inches = 'tight')
  pl.close()