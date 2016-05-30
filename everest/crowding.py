#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`crowding.py` - Crowding metrics
----------------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import matplotlib.pyplot as pl
import numpy as np
from scipy.optimize import fmin_l_bfgs_b
from scipy.special import erf
from .utils import PadWithZeros
from k2plr.config import KPLR_ROOT
from scipy.ndimage import zoom
import os

def Erf2D(x, y, xc, yc, a, sigma):
  '''

  '''
  
  c = np.sqrt(2) * sigma
  ex = (erf((x + 1 - xc) / c) - erf((x - xc) / c))
  ey = (erf((y + 1 - yc) / c) - erf((y - yc) / c))
  return a * ex * ey

def ChiSq(params, data, err, X, Y):
  '''
  
  '''
  
  y = data - Erf2D(X, Y, *params)
  result = np.nansum(y ** 2 / err ** 2)
  return result
  
def Fit(img, err, xc, yc):
  '''
  
  '''
  
  ny, nx = img.shape
  y, x = np.arange(0,ny), np.arange(0,nx)
  X, Y = np.meshgrid(x, y)
  mxd = np.nanmax(img)
  guess = [xc, yc, mxd, 0.5]
  bounds = [(xc - 2, xc + 2), (yc - 2, yc +  2),
            (0.01 * mxd, 100 * mxd), (0.25, 1.)]
  res = fmin_l_bfgs_b(ChiSq, guess, args = (img, err, X, Y), 
                      approx_grad = True,
                      bounds = bounds)
  params = res[0]
  return params

def Plot(img, params, apidx, C):
  '''
  
  '''
  
  # Compute
  ny, nx = img.shape
  y, x = np.arange(0,ny), np.arange(0,nx)
  X, Y = np.meshgrid(x, y)
  model = Erf2D(X, Y, *params)
  yhr, xhr = np.arange(0,ny,0.01), np.arange(0,nx,0.01)
  Xhr, Yhr = np.meshgrid(xhr, yhr)
  modhr = Erf2D(Xhr, Yhr, *params)
  vmin = 0
  vmax = max(np.nanmax(img), np.nanmax(model))
  
  # Contamination metric
  cont = np.zeros_like(img)
  cont[apidx] = np.abs(img - model)[apidx] / np.nansum(img[apidx])
  
  # Plot
  fig, ax = pl.subplots(2,3, figsize = (12,8))
  fig.subplots_adjust(top = 0.95, bottom = 0.05, left = 0.07, right = 0.9)
  ax[0,0].imshow(img, interpolation = 'nearest', vmin = vmin, vmax = vmax)
  ax[0,1].imshow(modhr, interpolation = 'nearest', extent = (-1,nx,ny,-1))
  ax[1,0].imshow(model, interpolation = 'nearest', vmin = vmin, vmax = vmax)
  ax[1,1].imshow(img - model, interpolation = 'nearest', vmin = vmin, vmax = vmax)
  im = ax[0,2].imshow(cont, interpolation = 'nearest')
  cax = fig.add_axes([0.915, 0.55, 0.015, 0.39])
  fig.colorbar(im, cax=cax, orientation='vertical')
  ax[1,2].plot(range(len(C)), C, 'b.')
  ax[1,2].axhline(np.nanmedian(C), color = 'k', ls = '--')
  ax[1,2].yaxis.tick_right()
  ax[1,2].set_xticks([0,len(C)-1])
  ax[1,2].set_xticklabels(['Start', 'End'])
  ax[1,2].set_aspect(1./ax[1,2].get_data_ratio() * ny / nx)
  
  # Labels
  ax[0,0].set_title('Data', fontsize = 20)
  ax[0,1].set_title('Model', fontsize = 20)
  ax[1,0].set_title('Model (binned)', fontsize = 20)
  ax[1,1].set_title('Data - model', fontsize = 20)
  ax[0,2].set_title('Contamination', fontsize = 20)
  ax[1,2].set_title('Contamination', fontsize = 20)

  # Apply the aperture contours
  contour = np.zeros((ny,nx))
  contour[apidx] = 1
  contour = np.lib.pad(contour, 1, PadWithZeros)
  highres = zoom(contour, 100, order=0, mode='nearest') 
  extent = np.array([0, nx, 0, ny]) + \
           np.array([0, -1, 0, -1]) + \
           np.array([-1, 1, -1, 1])
  for axis in ax.flatten()[:-1]:
    axis.contour(highres, levels=[0.5], extent=extent, 
                 origin='lower', colors='k', linewidths=3)
    axis.set_xlim(-0.7, nx - 0.3)
    axis.set_ylim(-0.7, ny - 0.3)

  return fig, ax

def Contamination(EPIC, fpix, perr, apidx, nearby, plot = False):
  '''
  
  '''
  
  # Get the source position according to MAST
  source = nearby[np.where([s.epic == EPIC for s in nearby])[0]]
  xc, yc = source.x - source.x0, source.y - source.y0  

  # Get the contamination evolution
  npts = 50
  C = np.zeros(npts)
  ny, nx = fpix[0].shape
  cont = np.zeros_like(fpix[0])
  y, x = np.arange(0,ny), np.arange(0,nx)
  X, Y = np.meshgrid(x, y)
  for i, n in enumerate(np.linspace(0, len(fpix) - 1, npts)):
    params = Fit(fpix[int(n)], perr[int(n)], xc, yc)
    model = Erf2D(X, Y, *params)
    cont[apidx] = np.abs(fpix[int(n)] - model)[apidx] / np.nansum(fpix[int(n)][apidx])
    C[i] = np.nanmax(cont)

  # Plot the first timestamp fit
  if plot:
    params = Fit(fpix[0], perr[0], xc, yc)
    fig, ax = Plot(fpix[0], params, apidx, C)
    fig.savefig(os.path.join(KPLR_ROOT, 'data', 'everest', str(EPIC), 'contamination.png'))
    pl.close()

  return np.nanmedian(C)