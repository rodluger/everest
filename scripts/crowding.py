#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
crowding.py
-----------

        
'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
EVEREST_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(1, EVEREST_ROOT)
import everest
import numpy as np
import matplotlib.pyplot as pl
from scipy.ndimage import zoom
from scipy.optimize import fmin

def AddApertureContour(ax, nx, ny, aperture):
  '''
  
  '''
  
  # Get the indices
  apidx = np.where(aperture & 1)
  
  # Create the array
  contour = np.zeros((ny,nx))
  contour[apidx] = 1

  # Add padding around the contour mask so that the edges get drawn
  contour = np.lib.pad(contour, 1, everest.utils.PadWithZeros)

  # Zoom in to make the contours look vertical/horizontal
  highres = zoom(contour, 100, order=0, mode='nearest') 
  extent = np.array([0, nx, 0, ny]) + \
           np.array([0, -1, 0, -1]) + \
           np.array([-1, 1, -1, 1])

  # Apply the contours
  ax.contour(highres, levels=[0.5], extent=extent, 
             origin='lower', colors='k', linewidths=3)  

  return

def PlotPostageStamp(EPIC, apnum = 15):
  '''
  
  '''
  
  # Grab the data
  data = everest.GetK2Data(EPIC)
  nearby = data.nearby
  aperture = data.apertures[apnum]
  fpix = data.fpix
  kepmag = data.kepmag
  _, ny, nx = fpix.shape
  total_flux = np.log10(np.nansum(fpix, axis = 0))

  # Plot the data and the aperture
  fig, ax = pl.subplots(1)
  ax.imshow(total_flux, interpolation = 'nearest', alpha = 0.75)
  AddApertureContour(ax, nx, ny, aperture)

  # Crop the image
  ax.set_xlim(-0.7, nx - 0.3)
  ax.set_ylim(-0.7, ny - 0.3)

  # Overplot nearby sources
  neighbors = []  
  def size(k):
    s = 1000 * 2 ** (kepmag - k)
    return s
  for source in nearby:
    ax.scatter(source.x - source.x0, source.y - source.y0, 
               s = size(source.kepmag), 
               c = ['g' if source.epic == EPIC else 'r'],
               alpha = 0.5,
               edgecolor = 'k')
    ax.scatter(source.x - source.x0, source.y - source.y0, 
               marker = r'$%.1f$' % source.kepmag, color = 'w',
               s = 500)
  ax.set_xticklabels([])
  ax.set_yticklabels([])
  
  pl.suptitle('EPIC %d' % EPIC, fontsize = 25, y = 0.98)
  pl.show()
  
  return

PlotPostageStamp(201367065, apnum = 15)
data = everest.GetK2Data(201367065)
total_flux = np.log10(np.nansum(fpix, axis = 0))

def GaussianFit(EPIC):
  data = everest.GetK2Data(EPIC)
  nearby = data.nearby
  aperture = data.apertures[apnum]
  fpix = data.fpix
  kepmag = data.kepmag
  _, ny, nx = fpix.shape

# generate 2 dimensional guassian plot
def Gaus2D(x, y, xc, yc, a, sigma):
  return a*np.exp((- (x - xc) ** 2 - (y - yc) ** 2) / (2 * sigma ** 2))

# define chi squared function
def ChiSq(params):
  xc, yc, a, sigma = params
  result = 0
  for x in range(100):
    for y in range(100):
      result += (Gaus2D(x,y,xc,yc,a,sigma) - data[x,y]) ** 2 / errors[x,y] ** 2
  return result

# Plot sample chi squared curve
fig1 = pl.figure()
sigma_arr = np.linspace(1,10,100)
pl.plot(sigma_arr, [ChiSq([50., 50., 1., s]) for s in sigma_arr])
pl.title('Chi squared', fontsize = 22)

  # Minimize chi squared to find true parameters
guess = [49., 51., 1.2, 7.]
res = fmin(ChiSq, guess)
print(res)

# Generate the model matrix
xc, yc, a, sigma = res
model = np.zeros((100, 100))
for x in range(100):
  for y in range(100):
    model[x,y] = Gaus2D(x, y, xc, yc, a, sigma)
