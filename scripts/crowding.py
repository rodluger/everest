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
sys.path.append('/Users/nks1994/Documents/Research/PyKE')
#reload(sys)
#sys.setdefaultencoding('utf-8')
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as pl
import everest
import scipy
from scipy.ndimage import zoom
from scipy.optimize import fmin
from scipy.special import erf
import k2plr
from k2plr.config import KPLR_ROOT
import glob
import pyfits
import kepio
import kepfunc


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
  # pl.show()
  
  return

PlotPostageStamp(201367065, apnum = 15)

data = everest.GetK2Data(201367065)
fpix = data.fpix
_, ny, nx = fpix.shape
kepmag = data.kepmag
total_flux = np.transpose(fpix[0])
errors = np.transpose(data.perr[0])

# generate 2D error function
def Erf2D(x, y, xc, yc, a, b, sigma):
  '''
  A two-dimensional error function (the integral of a Gaussian over a pixel).
  
  '''
  
  c = np.sqrt(2) * sigma
  ex = (erf((x + 1 - xc) / c) - erf((x - xc) / c))
  ey = (erf((y + 1 - yc) / c) - erf((y - yc) / c))
  return a * ex * ey + b

# define chi squared function
def ChiSq(params):
  xc, yc, a, b, sigma = params
  result = 0
  for x in range(nx):
    for y in range(ny):
      if np.isnan(total_flux[x,y]):
        continue
      else:
        result += (Erf2D(x,y,xc,yc,a,b,sigma) - total_flux[x,y]) ** 2 / errors[x,y] ** 2
  return result

# plot the data
fig1 = pl.figure()
pl.imshow(total_flux)
pl.title('Data', fontsize = 22)
pl.colorbar()
'''
# Plot sample chi squared curve
fig2 = pl.figure()
sigma_arr = np.linspace(1,10,100)
pl.plot(sigma_arr, [ChiSq([9.3, 7.0, 100000., 250., s]) for s in sigma_arr])
pl.xlim(0, 10)
pl.title('Chi squared', fontsize = 22)

  # Minimize chi squared to find true parameters
guess = [9.3, 7.0, 100000., 250., 1.3]
res = fmin(ChiSq, guess)
print(res)

# Generate the model matrix
xc, yc, a, b, sigma = res
model = np.zeros((nx, ny))
for x in range(nx):
  for y in range(ny):
    model[x,y] = Erf2D(x, y, xc, yc, a, b, sigma)

# plot the model
fig3 = pl.figure()
pl.imshow(model)
pl.colorbar()
pl.title('Model', fontsize = 22)

# Plot the difference
fig4 = pl.figure()
pl.imshow(total_flux - model)
pl.colorbar()
pl.title('Data - Model', fontsize = 22)
'''


'''
** KEPLER PRF **
'''

epic = 201367065

prfdir = '/Users/nks1994/Documents/Research/KeplerPRF'
logfile = 'test.log'
verbose = True

client = k2plr.API()
star = client.k2_star(epic)
tpf = star.get_target_pixel_files(fetch = True)[0]
ftpf = os.path.join(KPLR_ROOT, 'data', 'k2', 'target_pixel_files', '%d' % epic, tpf._filename)
tdim5 = pyfits.getheader(ftpf,1)['TDIM5']
xdim = int(tdim5.strip().strip('(').strip(')').split(',')[0])
ydim = int(tdim5.strip().strip('(').strip(')').split(',')[1])
module = pyfits.getheader(ftpf,0)['MODULE']
output = pyfits.getheader(ftpf,0)['OUTPUT']
column = pyfits.getheader(ftpf,2)['CRVAL1P']
row = pyfits.getheader(ftpf,2)['CRVAL2P']

if int(module) < 10:
  prefix = 'kplr0'
else:
  prefix = 'kplr'
prfglob = prfdir + '/' + prefix + str(module) + '.' + str(output) + '*' + '_prf.fits'
prffile = glob.glob(prfglob)[0]

prfn = [0,0,0,0,0]
crpix1p = np.zeros((5),dtype='float32')
crpix2p = np.zeros((5),dtype='float32')
crval1p = np.zeros((5),dtype='float32')
crval2p = np.zeros((5),dtype='float32')
cdelt1p = np.zeros((5),dtype='float32')
cdelt2p = np.zeros((5),dtype='float32')
for i in range(5):
    prfn[i], crpix1p[i], crpix2p[i], crval1p[i], crval2p[i], cdelt1p[i], cdelt2p[i], status \
        = kepio.readPRFimage(prffile,i+1,logfile,verbose) 
prfn = np.array(prfn)
PRFx = np.arange(0.5,np.shape(prfn[0])[1]+0.5)
PRFy = np.arange(0.5,np.shape(prfn[0])[0]+0.5)
PRFx = (PRFx - np.size(PRFx) / 2) * cdelt1p[0]
PRFy = (PRFy - np.size(PRFy) / 2) * cdelt2p[0]

prf = np.zeros(np.shape(prfn[0]), dtype='float32')
prfWeight = np.zeros((5), dtype='float32')
for i in range(5):
    prfWeight[i] = np.sqrt((column - crval1p[i])**2 + (row - crval2p[i])**2)
    if prfWeight[i] == 0.0:
        prfWeight[i] = 1.0e-6
    prf = prf + prfn[i] / prfWeight[i]
prf = prf / np.nansum(prf) / cdelt1p[0] / cdelt2p[0]

prfDimY = int(ydim / cdelt1p[0])
prfDimX = int(xdim / cdelt2p[0])
PRFy0 = (np.shape(prf)[0] - prfDimY) / 2
PRFx0 = (np.shape(prf)[1] - prfDimX) / 2

splineInterpolation = scipy.interpolate.RectBivariateSpline(PRFx,PRFy,prf)

pl.imshow(prf)
pl.show()

args = (DATx,DATy,DATimg,ERRimg,nsrc,splineInterpolation,float(x[0]),float(y[0]))
ans = fmin_powell(kepfunc.PRF,guess,args=args,xtol=xtol,
                  ftol=ftol,disp=False)