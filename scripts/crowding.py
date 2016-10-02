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
from scipy.optimize import fmin_powell
from scipy.special import erf
import k2plr
from k2plr.config import KPLR_ROOT
import glob
import pyfits
import kepio
import kepfunc
from numpy import empty

# define EPIC ID
epic = 205998445

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

PlotPostageStamp(epic, apnum = 15)

data = everest.GetK2Data(epic)
fpix = data.fpix
_, ny, nx = fpix.shape
kepmag = data.kepmag
total_flux = fpix[0]
mean_flux = np.sum(fpix[n] for n in range(len(fpix))) / len(fpix)
errors = data.perr[0]
mean_errs = np.sum(data.perr[n] for n in range(len(data.perr))) / len(data.perr)

'''
** KEPLER PRF **
'''

# set PRF directory
prfdir = '/Users/nks1994/Documents/Research/KeplerPRF'
logfile = 'test.log'
verbose = True

# read PRF header data
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
with pyfits.open(ftpf) as f:
    xpos = f[1].data['pos_corr1']
    ypos = f[1].data['pos_corr2']

if int(module) < 10:
  prefix = 'kplr0'
else:
  prefix = 'kplr'
prfglob = prfdir + '/' + prefix + str(module) + '.' + str(output) + '*' + '_prf.fits'
prffile = glob.glob(prfglob)[0]

# create PRF matrix
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

# generate PRF
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

# interpolate function over the PRF
splineInterpolation = scipy.interpolate.RectBivariateSpline(PRFx,PRFy,prf)

DATx = np.arange(column,column+xdim)
DATy = np.arange(row,row+ydim)

nsrc = 1
neighbors = []

# fit PRF model to pixel data
guess = [218000,28.76,126.67,100,28,121]
args = (DATx,DATy,total_flux,errors,nsrc,splineInterpolation,np.mean(DATx),np.mean(DATy))
ans = fmin_powell(kepfunc.PRF,guess,args=args,xtol=1.0e-4,ftol=1.0e-4,disp=True)

print("ans:" + str(ans))
print("Number of sources = " + str(nsrc))

prffit_guess = kepfunc.PRF2DET([ans[0]], [DATx[0] + ans[1]], [DATy[0] + ans[2]],DATx,DATy,1.0,1.0,0,splineInterpolation)
f=empty((nsrc));x=empty((nsrc));y=empty((nsrc))

for i in range(nsrc):
    src = (nsrc - 1)
    f[i] = ans[3 * i]
    x[i] = ans[3 * i + 1]
    y[i] = ans[3 * i + 2]
    if i == 0:
        print("Star 1 ans: (" + str(f[i])+", "+str(x[i])+", "+str(y[i])+")")
        prffit = kepfunc.PRF2DET([f[i]],[x[i]],[y[i]],DATx,DATy,1.0,1.0,0,splineInterpolation)
        fig4 = pl.figure()
        pl.imshow(prffit, interpolation='nearest')
        pl.title('Index 0')
    else:
        print("Star 2 ans: (" + str(f[i]) +", "+ str(x[i])+", "+ str(y[i])+")")
        prffit += kepfunc.PRF2DET([f[i]],[x[i]],[y[i]],DATx,DATy,1.0,1.0,0,splineInterpolation)
        fig5 = pl.figure()
        pl.imshow(prffit, interpolation='nearest')

# prffit = kepfunc.PRF2DET(f,x,y,DATx,DATy,1.0,1.0,0.0,splineInterpolation)
    print("Star " + str(i+1) + ":")
    print("Flux:    " + str(f[i]))
    print("Coords:  ({0:.2f}, {1:.2f})".format(x[i],y[i]))
    print("")


print('Guess:    %.3e' % np.nansum((total_flux - prffit_guess) ** 2))
print('Solution: %.3e' % np.nansum((total_flux - prffit) ** 2))

vmin = 0
vmax = max(np.nanmax(total_flux), np.nanmax(prffit))
fig, ax = pl.subplots(2,3, figsize = (12,8))

im1 = ax[0,0].imshow(total_flux, interpolation = 'nearest')
ax[0,0].set_title('Data', fontsize = 22)
ax[0,1].imshow(prffit, interpolation = 'nearest')
ax[0,1].set_title('PRF Fit', fontsize = 22)
ax[0,2].imshow(total_flux - prffit, interpolation = 'nearest')
ax[0,2].set_title('Data - Model', fontsize = 22)
ax[1,0].imshow(prf)
ax[1,0].set_title('Model', fontsize = 22)

nb_args = []
nearby = data.nearby
kepmag = data.kepmag
def size(k):
  s = 1000 * 2 ** (kepmag - k)
  return s

ax[1,1].set_xlim(-0.7, nx - 0.3)
ax[1,1].set_ylim(ny - 0.3, -0.7)

numsrc = 0
src_distance = []
for source in nearby:
    if np.sqrt((source.x - source.x0)**2 + (source.y - source.y0)**2) < 10:
        numsrc += 1
        src_distance.append(np.sqrt(source.x - source.x0) ** 2 + (source.y - source.y0) ** 2)
        nb_args.extend([source.x - source.x0, source.y - source.y0, data.fpix[0],
                           data.perr[0], nsrc, splineInterpolation,106,874])
        ax[1,1].scatter(source.x - source.x0, source.y - source.y0,
                   s = size(source.kepmag),
                   c = ['r'],
                   alpha = 0.5,
                   edgecolor = 'k')
        ax[1,1].scatter(source.x - source.x0, source.y - source.y0,
                   marker = r'$%.1f$' % source.kepmag, color = 'w',
                   s = 500)

    else:
       continue

print(numsrc, src_distance)

ax[1,1].set_title('Nearby Sources')


pl.show()
