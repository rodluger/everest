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
import prffunc
from numpy import empty

# define EPIC ID
# epic = 205998445
epic = 215796924
# epic = 215915109
# epic = int(sys.argv[1])

def AddApertureContour(ax, nx, ny, aperture):

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
  base_flux = np.log10(np.nansum(fpix, axis = 0))

  # Plot the data and the aperture
  fig, ax = pl.subplots(1)
  ax.imshow(base_flux, interpolation = 'nearest', alpha = 0.75)
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

  return

# PlotPostageStamp(epic, apnum = 15)

data = everest.GetK2Data(epic)
fpix = data.fpix
_, ny, nx = fpix.shape
kepmag = data.kepmag
base_flux = fpix[0]
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

# generate PRF dimensions
prfDimY = int(ydim / cdelt1p[0])
prfDimX = int(xdim / cdelt2p[0])
PRFy0 = (np.shape(prf)[0] - prfDimY) / 2
PRFx0 = (np.shape(prf)[1] - prfDimX) / 2

DATx = np.arange(column,column+xdim)
DATy = np.arange(row,row+ydim)

# interpolate function over the PRF
splineInterpolation = scipy.interpolate.RectBivariateSpline(PRFx,PRFy,prf)


# Include only the sources that are in the aperture, and sort them
# from brightest to faintest
def findNearby():

    nearby = data.nearby
    kepmag = data.kepmag
    nearby = np.array([source for source in nearby if (source.x >= DATx[0]) and (source.x <= DATx[-1]) and (source.y >= DATy[0]) and (source.y <= DATy[-1])])
    nearby = nearby[np.argsort([source.kepmag for source in nearby])]

    return nearby

# test nearby sources to determine if they improve the fit
# by minimizing Bayesian Information Criterion (BIC)
# Testing a maximum of 3 sources
def generateGuess():

    nearby = findNearby()
    C_mag = 17
    nsrc = 0
    X = [np.nansum([i**2 for i in base_flux/errors])]
    BIC = [X[0]]

    # contruct lists for f, x, and y for each star in field
    src_distance=[];fguess=[];xguess=[];yguess=[];

    for source in nearby[:3]:

        nsrc += 1
        src_distance.append(np.sqrt(source.x - source.x0) ** 2 + (source.y - source.y0) ** 2)

        # try new target parameters
        ftry = 10**(C_mag - source.kepmag)
        xtry = source.x
        ytry = source.y
        paramstry = fguess + [ftry] + xguess + [xtry] + yguess + [ytry]

        # calculate X^2 value for set, and input into BIC
        chisq = prffunc.PRF(paramstry,DATx,DATy,base_flux,errors,nsrc,splineInterpolation,np.mean(DATx),np.mean(DATy))
        X.append(chisq)
        BIC.append(chisq + len(paramstry) * np.log(len(fpix)))

        # Append the guess
        fguess.append(ftry)
        xguess.append(xtry)
        yguess.append(ytry)

    # return the BIC and guess array for the best fit
    for i in range(len(BIC)):
        if BIC[i] == np.min(BIC[1:]):
            nsrc = i
            return BIC, fguess[:nsrc] + xguess[:nsrc] + yguess[:nsrc]
        else:
            continue

# fit PRF model to pixel data

# concatenate guess array
BIC, guess = generateGuess()
nsrc = int(len(guess) / 3)
args = (DATx,DATy,base_flux,errors,nsrc,splineInterpolation,np.mean(DATx),np.mean(DATy))

# calculate solution array based on initial guess
ans = fmin_powell(prffunc.PRF,guess,args=args,xtol=1.0e-4,ftol=1.0e-4,disp=True)

# print guess and solution arrays, and number of sources
print("\nGuess:    " + str(['%.2f' % elem for elem in guess]))
print("Solution: " + str(['%.2f' % elem for elem in ans]))
print("Number of sources fit = " + str(nsrc))

# generate the prf fit for guess parameters
prffit_guess = prffunc.PRF2DET(guess[:nsrc], guess[nsrc:2*nsrc], guess[2*nsrc:], DATx, DATy, 1.0, 1.0, 0, splineInterpolation)

# populate arrays for f, x, and i with solution
f=empty((nsrc));x=empty((nsrc));y=empty((nsrc))
for i in range(nsrc):
    f[i] = ans[i]
    x[i] = ans[nsrc+i]
    y[i] = ans[nsrc*2+i]

# LUGER: This line is crucial. For some reason the PRF fit must be calculated
# simultaneously for both sources; you can't just compute them separately and add them
# as above. Not sure why -- something to do with the way it interpolates the PRF.
prffit = prffunc.PRF2DET(f,x,y,DATx,DATy,1.0,1.0,0.0,splineInterpolation)

# LUGER: Turns out prffunc.PRF actually returns the chi squared, so we can just call that instead
print('\nGuess X^2:    %.3e' % prffunc.PRF(guess, *args))
print('Solution X^2: %.3e' % prffunc.PRF(ans, *args))


# LUGER: I added vmin and vmax to the lines below so that everything is plotted on the same scale
vmin = 0
vmax = max(np.nanmax(base_flux), np.nanmax(prffit))
fig, ax = pl.subplots(2,3, figsize = (12,8))

im1 = ax[0,0].imshow(base_flux, interpolation = 'nearest', vmin=vmin, vmax=vmax)
ax[0,0].set_title('Data', fontsize = 22)
ax[0,1].imshow(prffit, interpolation = 'nearest', vmin=vmin, vmax=vmax)
ax[0,1].set_title('PRF Fit', fontsize = 22)
ax[0,2].imshow(np.abs(base_flux - prffit), interpolation = 'nearest', vmin=vmin, vmax=vmax)
ax[0,2].set_title('Data - Model', fontsize = 22)


# LUGER: Now showing what our guess looks like as well, for reference
ax[1,0].imshow(prffit_guess, interpolation = 'nearest', vmin=vmin, vmax=vmax)
ax[1,0].set_title('PRF Guess', fontsize = 22)

# display sources in field
def size(k):
  s = 1000 * 2 ** (kepmag - k)
  return s
for source in findNearby():
    if np.sqrt((source.x - source.x0)**2 + (source.y - source.y0)**2) < 10:
        ax[1,1].scatter(source.x - source.x0, source.y - source.y0,
            s = size(source.kepmag),
            c = ['k'],
            alpha = 0.5,
            edgecolor = 'k')
        ax[1,1].scatter(source.x - source.x0, source.y - source.y0 + np.log10(size(source.kepmag)) / 2.5,
            marker = r'$%.1f$' % source.kepmag, color = 'w',
            s = 500)
    else:
        continue
ax[1,1].set_xlim(-0.7, nx - 0.3)
ax[1,1].set_ylim(ny - 0.3, -0.7)
ax[1,1].set_title('Nearby Sources')
ax[1,1].set_axis_bgcolor('darkblue')

ax[1,2].plot([(i) for i in range(len(BIC))],BIC,'r-')
ax[1,2].set_xlabel('Number of Sources')
ax[1,2].set_title('BIC vs. Sources')
ax[1,2].set_yscale('log')

pl.tight_layout()
pl.show()
