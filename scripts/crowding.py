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
import prffunc
from numpy import empty

# define EPIC ID
# epic = 205998445
epic = 215796924
# epic = 215915109
# epic = int(sys.argv[1])



class CrowdingTarget(object):

    def __init__(self, epic):
        self.data = everest.GetK2Data(epic)
        self.nearby = self.data.nearby
        self.fpix = self.data.fpix
        self.kepmag = self.data.kepmag
        self._, self.ny, self.nx = self.fpix.shape
        self.base_flux = self.fpix[0]
        self.mean_flux = np.sum(self.fpix[n] for n in range(len(self.fpix))) / len(self.fpix)
        self.errors = self.data.perr[0]
        self.mean_errs = np.sum(self.data.perr[n] for n in range(len(self.data.perr))) / len(self.data.perr)
        self.BIC = []
        self.nsrc = 1
        self.epic = epic

        # current maximum number of nsrc being fit
        self.nsrc_fit = 3

        # current Kepler PRF directory
        self.prfdir = '/Users/nks1994/Documents/Research/KeplerPRF'


    # draw contour around aperture for plotting postage stamp
    def addApertureContour(self,ax, nx, ny, aperture):

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

    # plot region around target star with nearby neighbors marked
    def plotPostageStamp(self, apnum = 15):

        self.aperture = self.data.apertures[apnum]

        # Plot the data and the aperture
        fig, ax = pl.subplots(1)
        ax.imshow(self.base_flux, interpolation = 'nearest', alpha = 0.75)
        self.addApertureContour(ax, self.nx, self.ny, self.aperture)

        # Crop the image
        ax.set_xlim(-0.7, self.nx - 0.3)
        ax.set_ylim(-0.7, self.ny - 0.3)

        # Overplot nearby sources
        neighbors = []
        def size(k):
            s = 1000 * 2 ** (self.kepmag - k)
            return s

        for source in self.nearby:
            ax.scatter(source.x - source.x0, source.y - source.y0,
                       s = size(source.kepmag),
                       c = ['g' if source.epic == self.epic else 'r'],
                       alpha = 0.5,
                       edgecolor = 'k')
            ax.scatter(source.x - source.x0, source.y - source.y0,
                       marker = r'$%.1f$' % source.kepmag, color = 'w',
                       s = 500)
            ax.set_xticklabels([])
            ax.set_yticklabels([])

            pl.suptitle('EPIC %d' % self.epic, fontsize = 25, y = 0.98)

        return

    # create PRF for location of target on detector
    def generatePRF(self):

        # read PRF header data
        logfile = 'test.log'
        verbose = True
        client = k2plr.API()
        star = client.k2_star(self.epic)
        tpf = star.get_target_pixel_files(fetch = True)[0]
        ftpf = os.path.join(KPLR_ROOT, 'data', 'k2', 'target_pixel_files', '%d' % self.epic, tpf._filename)
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
            prfglob = self.prfdir + '/' + prefix + str(module) + '.' + str(output) + '*' + '_prf.fits'
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
                = prffunc.readPRFimage(prffile,i+1,logfile,verbose)
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

        # generate PRF dimensions
        prfDimY = int(ydim / cdelt1p[0])
        prfDimX = int(xdim / cdelt2p[0])
        PRFy0 = (np.shape(prf)[0] - prfDimY) / 2
        PRFx0 = (np.shape(prf)[1] - prfDimX) / 2

        DATx = np.arange(column,column+xdim)
        DATy = np.arange(row,row+ydim)

        # returns data dimensions, PRF dimensions, and the prf at the target location
        return DATx,DATy,PRFx,PRFy,prf

    # interpolate function over the PRF
    def interpolate(self, PRFx, PRFy, prf):

        return scipy.interpolate.RectBivariateSpline(PRFx,PRFy,prf)

    # include only the sources that are in the aperture
    # and sort them from brightest to faintest
    def findNearby(self,DATx,DATy):

        data = self.data

        nearby = data.nearby
        kepmag = data.kepmag
        nearby = np.array([source for source in nearby if (source.x >= DATx[0])
                          and (source.x <= DATx[-1]) and (source.y >= DATy[0]) and (source.y <= DATy[-1])])
        nearby = nearby[np.argsort([source.kepmag for source in nearby])]

        return nearby

    # test nearby sources to determine if they improve the fit
    # by minimizing Bayesian Information Criterion (BIC)
    # Testing a maximum of 3 sources
    def generateGuess(self,PRFx,PRFy,prf,DATx,DATy,splineInterpolation,nearby):

        C_mag = 17
        nsrc = 0

        # contruct lists for f, x, and y for each star in field
        src_distance=[];fguess=[];xguess=[];yguess=[];

        for source in nearby[:self.nsrc_fit]:

            nsrc += 1
            src_distance.append(np.sqrt(source.x - source.x0) ** 2 + (source.y - source.y0) ** 2)

            # Append the guess
            fguess.append(10**(C_mag - source.kepmag))
            xguess.append(source.x)
            yguess.append(source.y)

        return fguess + xguess + yguess

    # minimize residuals to find array with best parameters
    def findSolution(self, guess, DATx, DATy,splineInterpolation):

        X = [np.nansum([i**2 for i in self.base_flux/self.errors])]
        self.BIC = [X[0]]
        self.nsrc = int(len(guess) / 3)
        nsrc = self.nsrc
        args = (DATx,DATy,self.base_flux,self.errors,nsrc,splineInterpolation,np.mean(DATx),np.mean(DATy))
        f = guess[:nsrc]; x = guess[nsrc:2*nsrc]; y = guess[2*nsrc:];

        # initial guess for first target
        paramstry = np.concatenate((f[:1],x[:1],y[:1]),axis=0);status=1;

        # return the BIC and guess array for the best fit
        for i in range(nsrc):

            # create temporary parameters and arguments
            if status == 1:
                status = 0

            # for nsrc > 1, add the new guess parameters to the solved for parameters
            else:
                ftry = np.concatenate((paramstry[:i], [f[i]]),axis=0)
                xtry = np.concatenate((paramstry[i:2*i], [x[i]]),axis=0)
                ytry = np.concatenate((paramstry[2*i:], [y[i]]),axis=0)

                paramstry = np.concatenate((ftry,xtry,ytry),axis=0)


            args = (DATx,DATy,self.base_flux,self.errors,i+1,splineInterpolation,np.mean(DATx),np.mean(DATy))

            # calculate best parameters for PRF
            ans = fmin_powell(prffunc.PRF,paramstry,args=args,xtol=1.0e-4,ftol=1.0e-4,disp=False)
            paramstry = ans

            # calculate X^2 value for set, and input into BIC
            chisq = prffunc.PRF(ans,DATx,DATy,self.base_flux,self.errors,i+1,splineInterpolation,np.mean(DATx),np.mean(DATy))
            X.append(chisq)

            # BIC = prffunc.PRF + k + ln(n)
            self.BIC.append(chisq + len(paramstry) * np.log(len(self.fpix)))

        for i in range(nsrc):

            # set nsrc to index (>0) of minimum BIC
            if self.BIC[i] == np.min(self.BIC[1:]):
                self.nsrc = int(i)

                ans = np.concatenate((ans[:i], ans[nsrc:nsrc+i], ans[2*nsrc:2*nsrc+i]),axis=0)
                guess = np.concatenate((guess[:i], guess[nsrc:nsrc+i], guess[2*nsrc:2*nsrc+i]),axis=0)
            else:
                continue

        nsrc = self.nsrc
        # create final arguments and parameters for BIC-minimized PRF fit
        args = (DATx,DATy,self.base_flux,self.errors,nsrc,splineInterpolation,np.mean(DATx),np.mean(DATy))

        # calculate solution array based on initial guess

        # print guess and solution arrays, and number of sources
        print("\nGuess:    " + str(['%.2f' % elem for elem in guess]))
        print("Solution: " + str(['%.2f' % elem for elem in ans]))
        print("Number of sources fit = " + str(nsrc))

        print('\nGuess X^2:    %.3e' % prffunc.PRF(guess, *args))
        print('Solution X^2: %.3e\n' % prffunc.PRF(ans, *args))

        return ans

    # create the guess prf fit
    def createGuessFit(self,guess,DATx,DATy,splineInterpolation):

        nsrc = self.nsrc
        # generate the prf fit for guess parameters
        return prffunc.PRF2DET(guess[:nsrc], guess[nsrc:2*nsrc], guess[2*nsrc:], DATx, DATy, 1.0, 1.0, 0, splineInterpolation)

    # create the prf fit for the solution array
    def createFit(self,ans,DATx,DATy,splineInterpolation):

        nsrc = self.nsrc

        f=empty((nsrc));x=empty((nsrc));y=empty((nsrc))
        for i in range(nsrc):
            f[i] = ans[i]
            x[i] = ans[nsrc+i]
            y[i] = ans[nsrc*2+i]

        return prffunc.PRF2DET(f,x,y,DATx,DATy,1.0,1.0,0.0,splineInterpolation)

    def timeSeries(self):

        pass

    # plot data, model, residuals, prf model, nearby neighbors, and the BIC
    def plotResults(self,prffit,guessfit,nearby):

        # LUGER: I added vmin and vmax to the lines below so that everything is plotted on the same scale

        vmin = 0
        vmax = max(np.nanmax(self.base_flux), np.nanmax(prffit))
        fig, ax = pl.subplots(2,3, figsize = (12,8))

        im1 = ax[0,0].imshow(self.base_flux, interpolation = 'nearest', vmin=vmin, vmax=vmax)
        ax[0,0].set_title('Data', fontsize = 22)
        ax[0,1].imshow(prffit, interpolation = 'nearest', vmin=vmin, vmax=vmax)
        ax[0,1].set_title('PRF Fit', fontsize = 22)
        ax[0,2].imshow(np.abs(self.base_flux - prffit), interpolation = 'nearest', vmin=vmin, vmax=vmax)
        ax[0,2].set_title('Data - Model', fontsize = 22)


        # LUGER: Now showing what our guess looks like as well, for reference
        ax[1,0].imshow(guessfit, interpolation = 'nearest', vmin=vmin, vmax=vmax)
        ax[1,0].set_title('PRF Guess', fontsize = 22)

        # display sources in field
        def size(k):
          s = 1000 * 2 ** (self.kepmag - k)
          return s
        for source in nearby:
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
        ax[1,1].set_xlim(-0.7, self.nx - 0.3)
        ax[1,1].set_ylim(self.ny - 0.3, -0.7)
        ax[1,1].set_title('Nearby Sources')
        ax[1,1].set_axis_bgcolor('darkblue')

        ax[1,2].plot([(i) for i in range(len(self.BIC))],self.BIC,'r-')
        ax[1,2].set_xlabel('Number of Sources')
        ax[1,2].set_title('BIC vs. Sources')
        ax[1,2].set_yscale('log')

        pl.tight_layout()
        pl.show()

    # calls functions
    def runCrowding(self):

        DATx,DATy,PRFx,PRFy,prf = self.generatePRF()
        splineInterpolation = self.interpolate(PRFx,PRFy,prf)
        nearby = self.findNearby(DATx,DATy)

        # concatenate guess array
        guess = self.generateGuess(PRFx,PRFy,prf,DATx,DATy,splineInterpolation,nearby)
        ans = self.findSolution(guess, DATx, DATy,splineInterpolation)

        # generate the prf fit for guess parameters
        guessfit = self.createGuessFit(guess,DATx,DATy,splineInterpolation)
        prffit = self.createFit(ans,DATx,DATy,splineInterpolation)

        self.plotResults(prffit,guessfit,nearby)

CrowdingTarget(epic).runCrowding()
