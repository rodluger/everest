#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
crowding.py
-----------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from everest.missions.k2 import GetData, Season, TargetDirectory
from everest.config import KEPPRF_DIR
from everest.math import SavGol
import os, sys, glob
import k2plr
from k2plr.config import KPLR_ROOT
import numpy as np
import matplotlib.pyplot as pl
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import fmin_powell
from scipy.ndimage import zoom
try:
  import pyfits
except ImportError:
  try:
    import astropy.io.fits as pyfits
  except ImportError:
    raise Exception('Please install the `pyfits` package.')

def PadWithZeros(vector, pad_width, iaxis, kwargs):
    '''

    '''

    vector[:pad_width[0]] = 0
    vector[-pad_width[1]:] = 0
    return vector

class Fit(object):
    '''

    '''

    def __init__(self, parent, params):
        '''

        '''

        # Convert the list ``params`` into individual variables
        i = (len(params) - 1) // 5
        self.f = params[:i]
        self.x = params[i:2 * i]
        self.y = params[2 * i:3 * i]
        self.wx = params[3 * i:4 * i]
        self.wy = params[4 * i:5 * i]
        self.a = params[-1]
        self.fit = parent.PRF2DET(self.f, self.x, self.y, self.wx, self.wy, self.a)
        self._parent = parent

    @property
    def fit0(self):
        '''
        Returns the PRF fit for *only* the desired target.

        '''

        i = np.argmax(self._parent.IDs == self._parent.ID)
        if i == 0 and self._parent.IDs[0] != self._parent.ID:
          # This shouldn't happen. If it does, the target
          # is somehow missing from the `nearby` list.
          raise ValueError("Target not found in list!")
        return self._parent.PRF2DET([self.f[i]], [self.x[i]], [self.y[i]], [self.wx[i]], [self.wy[i]], self.a)
    
    @property
    def fits(self):
        '''
        Returns a list of the PRF fits for each of the targets individually.

        '''

        return [self._parent.PRF2DET([self.f[i]], [self.x[i]], 
                                     [self.y[i]], [self.wx[i]], 
                                     [self.wy[i]], self.a) for i in range(len(self.f))]

class CrowdingTarget(object):
    '''

    '''

    def __init__(self, ID, maxsrc = 3, xtol = 1.e-4, ftol = 1.e-4, dtol = 0.01, aperture = 'k2sff_15'):
        '''

        '''

        # Target ID (EPIC number)
        self.ID = ID

        # Maximum number of sources to fit
        self.nsrc = maxsrc

        # Powell minimization params
        self.xtol = xtol
        self.ftol = ftol

        # The distance tolerance when fitting sources (in pixels)
        self.dtol = dtol

        # Name of the aperture we're using
        self.aperture = aperture

        # Get the data
        self.GetK2Data()

        # Generate the PRF and clip sources outside the aperture
        self.generatePRF()
        self.clipNearby()

        # Initialize arrays
        self.index = 100
        self.bic = []
        self.chisq = []
        self.guesses = []
        self.answers = []

    def GetK2Data(self):
        '''

        '''

        # Call `GetData()` to download the TPF
        GetData(self.ID, download_only = True)
        data = np.load(os.path.join(TargetDirectory(self.ID, Season(self.ID)), 'data.npz'))

        # Load the raw data
        fpix = data['fpix']
        ferr = data['fpix_err']
        time = data['time']
        qual = data['qual']
        pc1 = data['pc1']
        pc2 = data['pc2']
        flux = np.nansum(fpix, axis = (1,2))

        # Remove bad timestamps
        badmask = list(np.where(np.isnan(time) | np.isnan(flux) | (flux == 0))[0])
        for b in [1,2,3,4,5,6,7,8,9,11,12,13,14,16,17]:
            badmask += list(np.where(qual & 2 ** (b - 1))[0])
        badmask = np.array(list(set(badmask)))
        time = np.delete(time, badmask)
        fpix = np.delete(fpix, badmask, axis = 0)
        ferr = np.delete(ferr, badmask, axis = 0)
        flux = np.delete(flux, badmask)
        pc1 = np.delete(pc1, badmask)
        pc2 = np.delete(pc2, badmask)

        # Get the nearby targets
        nearby_dict = data['nearby']
        nearby = []
        for source in nearby_dict:
          source['flux'] = (10 ** (17. - source['mag']))
          nearby.append(type('Source', (object,), source))

        # Get the kepler magnitude
        mag = data['fitsheader'][()][0]['KEPMAG'][1]

        # Save
        self.fpix = fpix
        _, self.ny, self.nx = self.fpix.shape
        self.ferr = ferr
        self.nearby = nearby
        self.mag = mag
        self.aperture = data['apertures'][()][self.aperture]
        self.X = pc1
        self.Y = pc2

    def generatePRF(self):
        '''
        Create PRF for location of target on detector

        '''

        # read PRF header data
        client = k2plr.API()
        star = client.k2_star(self.ID)
        tpf = star.get_target_pixel_files(fetch = True)[0]
        ftpf = os.path.join(KPLR_ROOT, 'data', 'k2', 'target_pixel_files', '%d' % self.ID, tpf._filename)
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

        # Get PRF file name
        if int(module) < 10:
            prefix = 'kplr0'
        else:
            prefix = 'kplr'
        prfglob = os.path.join(KEPPRF_DIR, prefix + str(module) + '.' + str(output) + '*' + '_prf.fits')
        prffile = glob.glob(prfglob)[0]

        # create PRF matrix
        prfn = [0,0,0,0,0]
        crpix1p = np.zeros((5),dtype='float32')
        crpix2p = np.zeros((5),dtype='float32')
        crval1p = np.zeros((5),dtype='float32')
        crval2p = np.zeros((5),dtype='float32')
        cdelt1p = np.zeros((5),dtype='float32')
        cdelt2p = np.zeros((5),dtype='float32')

        # Read in the PRF
        for i in range(5):
          with pyfits.open(prffile, mode = 'readonly') as prf:
              prfn[i] = prf[i+1].data
              crpix1p[i] = prf[i+1].header['CRPIX1P']
              crpix2p[i] = prf[i+1].header['CRPIX2P']
              crval1p[i] = prf[i+1].header['CRVAL1P']
              crval2p[i] = prf[i+1].header['CRVAL2P']
              cdelt1p[i] = prf[i+1].header['CDELT1P']
              cdelt2p[i] = prf[i+1].header['CDELT2P']
        prfn = np.array(prfn)

        # PRF dimensions
        PRFx = np.arange(0.5,np.shape(prfn[0])[1]+0.5)
        PRFy = np.arange(0.5,np.shape(prfn[0])[0]+0.5)
        self.PRFx = (PRFx - np.size(PRFx) / 2) * cdelt1p[0]
        self.PRFy = (PRFy - np.size(PRFy) / 2) * cdelt2p[0]

        # Combine the PRFs
        prf = np.zeros(np.shape(prfn[0]), dtype='float32')
        prfWeight = np.zeros((5), dtype='float32')
        for i in range(5):
            prfWeight[i] = np.sqrt((column - crval1p[i])**2 + (row - crval2p[i])**2)
            if prfWeight[i] == 0.0:
                prfWeight[i] = 1.0e-6
            prf = prf + prfn[i] / prfWeight[i]
        self.prf = prf / np.nansum(prf) / cdelt1p[0] / cdelt2p[0]
        self.DATx = np.arange(column,column+xdim)
        self.DATy = np.arange(row,row+ydim)

        # Generate interpolant
        self.splineInterpolation = RectBivariateSpline(self.PRFx, self.PRFy, self.prf)

    def findSolution(self, index = 100):
        '''
        Minimize residuals to find array with best parameters.

        :param int index: The index of the timeseries to fit the PRF to. Default 100

        '''

        # This is the index of the timeseries we're going to fit
        self.index = index

        # Initialize the BIC and the CHISQ. The first term
        # corresponds to zero fit parameters (i.e., model = 0.)
        X = np.nansum((self.fpix[self.index] / self.ferr) ** 2)
        self.bic = [X]
        self.chisq = [X]
        self.guesses = [None]
        self.answers = [None]

        # return the BIC and guess array for the best fit
        for i in range(1, self.nsrc + 1):

            # Generate our guess parameter array
            f = [source.flux for source in self.nearby[:i]]
            x = [source.x for source in self.nearby[:i]]
            y = [source.y for source in self.nearby[:i]]
            wx = [1.0 for source in self.nearby[:i]]      # wx and wy are the
            wy = [1.0 for source in self.nearby[:i]]      # focusing parameters
            a = [0.]                                      # a is the angle of rotation
            guess = np.concatenate([f, x, y, wx, wy, a])

            # calculate best parameters for PRF
            answer, chisq, _, iter, funcalls, warn = fmin_powell(self.PRF, guess, xtol = self.xtol, ftol = self.ftol,
                                                                 disp = False, full_output = True)

            # Update chisq and bic
            self.chisq.append(chisq)
            self.bic.append(chisq + len(answer) * np.log(len(self.fpix)))

            # Save the guess fit and the answer fit
            self.guesses.append(Fit(self, guess))
            self.answers.append(Fit(self, answer))

    def PRF(self, params):
        '''
        Returns the residuals of the PRF fit for ``nsrc`` sources.

        '''

        # calculate PRF model binned to the detector pixel size
        fit = Fit(self, params)
        x = fit.x
        y = fit.y
        f = fit.f
        wx = fit.wx
        wy = fit.wy
        a = fit.a
        PRFfit = fit.fit

        # calculate the sum squared difference between data and model
        PRFres = np.nansum(((self.fpix[self.index] - PRFfit) / self.ferr) ** 2)

        # Prior likelihood
        x0 = np.array([source.x for source in self.nearby[:len(x)]])
        y0 = np.array([source.y for source in self.nearby[:len(y)]])
        PRFres += np.nansum(((x - x0) / self.dtol) ** 2 + ((y - y0) / self.dtol) ** 2)

        # keep the fit centered
        if max(np.abs(np.mean(self.DATx) - x[0]), np.abs(np.mean(self.DATy) - y[0])) > 10.0:
            PRFres = 1.0e300

        # Reject negative fluxes
        for elem in f:
            if elem < 0:
                PRFres = 1.0e300
        
        # Reject large/small focus factors
        if np.any(wx > 1.15) or np.any(wx < 0.85):
            PRFres = 1.0e30
        if np.any(wy > 1.15) or np.any(wy < 0.85):
            PRFres = 1.0e30
        
        return PRFres

    def PRF2DET(self, flux, OBJx, OBJy, wx, wy, a):
        '''

        '''

        # Rotation angle
        cosa = np.cos(np.radians(a))
        sina = np.sin(np.radians(a))

        # Loop over each of the sources
        PRFfit = np.zeros((np.size(self.DATy), np.size(self.DATx)))
        for i in range(len(flux)):
            FRCx,INTx = np.modf(OBJx[i])
            FRCy,INTy = np.modf(OBJy[i])
            if FRCx > 0.5:
                FRCx -= 1.0
                INTx += 1.0
            if FRCy > 0.5:
                FRCy -= 1.0
                INTy += 1.0
            FRCx = -FRCx
            FRCy = -FRCy

            # Construct model PRF in detector coordinates
            for j, y in enumerate(self.DATy):
                for k, x in enumerate(self.DATx):
                    xx = x - INTx + FRCx
                    yy = y - INTy + FRCy
                    dx = xx * cosa - yy * sina
                    dy = xx * sina + yy * cosa
                    PRFfit[j,k] = PRFfit[j,k] + self.splineInterpolation(dy * wy[i], dx * wx[i]) * flux[i]

        return PRFfit

    def clipNearby(self):
        '''

        '''

        self.nearby = np.array([source for source in self.nearby if
                               (source.x >= self.DATx[0]) and
                               (source.x <= self.DATx[-1]) and
                               (source.y >= self.DATy[0]) and
                               (source.y <= self.DATy[-1])])
        self.nearby = self.nearby[np.argsort([source.mag for source in self.nearby])]
        self.nearby = self.nearby[:self.nsrc]
        self.nsrc = len(self.nearby)
        self.IDs = [source.ID for source in self.nearby]

    def findCrowding(self):
        '''
        Returns a crowding parameter defined by
            C = F_star / F_total
        '''

        # crowding parameter for entire postage stamp
        self.c_postage = np.nansum(self.answers[-1].fit0) / np.nansum(self.answers[-1].fit)

        # add aperture
        contour = np.zeros((self.ny, self.nx))
        contour[np.where(self.aperture)] = 1

        ap_total = []
        ap_target = []

        # crowding parameter for each pixel
        self.c_pixel = np.zeros((self.ny,self.nx))

        for i in range(self.ny):
            for j in range(self.nx):

                self.c_pixel[i][j] = self.answers[-1].fit0[i][j] / self.answers[-1].fit[i][j]

                if contour[i][j] == 1:
                    ap_total.append(self.answers[-1].fit[i][j])
                    ap_target.append(self.answers[-1].fit0[i][j])
        
        import pdb; pdb.set_trace()
        
        # crowding parameter for aperture
        self.c_aperture = np.nansum(ap_target) / np.nansum(ap_total)

    def plot(self):
        '''

        '''

        rdbu = pl.get_cmap('RdBu_r')

        fig, ax = pl.subplots(2, 2, figsize = (8, 8))
        vmax = np.max([np.nanmax(self.fpix[self.index])] + [np.nanmax(a.fit) for a in self.answers[1:]])
        vmin = np.min([np.nanmin(self.fpix[self.index])] + [np.nanmin(a.fit) for a in self.answers[1:]])
        vmax = max(vmax, -vmin)
        vmin = -vmax

        # Show the data
        ax[0,0].imshow(self.fpix[self.index], interpolation = 'nearest', vmin = vmin, vmax = vmax, cmap = rdbu)
        ax[0,0].set_title('Data')

        # Show the fit
        ax[0,1].imshow(self.answers[-1].fit, interpolation = 'nearest', vmin = vmin, vmax = vmax, cmap = rdbu)
        ax[0,1].set_title('Fit')

        ax[0,1].annotate(r'$\mathrm{Sources}: %i$' % self.nsrc,
                         xy = (0.025, 0.95), xycoords = 'axes fraction',
                         ha = 'left', va = 'top', color = 'k', fontsize = 12)

        # Add aperture contour
        contour = np.zeros((self.ny, self.nx))
        contour[np.where(self.aperture)] = 1
        contour = np.lib.pad(contour, 1, PadWithZeros)
        highres = zoom(contour, 100, order = 0, mode='nearest')
        extent = np.array([-1, self.nx, -1, self.ny])

        # Show aperture and set bounds
        for i in range(2):
            for j in range(2):
                ax[i,j].contour(highres, levels=[0.5], extent=extent, origin='lower', colors='r', linewidths=1)
                ax[i,j].set_xlim(-0.5, self.nx - 0.5)
                ax[i,j].set_ylim(self.ny - 0.5, -0.5)

        for n in range(self.nsrc):

            # Display the catalog positions
            ax[0,1].plot(self.nearby[n].x - self.nearby[n].x0,
                         self.nearby[n].y - self.nearby[n].y0,
                         'ko', markeredgecolor = 'none', alpha = 0.5,
                         markersize = 10)

            # Display the solution positions
            ax[0,1].plot(self.answers[n+1].x - self.nearby[n].x0,
                         self.answers[n+1].y - self.nearby[n].y0,
                         'ro', markeredgecolor = 'none',
                         markersize = 3)

        # Calculate the fractional error in the fit
        err = np.sqrt(np.nansum((self.fpix[self.index] - self.answers[-1].fit) ** 2) / np.nansum(self.fpix[self.index] ** 2))

        # Display fit info
        ax[1,0].annotate(r'$\log(\chi^2) = %.3f$' % np.log10(self.chisq[-1]),
                         xy = (0.025, 0.95), xycoords = 'axes fraction',
                         ha = 'left', va = 'top', color = 'k', fontsize = 12)
        ax[1,0].annotate(r'$\mathrm{ERROR} = %.1f$' % (100 * err) + r'$\%$',
                         xy = (0.975, 0.05), xycoords = 'axes fraction',
                         ha = 'right', va = 'bottom', color = 'k', fontsize = 12)

        # Show the residuals
        ax[1,0].imshow(self.fpix[self.index] - self.answers[-1].fit, interpolation = 'nearest', vmin = vmin, vmax = vmax, cmap = rdbu)
        ax[1,0].set_title('Residuals')

        # Show the crowding
        ax[1,1].imshow(self.c_pixel, interpolation='nearest', cmap=rdbu)
        ax[1,1].set_title('Crowding')
        ax[1,1].annotate(r'$C_{total} = %.3f$' % self.c_postage,
                         xy = (0.025, 0.95), xycoords = 'axes fraction',
                         ha = 'left', va = 'top', color = 'k', fontsize = 12)
        ax[1,1].annotate(r'$C_{aperture} = %.3f$' % self.c_aperture,
                         xy = (0.975, 0.05), xycoords = 'axes fraction',
                         ha = 'right', va = 'bottom', color = 'k', fontsize = 12)

        pl.show()

c = CrowdingTarget(215796924)
c.findSolution(3000)

# Find crowding parameter for postage stamp, aperture, and pixel
c.findCrowding()

c.plot()
