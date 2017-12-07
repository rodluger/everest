#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`user.py` - User Python routines
----------------------------------------

This is the gateway to the :py:obj:`everest` catalog, containing
all of the user-facing code.

- :py:class:`Everest` is the main user-facing class for
  interfacing with the catalog
- :py:func:`DVS` downloads and plots the data validation
  summary for a given target

Instantiating an :py:class:`Everest` class automatically downloads
the light curve from the online MAST catalog. So, to get started,
all you need to do is run

.. code-block :: python

   import everest
   star = everest.Everest(201367065)

'''

from __future__ import division, print_function, absolute_import, \
     unicode_literals
from . import __version__ as EVEREST_VERSION
from . import missions
from .basecamp import Basecamp
from .detrender import pPLD
from .gp import GetCovariance, GP
from .config import QUALITY_BAD, QUALITY_NAN, QUALITY_OUT, QUALITY_REC, \
     QUALITY_TRN, EVEREST_DEV, EVEREST_FITS, EVEREST_MAJOR_MINOR
from .utils import InitLog, Formatter
import george
import os
import sys
import platform
import numpy as np
import matplotlib.pyplot as pl
try:
    import pyfits
except ImportError:
    try:
        import astropy.io.fits as pyfits
    except ImportError:
        raise Exception('Please install the `pyfits` package.')
import subprocess
import six
from six.moves import urllib
from tempfile import NamedTemporaryFile
import shutil
from distutils.version import LooseVersion
import k2plr
k2plr_client = k2plr.API()
import logging
log = logging.getLogger(__name__)


def Search(ID, mission='k2'):
    """Why is my target not in the EVEREST database?"""
    # Only K2 supported for now
    assert mission == 'k2', "Only the K2 mission is supported for now."
    print("Searching for target %d..." % ID)

    # First check if it is in the database
    season = missions.k2.Season(ID)
    if season in [91, 92, [91, 92]]:
        print("Campaign 9 is currently not part of the EVEREST catalog.")
        return
    elif season == 101:
        print("The first half of campaign 10 is not currently part of " +
              "the EVEREST catalog.")
        return
    elif season is not None:
        print("Target is in campaign %d of the EVEREST catalog." % season)
        return

    # Get the kplr object
    star = k2plr_client.k2_star(ID)

    # First check if this is a star
    if star.objtype.lower() != "star":
        print("Target is of type %s, not STAR, " % star.objtype +
              "and is therefore not included in the EVEREST catalog.")
        return

    # Let's try to download the pixel data and see what happens
    try:
        tpf = star.get_target_pixel_files()
    except:
        print("Unable to download the raw pixel files for this target.")
        return
    if len(tpf) == 0:
        print("Raw pixel files are not available for this target. Looks like " +
              "data may not have been collected for it.")
        return

    # Perhaps it's in a campaign we haven't gotten to yet
    if tpf[0].sci_campaign not in missions.k2.SEASONS:
        print("Targets for campaign %d are not yet available."
              % tpf[0].sci_campaign)
        return

    # Let's try to download the K2SFF data
    try:
        k2sff = k2plr.K2SFF(ID)
    except:
        print("Error downloading the K2SFF light curve for this target. " +
              "Currently, EVEREST uses the K2SFF apertures to perform " +
              "photometry. This is likely to change in the next version.")
        return

    # Let's try to get the aperture
    try:
        assert np.count_nonzero(k2sff.apertures[15]), "Invalid aperture."
    except:
        print("Unable to retrieve the K2SFF aperture for this target. " +
              "Currently, EVEREST uses the K2SFF apertures to perform " +
              "photometry. This is likely to change in the next version.")
        return

    # Perhaps the star is *super* saturated and we didn't bother
    # de-trending it?
    if star.kp < 8:
        print("Target has Kp = %.1f and is too saturated " +
              "for proper de-trending with EVEREST.")
        return

    # I'm out of ideas
    print("I'm not sure why this target isn't in the EVEREST catalog." +
          "You can try de-trending it yourself:")
    print("http://faculty.washington.edu/rodluger/everest/pipeline.html")
    return

def DownloadFile(ID, season=None, mission='k2', cadence='lc',
                 filename=None, clobber=False):
    '''
    Download a given :py:mod:`everest` file from MAST.

    :param str mission: The mission name. Default `k2`
    :param str cadence: The light curve cadence. Default `lc`
    :param str filename: The name of the file to download. Default \
           :py:obj:`None`, in which case the default \
           FITS file is retrieved.
    :param bool clobber: If :py:obj:`True`, download and overwrite \
           existing files. Default :py:obj:`False`

    '''

    # Get season
    if season is None:
        season = getattr(missions, mission).Season(ID)
    if hasattr(season, '__len__'):
        raise AttributeError(
            "Please choose a `season` for this target: %s." % season)
    if season is None:
        if getattr(missions, mission).ISTARGET(ID):
            raise ValueError('Target not found in local database. ' +
                             'Run `everest.Search(%d)` for more information.'
                             % ID)
        else:
            raise ValueError('Invalid target ID.')
    path = getattr(missions, mission).TargetDirectory(ID, season)
    relpath = getattr(missions, mission).TargetDirectory(
        ID, season, relative=True)
    if filename is None:
        filename = getattr(missions, mission).FITSFile(ID, season, cadence)

    # Check if file exists
    if not os.path.exists(path):
        os.makedirs(path)
    elif os.path.exists(os.path.join(path, filename)) and not clobber:
        log.info('Found cached file.')
        return os.path.join(path, filename)

    # Get file URL
    log.info('Downloading the file...')
    fitsurl = getattr(missions, mission).FITSUrl(ID, season)
    if not fitsurl.endswith('/'):
        fitsurl += '/'

    # Download the data
    r = urllib.request.Request(fitsurl + filename)
    try:
        handler = urllib.request.urlopen(r)
        code = handler.getcode()
    except (urllib.error.HTTPError, urllib.error.URLError):
        code = 0
    if int(code) == 200:

        # Read the data
        data = handler.read()

        # Atomically save to disk
        f = NamedTemporaryFile("wb", delete=False)
        f.write(data)
        f.flush()
        os.fsync(f.fileno())
        f.close()
        shutil.move(f.name, os.path.join(path, filename))

    else:

        # Something went wrong!
        log.error("Error code {0} for URL '{1}'".format(
            code, fitsurl + filename))

        # If the files can be accessed by `ssh`, let's try that
        # (development version only!)
        if EVEREST_FITS is None:
            raise Exception("Unable to locate the file.")

        # Get the url
        inpath = os.path.join(EVEREST_FITS, relpath, filename)
        outpath = os.path.join(path, filename)

        # Download the data
        log.info("Accessing file via `scp`...")
        subprocess.call(['scp', inpath, outpath])

    # Success?
    if os.path.exists(os.path.join(path, filename)):
        return os.path.join(path, filename)
    else:
        raise Exception("Unable to download the file." +
                        "Run `everest.Search(%d)` to troubleshoot." % ID)


def DVS(ID, season=None, mission='k2', clobber=False,
        cadence='lc', model='nPLD'):
    '''
    Show the data validation summary (DVS) for a given target.

    :param str mission: The mission name. Default `k2`
    :param str cadence: The light curve cadence. Default `lc`
    :param bool clobber: If :py:obj:`True`, download and overwrite \
           existing files. Default :py:obj:`False`

    '''

    # Get season
    if season is None:
        season = getattr(missions, mission).Season(ID)
    if hasattr(season, '__len__'):
        raise AttributeError(
            "Please choose a `season` for this target: %s." % season)

    # Get file name
    if model == 'nPLD':
        filename = getattr(missions, mission).DVSFile(ID, season, cadence)
    else:
        if cadence == 'sc':
            filename = model + '.sc.pdf'
        else:
            filename = model + '.pdf'

    file = DownloadFile(ID, season=season,
                        mission=mission,
                        filename=filename,
                        clobber=clobber)

    try:
        if platform.system().lower().startswith('darwin'):
            subprocess.call(['open', file])
        elif os.name == 'nt':
            os.startfile(file)
        elif os.name == 'posix':
            subprocess.call(['xdg-open', file])
        else:
            raise Exception("")
    except:
        log.info("Unable to open the pdf. Try opening it manually:")
        log.info(file)


class Everest(Basecamp):
    '''
    The main user-accessible :py:mod:`everest` class for interfacing with the
    light curves stored on MAST. Instantiating this class downloads the current
    :py:mod:`everest` FITS file for the requested target and populates the
    class instance with the light curve data and attributes. Many of the
    methods are inherited from :py:class:`everest.Basecamp`.

    :param int ID: The target ID. For `k2`, this is the `EPIC` \
           number of the star.
    :param str mission: The mission name. Default `k2`
    :param bool quiet: Suppress :py:obj:`stdout` messages? \
           Default :py:obj:`False`
    :param str cadence: The light curve cadence. Default `lc`
    :param bool clobber: If :py:obj:`True`, download and overwrite existing \
           files. Default :py:obj:`False`

    '''

    def __init__(self, ID, season=None, mission='k2', quiet=False,
                 clobber=False, cadence='lc', **kwargs):
        '''

        '''

        # Read kwargs
        self.ID = ID
        self.mission = mission
        self.clobber = clobber
        if season is not None:
            self._season = season

        # Initialize preliminary logging
        if not quiet:
            screen_level = logging.DEBUG
        else:
            screen_level = logging.CRITICAL
        InitLog(None, logging.DEBUG, screen_level, False)

        # Check the cadence
        if cadence not in ['lc', 'sc']:
            raise ValueError("Invalid cadence selected.")
        self.cadence = cadence

        # Download the FITS file if necessary
        self.fitsfile = DownloadFile(
            ID, season=season, mission=mission, clobber=clobber,
            cadence=cadence)
        self.model_name = pyfits.getheader(self.fitsfile, 1)['MODEL']
        self._weights = None

        # Check the pipeline version. Do we need to upgrade?
        subversion = pyfits.getheader(self.fitsfile, 1).get('SUBVER', None)
        if subversion is not None:
            if LooseVersion(subversion) > LooseVersion(EVEREST_VERSION):
                raise Exception("Desired light curve was generated with " +
                                "EVEREST version %s, but current version " +
                                "is %s.\n" % (subversion, EVEREST_VERSION) +
                                "Please upgrade EVEREST by running " +
                                "`pip install everest-pipeline --upgrade`.")

        # Load the FITS file
        self.load_fits()

    def __repr__(self):
        '''

        '''

        return "<everest.Everest(%d)>" % self.ID

    @property
    def name(self):
        '''
        Returns the name of the :py:mod:`everest` model used
        to generate this light curve.

        '''

        return self.model_name

    def reset(self):
        '''
        Re-loads the FITS file from disk.

        '''

        self.load_fits()
        self._weights = None

    def compute(self):
        '''
        Re-compute the :py:mod:`everest` model for the given
        value of :py:obj:`lambda`.
        For long cadence `k2` light curves, this should take several
        seconds. For short cadence `k2` light curves, it may take a
        few minutes. Note that this is a simple wrapper around
        :py:func:`everest.Basecamp.compute`.

        '''

        # If we're doing iterative PLD, get the normalization
        if self.model_name == 'iPLD':
            self._get_norm()

        # Compute as usual
        super(Everest, self).compute()

        # Make NaN cadences NaNs
        self.flux[self.nanmask] = np.nan

    def _get_norm(self):
        '''
        Computes the PLD flux normalization array.

        ..note :: `iPLD` model **only**.

        '''

        log.info('Computing the PLD normalization...')

        # Loop over all chunks
        mod = [None for b in self.breakpoints]
        for b, brkpt in enumerate(self.breakpoints):

            # Unmasked chunk
            c = self.get_chunk(b)

            # Masked chunk (original mask plus user transit mask)
            inds = np.array(
                list(set(np.concatenate([self.transitmask,
                                         self.recmask]))), dtype=int)
            M = np.delete(np.arange(len(self.time)), inds, axis=0)
            if b > 0:
                m = M[(M > self.breakpoints[b - 1] - self.bpad)
                      & (M <= self.breakpoints[b] + self.bpad)]
            else:
                m = M[M <= self.breakpoints[b] + self.bpad]

            # This block of the masked covariance matrix
            mK = GetCovariance(self.kernel, self.kernel_params,
                               self.time[m], self.fraw_err[m])

            # Get median
            med = np.nanmedian(self.fraw[m])

            # Normalize the flux
            f = self.fraw[m] - med

            # The X^2 matrices
            A = np.zeros((len(m), len(m)))
            B = np.zeros((len(c), len(m)))

            # Loop over all orders
            for n in range(self.pld_order):
                XM = self.X(n, m)
                XC = self.X(n, c)
                A += self.reclam[b][n] * np.dot(XM, XM.T)
                B += self.reclam[b][n] * np.dot(XC, XM.T)
                del XM, XC

            W = np.linalg.solve(mK + A, f)
            mod[b] = np.dot(B, W)
            del A, B, W

        # Join the chunks after applying the correct offset
        if len(mod) > 1:

            # First chunk
            model = mod[0][:-self.bpad]

            # Center chunks
            for m in mod[1:-1]:
                offset = model[-1] - m[self.bpad - 1]
                model = np.concatenate(
                    [model, m[self.bpad:-self.bpad] + offset])

            # Last chunk
            offset = model[-1] - mod[-1][self.bpad - 1]
            model = np.concatenate([model, mod[-1][self.bpad:] + offset])

        else:

            model = mod[0]

        # Subtract the global median
        model -= np.nanmedian(model)

        # Save the norm
        self._norm = self.fraw - model

    def load_fits(self):
        '''
        Load the FITS file from disk and populate the
        class instance with its data.

        '''

        log.info("Loading FITS file for %d." % (self.ID))
        with pyfits.open(self.fitsfile) as f:

            # Params and long cadence data
            self.loaded = True
            self.is_parent = False
            try:
                self.X1N = f[2].data['X1N']
            except KeyError:
                self.X1N = None
            self.aperture = f[3].data
            self.aperture_name = f[1].header['APNAME']
            try:
                self.bkg = f[1].data['BKG']
            except KeyError:
                self.bkg = 0.
            self.bpad = f[1].header['BPAD']
            self.cbv_minstars = []
            self.cbv_num = f[1].header.get('CBVNUM', 1)
            self.cbv_niter = f[1].header['CBVNITER']
            self.cbv_win = f[1].header['CBVWIN']
            self.cbv_order = f[1].header['CBVORD']
            self.cadn = f[1].data['CADN']
            self.cdivs = f[1].header['CDIVS']
            self.cdpp = f[1].header['CDPP']
            self.cdppr = f[1].header['CDPPR']
            self.cdppv = f[1].header['CDPPV']
            self.cdppg = f[1].header['CDPPG']
            self.cv_min = f[1].header['CVMIN']
            self.fpix = f[2].data['FPIX']
            self.pixel_images = [f[4].data['STAMP1'],
                                 f[4].data['STAMP2'], f[4].data['STAMP3']]
            self.fraw = f[1].data['FRAW']
            self.fraw_err = f[1].data['FRAW_ERR']
            self.giter = f[1].header['GITER']
            self.gmaxf = f[1].header.get('GMAXF', 200)
            self.gp_factor = f[1].header['GPFACTOR']
            try:
                self.hires = f[5].data
            except:
                self.hires = None
            self.kernel_params = np.array([f[1].header['GPWHITE'],
                                           f[1].header['GPRED'],
                                           f[1].header['GPTAU']])
            try:
                self.kernel = f[1].header['KERNEL']
                self.kernel_params = np.append(
                    self.kernel_params,
                    [f[1].header['GPGAMMA'],
                     f[1].header['GPPER']])
            except KeyError:
                self.kernel = 'Basic'
            self.pld_order = f[1].header['PLDORDER']
            self.lam_idx = self.pld_order
            self.leps = f[1].header['LEPS']
            self.mag = f[0].header['KEPMAG']
            self.max_pixels = f[1].header['MAXPIX']
            self.model = self.fraw - f[1].data['FLUX']
            self.nearby = []
            for i in range(99):
                try:
                    ID = f[1].header['NRBY%02dID' % (i + 1)]
                    x = f[1].header['NRBY%02dX' % (i + 1)]
                    y = f[1].header['NRBY%02dY' % (i + 1)]
                    mag = f[1].header['NRBY%02dM' % (i + 1)]
                    x0 = f[1].header['NRBY%02dX0' % (i + 1)]
                    y0 = f[1].header['NRBY%02dY0' % (i + 1)]
                    self.nearby.append(
                        {'ID': ID, 'x': x, 'y': y,
                         'mag': mag, 'x0': x0, 'y0': y0})
                except KeyError:
                    break
            self.neighbors = []
            for c in range(99):
                try:
                    self.neighbors.append(f[1].header['NEIGH%02d' % (c + 1)])
                except KeyError:
                    break
            self.oiter = f[1].header['OITER']
            self.optimize_gp = f[1].header['OPTGP']
            self.osigma = f[1].header['OSIGMA']
            self.planets = []
            for i in range(99):
                try:
                    t0 = f[1].header['P%02dT0' % (i + 1)]
                    per = f[1].header['P%02dPER' % (i + 1)]
                    dur = f[1].header['P%02dDUR' % (i + 1)]
                    self.planets.append((t0, per, dur))
                except KeyError:
                    break
            self.quality = f[1].data['QUALITY']
            self.saturated = f[1].header['SATUR']
            self.saturation_tolerance = f[1].header['SATTOL']
            self.time = f[1].data['TIME']
            self._norm = np.array(self.fraw)

            # Chunk arrays
            self.breakpoints = []
            self.cdpp_arr = []
            self.cdppv_arr = []
            self.cdppr_arr = []
            for c in range(99):
                try:
                    self.breakpoints.append(f[1].header['BRKPT%02d' % (c + 1)])
                    self.cdpp_arr.append(f[1].header['CDPP%02d' % (c + 1)])
                    self.cdppr_arr.append(f[1].header['CDPPR%02d' % (c + 1)])
                    self.cdppv_arr.append(f[1].header['CDPPV%02d' % (c + 1)])
                except KeyError:
                    break
            self.lam = [[f[1].header['LAMB%02d%02d' % (c + 1, o + 1)]
                        for o in range(self.pld_order)]
                        for c in range(len(self.breakpoints))]
            if self.model_name == 'iPLD':
                self.reclam = [[f[1].header['RECL%02d%02d' % (c + 1, o + 1)]
                               for o in range(self.pld_order)]
                               for c in range(len(self.breakpoints))]

            # Masks
            self.badmask = np.where(self.quality & 2 ** (QUALITY_BAD - 1))[0]
            self.nanmask = np.where(self.quality & 2 ** (QUALITY_NAN - 1))[0]
            self.outmask = np.where(self.quality & 2 ** (QUALITY_OUT - 1))[0]
            self.recmask = np.where(self.quality & 2 ** (QUALITY_REC - 1))[0]
            self.transitmask = np.where(
                self.quality & 2 ** (QUALITY_TRN - 1))[0]

            # CBVs
            self.XCBV = np.empty((len(self.time), 0))
            for i in range(99):
                try:
                    self.XCBV = np.hstack(
                        [self.XCBV,
                         f[1].data['CBV%02d' % (i + 1)].reshape(-1, 1)])
                except KeyError:
                    break

        # These are not stored in the fits file; we don't need them
        self.saturated_aperture_name = None
        self.apertures = None
        self.Xpos = None
        self.Ypos = None
        self.fpix_err = None
        self.parent_model = None
        self.lambda_arr = None
        self.meta = None
        self._transit_model = None
        self.transit_depth = None

    def plot_aperture(self, show=True):
        '''
        Plot sample postage stamps for the target with the aperture
        outline marked, as well as a high-res target image (if available).

        :param bool show: Show the plot or return the `(fig, ax)` instance? \
               Default :py:obj:`True`

        '''

        # Set up the axes
        fig, ax = pl.subplots(2, 2, figsize=(6, 8))
        fig.subplots_adjust(top=0.975, bottom=0.025, left=0.05,
                            right=0.95, hspace=0.05, wspace=0.05)
        ax = ax.flatten()
        fig.canvas.set_window_title(
            '%s %d' % (self._mission.IDSTRING, self.ID))
        super(Everest, self).plot_aperture(ax, labelsize=12)

        if show:
            pl.show()
            pl.close()
        else:
            return fig, ax

    def plot(self, show=True, plot_raw=True, plot_gp=True,
             plot_bad=True, plot_out=True, plot_cbv=True,
             simple=False):
        '''
        Plots the final de-trended light curve.

        :param bool show: Show the plot or return the `(fig, ax)` instance? \
               Default :py:obj:`True`
        :param bool plot_raw: Show the raw light curve? Default :py:obj:`True`
        :param bool plot_gp: Show the GP model prediction? \
               Default :py:obj:`True`
        :param bool plot_bad: Show and indicate the bad data points? \
               Default :py:obj:`True`
        :param bool plot_out: Show and indicate the outliers? \
               Default :py:obj:`True`
        :param bool plot_cbv: Plot the CBV-corrected light curve? \
               Default :py:obj:`True`. If :py:obj:`False`, plots the \
               de-trended but uncorrected light curve.

        '''

        log.info('Plotting the light curve...')

        # Set up axes
        if plot_raw:
            fig, axes = pl.subplots(2, figsize=(13, 9), sharex=True)
            fig.subplots_adjust(hspace=0.1)
            axes = [axes[1], axes[0]]
            if plot_cbv:
                fluxes = [self.fcor, self.fraw]
            else:
                fluxes = [self.flux, self.fraw]
            labels = ['EVEREST Flux', 'Raw Flux']
        else:
            fig, axes = pl.subplots(1, figsize=(13, 6))
            axes = [axes]
            if plot_cbv:
                fluxes = [self.fcor]
            else:
                fluxes = [self.flux]
            labels = ['EVEREST Flux']
        fig.canvas.set_window_title('EVEREST Light curve')

        # Set up some stuff
        time = self.time
        badmask = self.badmask
        nanmask = self.nanmask
        outmask = self.outmask
        transitmask = self.transitmask
        fraw_err = self.fraw_err
        breakpoints = self.breakpoints
        if self.cadence == 'sc':
            ms = 2
        else:
            ms = 4

        # Get the cdpps
        cdpps = [[self.get_cdpp(self.flux), self.get_cdpp_arr(self.flux)],
                 [self.get_cdpp(self.fraw), self.get_cdpp_arr(self.fraw)]]
        self.cdpp = cdpps[0][0]
        self.cdpp_arr = cdpps[0][1]

        for n, ax, flux, label, c in zip([0, 1], axes, fluxes, labels, cdpps):

            # Initialize CDPP
            cdpp = c[0]
            cdpp_arr = c[1]

            # Plot the good data points
            ax.plot(self.apply_mask(time), self.apply_mask(flux),
                    ls='none', marker='.', color='k', markersize=ms, alpha=0.5)

            # Plot the outliers
            bnmask = np.array(
                list(set(np.concatenate([badmask, nanmask]))), dtype=int)
            bmask = [i for i in self.badmask if i not in self.nanmask]

            def O1(x): return x[outmask]

            def O2(x): return x[bmask]

            def O3(x): return x[transitmask]
            if plot_out:
                ax.plot(O1(time), O1(flux), ls='none', color="#777777",
                        marker='.', markersize=ms, alpha=0.5)
            if plot_bad:
                ax.plot(O2(time), O2(flux), 'r.', markersize=ms, alpha=0.25)
            ax.plot(O3(time), O3(flux), 'b.', markersize=ms, alpha=0.25)

            # Plot the GP
            if n == 0 and plot_gp and self.cadence != 'sc':
                gp = GP(self.kernel, self.kernel_params)
                gp.compute(self.apply_mask(time), self.apply_mask(fraw_err))
                med = np.nanmedian(self.apply_mask(flux))
                y, _ = gp.predict(self.apply_mask(flux) - med, time)
                y += med
                ax.plot(self.apply_mask(time), self.apply_mask(
                    y), 'r-', lw=0.5, alpha=0.5)

            # Appearance
            if n == 0:
                ax.set_xlabel('Time (%s)' %
                              self._mission.TIMEUNITS, fontsize=18)
            ax.set_ylabel(label, fontsize=18)
            for brkpt in breakpoints[:-1]:
                ax.axvline(time[brkpt], color='r', ls='--', alpha=0.25)
            if len(cdpp_arr) == 2:
                ax.annotate('%.2f ppm' % cdpp_arr[0], xy=(0.02, 0.975),
                            xycoords='axes fraction',
                            ha='left', va='top', fontsize=12, color='r',
                            zorder=99)
                ax.annotate('%.2f ppm' % cdpp_arr[1], xy=(0.98, 0.975),
                            xycoords='axes fraction',
                            ha='right', va='top', fontsize=12,
                            color='r', zorder=99)
            elif len(cdpp_arr) < 6:
                for n in range(len(cdpp_arr)):
                    if n > 0:
                        x = (self.time[self.breakpoints[n - 1]] - self.time[0]
                             ) / (self.time[-1] - self.time[0]) + 0.02
                    else:
                        x = 0.02
                    ax.annotate('%.2f ppm' % cdpp_arr[n], xy=(x, 0.975),
                                xycoords='axes fraction',
                                ha='left', va='top', fontsize=10,
                                zorder=99, color='r')
            else:
                ax.annotate('%.2f ppm' % cdpp, xy=(0.02, 0.975),
                            xycoords='axes fraction',
                            ha='left', va='top', fontsize=12,
                            color='r', zorder=99)
            ax.margins(0.01, 0.1)

            # Get y lims that bound 99% of the flux
            f = np.concatenate([np.delete(f, bnmask) for f in fluxes])
            N = int(0.995 * len(f))
            hi, lo = f[np.argsort(f)][[N, -N]]
            pad = (hi - lo) * 0.1
            ylim = (lo - pad, hi + pad)
            ax.set_ylim(ylim)
            ax.get_yaxis().set_major_formatter(Formatter.Flux)

            # Indicate off-axis outliers
            for i in np.where(flux < ylim[0])[0]:
                if i in bmask:
                    color = "#ffcccc"
                    if not plot_bad:
                        continue
                elif i in outmask:
                    color = "#cccccc"
                    if not plot_out:
                        continue
                elif i in nanmask:
                    continue
                else:
                    color = "#ccccff"
                ax.annotate('', xy=(time[i], ylim[0]), xycoords='data',
                            xytext=(0, 15), textcoords='offset points',
                            arrowprops=dict(arrowstyle="-|>", color=color))
            for i in np.where(flux > ylim[1])[0]:
                if i in bmask:
                    color = "#ffcccc"
                    if not plot_bad:
                        continue
                elif i in outmask:
                    color = "#cccccc"
                    if not plot_out:
                        continue
                elif i in nanmask:
                    continue
                else:
                    color = "#ccccff"
                ax.annotate('', xy=(time[i], ylim[1]), xycoords='data',
                            xytext=(0, -15), textcoords='offset points',
                            arrowprops=dict(arrowstyle="-|>", color=color))

        # Show total CDPP improvement
        pl.figtext(0.5, 0.94, '%s %d' % (self._mission.IDSTRING, self.ID),
                   fontsize=18, ha='center', va='bottom')
        pl.figtext(0.5, 0.905,
                   r'$%.2f\ \mathrm{ppm} \rightarrow %.2f\ \mathrm{ppm}$' %
                   (self.cdppr, self.cdpp), fontsize=14,
                   ha='center', va='bottom')

        if show:
            pl.show()
            pl.close()
        else:
            if plot_raw:
                return fig, axes
            else:
                return fig, axes[0]

    def dvs(self):
        '''
        Shows the data validation summary (DVS) for the target.

        '''

        DVS(self.ID, season=self.season, mission=self.mission,
            model=self.model_name, clobber=self.clobber)

    def plot_pipeline(self, pipeline, *args, **kwargs):
        '''
        Plots the light curve for the target de-trended with a given pipeline.

        :param str pipeline: The name of the pipeline (lowercase). Options \
               are 'everest2', 'everest1', and other mission-specific \
               pipelines. For `K2`, the available pipelines are 'k2sff' \
               and 'k2sc'.

        Additional :py:obj:`args` and :py:obj:`kwargs` are passed directly to
        the :py:func:`pipelines.plot` function of the mission.

        '''

        if pipeline != 'everest2':
            return getattr(missions, self.mission).pipelines.plot(self.ID,
                                                                  pipeline,
                                                                  *args,
                                                                  **kwargs)

        else:

            # We're going to plot the everest 2 light curve like we plot
            # the other pipelines for easy comparison
            plot_raw = kwargs.get('plot_raw', False)
            plot_cbv = kwargs.get('plot_cbv', True)
            show = kwargs.get('show', True)

            if plot_raw:
                y = self.fraw
                ylabel = 'Raw Flux'
            elif plot_cbv:
                y = self.fcor
                ylabel = "EVEREST2 Flux"
            else:
                y = self.flux
                ylabel = "EVEREST2 Flux"

            # Remove nans
            bnmask = np.concatenate([self.nanmask, self.badmask])
            time = np.delete(self.time, bnmask)
            flux = np.delete(y, bnmask)

            # Plot it
            fig, ax = pl.subplots(1, figsize=(10, 4))
            fig.subplots_adjust(bottom=0.15)
            ax.plot(time, flux, "k.", markersize=3, alpha=0.5)

            # Axis limits
            N = int(0.995 * len(flux))
            hi, lo = flux[np.argsort(flux)][[N, -N]]
            pad = (hi - lo) * 0.1
            ylim = (lo - pad, hi + pad)
            ax.set_ylim(ylim)

            # Plot bad data points
            ax.plot(self.time[self.badmask], y[self.badmask],
                    "r.", markersize=3, alpha=0.2)

            # Show the CDPP
            ax.annotate('%.2f ppm' % self._mission.CDPP(flux),
                        xy=(0.98, 0.975), xycoords='axes fraction',
                        ha='right', va='top', fontsize=12, color='r',
                        zorder=99)

            # Appearance
            ax.margins(0, None)
            ax.set_xlabel("Time (%s)" % self._mission.TIMEUNITS, fontsize=16)
            ax.set_ylabel(ylabel, fontsize=16)
            fig.canvas.set_window_title("EVEREST2: EPIC %d" % (self.ID))

            if show:
                pl.show()
                pl.close()
            else:
                return fig, ax

    def get_pipeline(self, *args, **kwargs):
        '''
        Returns the `time` and `flux` arrays for the target obtained by a given
        pipeline.

        Options :py:obj:`args` and :py:obj:`kwargs` are passed directly to
        the :py:func:`pipelines.get` function of the mission.

        '''

        return getattr(missions, self.mission).pipelines.get(self.ID, *args,
                                                             **kwargs)

    def mask_planet(self, t0, period, dur=0.2):
        '''
        Mask all of the transits/eclipses of a given planet/EB. After calling
        this method, you must re-compute the model by calling
        :py:meth:`compute` in order for the mask to take effect.

        :param float t0: The time of first transit (same units as light curve)
        :param float period: The period of the planet in days
        :param foat dur: The transit duration in days. Default 0.2

        '''

        mask = []
        t0 += np.ceil((self.time[0] - dur - t0) / period) * period
        for t in np.arange(t0, self.time[-1] + dur, period):
            mask.extend(np.where(np.abs(self.time - t) < dur / 2.)[0])
        self.transitmask = np.array(
            list(set(np.concatenate([self.transitmask, mask]))))

    def _plot_weights(self, show=True):
        '''
        .. warning:: Untested!

        '''

        # Set up the axes
        fig = pl.figure(figsize=(12, 12))
        fig.subplots_adjust(top=0.95, bottom=0.025, left=0.1, right=0.92)
        fig.canvas.set_window_title(
            '%s %d' % (self._mission.IDSTRING, self.ID))
        ax = [pl.subplot2grid((80, 130), (20 * j, 25 * i), colspan=23,
                              rowspan=18)
              for j in range(len(self.breakpoints) * 2)
              for i in range(1 + 2 * (self.pld_order - 1))]
        cax = [pl.subplot2grid((80, 130),
                               (20 * j, 25 * (1 + 2 * (self.pld_order - 1))),
                               colspan=4, rowspan=18)
               for j in range(len(self.breakpoints) * 2)]
        ax = np.array(ax).reshape(2 * len(self.breakpoints), -1)
        cax = np.array(cax)

        # Check number of segments
        if len(self.breakpoints) > 3:
            log.error('Cannot currently plot weights for light ' +
                      'curves with more than 3 segments.')
            return

        # Loop over all PLD orders and over all chunks
        npix = len(self.fpix[1])
        ap = self.aperture.flatten()
        ncol = 1 + 2 * (len(self.weights[0]) - 1)
        raw_weights = np.zeros(
            (len(self.breakpoints), ncol, self.aperture.shape[0],
             self.aperture.shape[1]), dtype=float)
        scaled_weights = np.zeros(
            (len(self.breakpoints), ncol, self.aperture.shape[0],
             self.aperture.shape[1]), dtype=float)

        # Loop over orders
        for o in range(len(self.weights[0])):
            if o == 0:
                oi = 0
            else:
                oi = 1 + 2 * (o - 1)

            # Loop over chunks
            for b in range(len(self.weights)):

                c = self.get_chunk(b)
                rw_ii = np.zeros(npix)
                rw_ij = np.zeros(npix)
                sw_ii = np.zeros(npix)
                sw_ij = np.zeros(npix)
                X = np.nanmedian(self.X(o, c), axis=0)

                # Compute all sets of pixels at this PLD order, then
                # loop over them and assign the weights to the correct pixels
                sets = np.array(list(multichoose(np.arange(npix).T, o + 1)))
                for i, s in enumerate(sets):
                    if (o == 0) or (s[0] == s[1]):
                        # Not the cross-terms
                        j = s[0]
                        rw_ii[j] += self.weights[b][o][i]
                        sw_ii[j] += X[i] * self.weights[b][o][i]
                    else:
                        # Cross-terms
                        for j in s:
                            rw_ij[j] += self.weights[b][o][i]
                            sw_ij[j] += X[i] * self.weights[b][o][i]

                # Make the array 2D and plot it
                rw = np.zeros_like(ap, dtype=float)
                sw = np.zeros_like(ap, dtype=float)
                n = 0
                for i, a in enumerate(ap):
                    if (a & 1):
                        rw[i] = rw_ii[n]
                        sw[i] = sw_ii[n]
                        n += 1
                raw_weights[b][oi] = rw.reshape(*self.aperture.shape)
                scaled_weights[b][oi] = sw.reshape(*self.aperture.shape)

                if o > 0:
                    # Make the array 2D and plot it
                    rw = np.zeros_like(ap, dtype=float)
                    sw = np.zeros_like(ap, dtype=float)
                    n = 0
                    for i, a in enumerate(ap):
                        if (a & 1):
                            rw[i] = rw_ij[n]
                            sw[i] = sw_ij[n]
                            n += 1
                    raw_weights[b][oi + 1] = rw.reshape(*self.aperture.shape)
                    scaled_weights[b][oi +
                                      1] = sw.reshape(*self.aperture.shape)

        # Plot the images
        log.info('Plotting the PLD weights...')
        rdbu = pl.get_cmap('RdBu_r')
        rdbu.set_bad('k')
        for b in range(len(self.weights)):
            rmax = max([-raw_weights[b][o].min() for o in range(ncol)] +
                       [raw_weights[b][o].max() for o in range(ncol)])
            smax = max([-scaled_weights[b][o].min() for o in range(ncol)] +
                       [scaled_weights[b][o].max() for o in range(ncol)])
            for o in range(ncol):
                imr = ax[2 * b, o].imshow(raw_weights[b][o], aspect='auto',
                                          interpolation='nearest', cmap=rdbu,
                                          origin='lower', vmin=-rmax,
                                          vmax=rmax)
                ims = ax[2 * b + 1, o].imshow(scaled_weights[b][o],
                                              aspect='auto',
                                              interpolation='nearest',
                                              cmap=rdbu, origin='lower',
                                              vmin=-smax, vmax=smax)

            # Colorbars
            def fmt(x, pos):
                a, b = '{:.0e}'.format(x).split('e')
                b = int(b)
                if float(a) > 0:
                    a = r'+' + a
                elif float(a) == 0:
                    return ''
                return r'${} \times 10^{{{}}}$'.format(a, b)
            cbr = pl.colorbar(imr, cax=cax[2 * b], format=FuncFormatter(fmt))
            cbr.ax.tick_params(labelsize=8)
            cbs = pl.colorbar(
                ims, cax=cax[2 * b + 1], format=FuncFormatter(fmt))
            cbs.ax.tick_params(labelsize=8)

        # Plot aperture contours
        def PadWithZeros(vector, pad_width, iaxis, kwargs):
            vector[:pad_width[0]] = 0
            vector[-pad_width[1]:] = 0
            return vector
        ny, nx = self.aperture.shape
        contour = np.zeros((ny, nx))
        contour[np.where(self.aperture)] = 1
        contour = np.lib.pad(contour, 1, PadWithZeros)
        highres = zoom(contour, 100, order=0, mode='nearest')
        extent = np.array([-1, nx, -1, ny])
        for axis in ax.flatten():
            axis.contour(highres, levels=[
                         0.5], extent=extent, origin='lower', colors='r',
                         linewidths=1)

            # Check for saturated columns
            for x in range(self.aperture.shape[0]):
                for y in range(self.aperture.shape[1]):
                    if self.aperture[x][y] == AP_SATURATED_PIXEL:
                        axis.fill([y - 0.5, y + 0.5, y + 0.5, y - 0.5],
                                  [x - 0.5, x - 0.5, x + 0.5, x + 0.5],
                                  fill=False, hatch='xxxxx', color='r', lw=0)

            axis.set_xlim(-0.5, nx - 0.5)
            axis.set_ylim(-0.5, ny - 0.5)
            axis.set_xticks([])
            axis.set_yticks([])

        # Labels
        titles = [r'$1^{\mathrm{st}}$',
                  r'$2^{\mathrm{nd}}\ (i = j)$',
                  r'$2^{\mathrm{nd}}\ (i \neq j)$',
                  r'$3^{\mathrm{rd}}\ (i = j)$',
                  r'$3^{\mathrm{rd}}\ (i \neq j)$'] + ['' for i in range(10)]
        for i, axis in enumerate(ax[0]):
            axis.set_title(titles[i], fontsize=12)
        for j in range(len(self.weights)):
            ax[2 * j, 0].text(-0.55, -0.15, r'$%d$' % (j + 1),
                              fontsize=16, transform=ax[2 * j, 0].transAxes)
            ax[2 * j, 0].set_ylabel(r'$w_{ij}$', fontsize=18)
            ax[2 * j + 1,
                0].set_ylabel(r'$\bar{X}_{ij} \cdot w_{ij}$', fontsize=18)

        if show:
            pl.show()
            pl.close()
        else:
            return fig, ax, cax

    def _plot_chunks(self, show=True, plot_bad=True, plot_out=True):
        '''

        '''

        log.info('Plotting the light curve...')

        # Set up axes
        fig, ax = pl.subplots(len(self.breakpoints), figsize=(10, 8))
        fig.canvas.set_window_title('EVEREST Light curve')
        if self.cadence == 'sc':
            ms = 2
        else:
            ms = 4

        # Calculate the fluxes and cdpps
        fluxes = [None for i in range(len(self.breakpoints))]
        cdpps = [None for i in range(len(self.breakpoints))]
        for b in range(len(self.breakpoints)):
            m = self.get_masked_chunk(b)
            c = np.arange(len(self.time))
            mK = GetCovariance(self.kernel, self.kernel_params,
                               self.time[m], self.fraw_err[m])
            med = np.nanmedian(self.fraw[m])
            f = self.fraw[m] - med
            A = np.zeros((len(m), len(m)))
            B = np.zeros((len(c), len(m)))
            for n in range(self.pld_order):
                if (self.lam_idx >= n) and (self.lam[b][n] is not None):
                    XM = self.X(n, m)
                    XC = self.X(n, c)
                    A += self.lam[b][n] * np.dot(XM, XM.T)
                    B += self.lam[b][n] * np.dot(XC, XM.T)
                    del XM, XC
            W = np.linalg.solve(mK + A, f)
            model = np.dot(B, W)
            del A, B, W
            fluxes[b] = self.fraw - model + np.nanmedian(model)
            cdpps[b] = self.get_cdpp_arr(fluxes[b])

        # Loop over all chunks
        for i in range(len(self.breakpoints)):

            # Get current flux/cdpp
            flux = fluxes[i]
            cdpp_arr = cdpps[i]

            # Plot the good data points
            ax[i].plot(self.apply_mask(self.time), self.apply_mask(
                flux), ls='none', marker='.', color='k', markersize=ms,
                alpha=0.5)

            # Plot the outliers
            bnmask = np.array(
                list(set(np.concatenate([self.badmask, self.nanmask]))),
                dtype=int)

            def O1(x): return x[self.outmask]

            def O2(x): return x[bnmask]

            def O3(x): return x[self.transitmask]
            if plot_out:
                ax[i].plot(O1(self.time), O1(flux), ls='none',
                           color="#777777", marker='.', markersize=ms,
                           alpha=0.5)
            if plot_bad:
                ax[i].plot(O2(self.time), O2(flux), 'r.',
                           markersize=ms, alpha=0.25)
            ax[i].plot(O3(self.time), O3(flux), 'b.',
                       markersize=ms, alpha=0.25)

            # Appearance
            if i == len(self.breakpoints) - 1:
                ax[i].set_xlabel('Time (%s)' %
                                 self._mission.TIMEUNITS, fontsize=18)
            ax[i].set_ylabel('Flux %d' % (i + 1), fontsize=18)
            for brkpt in self.breakpoints[:-1]:
                ax[i].axvline(self.time[brkpt], color='r', ls='--', alpha=0.25)
            if len(self.breakpoints) == 2:
                ax[i].annotate('%.2f ppm' % cdpp_arr[0], xy=(0.02, 0.975),
                               xycoords='axes fraction',
                               ha='left', va='top', fontsize=12, color='r',
                               zorder=99)
                ax[i].annotate('%.2f ppm' % cdpp_arr[1], xy=(0.98, 0.975),
                               xycoords='axes fraction',
                               ha='right', va='top', fontsize=12,
                               color='r', zorder=99)
            elif len(self.breakpoints) < 6:
                for n in range(len(self.breakpoints)):
                    if n > 0:
                        x = (self.time[self.breakpoints[n - 1]] - self.time[0]
                             ) / (self.time[-1] - self.time[0]) + 0.02
                    else:
                        x = 0.02
                    ax[i].annotate('%.2f ppm' % cdpp_arr[n], xy=(x, 0.975),
                                   xycoords='axes fraction',
                                   ha='left', va='top', fontsize=10,
                                   zorder=99, color='r')
            else:
                ax[i].annotate('%.2f ppm' % cdpp_arr[0], xy=(0.02, 0.975),
                               xycoords='axes fraction',
                               ha='left', va='top', fontsize=12,
                               color='r', zorder=99)
            ax[i].margins(0.01, 0.1)
            if i == 0:
                a = self.time[0]
            else:
                a = self.time[self.breakpoints[i - 1]]
            b = self.time[self.breakpoints[i]]
            ax[i].axvspan(a, b, color='b', alpha=0.1, zorder=-99)

            # Get y lims that bound 99% of the flux
            f = np.concatenate([np.delete(f, bnmask) for f in fluxes])
            N = int(0.995 * len(f))
            hi, lo = f[np.argsort(f)][[N, -N]]
            pad = (hi - lo) * 0.1
            ylim = (lo - pad, hi + pad)
            ax[i].set_ylim(ylim)
            ax[i].get_yaxis().set_major_formatter(Formatter.Flux)

            # Indicate off-axis outliers
            for j in np.where(flux < ylim[0])[0]:
                if j in bnmask:
                    color = "#ffcccc"
                    if not plot_bad:
                        continue
                elif j in self.outmask:
                    color = "#cccccc"
                    if not plot_out:
                        continue
                else:
                    color = "#ccccff"
                ax[i].annotate('', xy=(self.time[j], ylim[0]), xycoords='data',
                               xytext=(0, 15), textcoords='offset points',
                               arrowprops=dict(arrowstyle="-|>", color=color))
            for j in np.where(flux > ylim[1])[0]:
                if j in bnmask:
                    color = "#ffcccc"
                    if not plot_bad:
                        continue
                elif j in self.outmask:
                    color = "#cccccc"
                    if not plot_out:
                        continue
                else:
                    color = "#ccccff"
                ax[i].annotate('', xy=(self.time[j], ylim[1]), xycoords='data',
                               xytext=(0, -15), textcoords='offset points',
                               arrowprops=dict(arrowstyle="-|>", color=color))

        if show:
            pl.show()
            pl.close()
        else:
            return fig, axes

    def _save_npz(self):
        '''
        Saves all of the de-trending information to disk in an `npz` file

        '''

        # Save the data
        d = dict(self.__dict__)
        d.pop('_weights', None)
        d.pop('_A', None)
        d.pop('_B', None)
        d.pop('_f', None)
        d.pop('_mK', None)
        d.pop('K', None)
        d.pop('dvs', None)
        d.pop('clobber', None)
        d.pop('clobber_tpf', None)
        d.pop('_mission', None)
        d.pop('debug', None)
        np.savez(os.path.join(self.dir, self.name + '.npz'), **d)

    def optimize(self, piter=3, pmaxf=300, ppert=0.1):
        '''
        Runs :py:obj:`pPLD` on the target in an attempt to further optimize the
        values of the PLD priors. See :py:class:`everest.detrender.pPLD`.

        '''

        self._save_npz()
        optimized = pPLD(self.ID, piter=piter, pmaxf=pmaxf,
                         ppert=ppert, debug=True, clobber=True)
        optimized.publish()
        self.reset()

    def plot_folded(self, t0, period, dur=0.2):
        '''
        Plot the light curve folded on a given `period` and centered at `t0`.
        When plotting folded transits, please mask them using
        :py:meth:`mask_planet` and re-compute the model using
        :py:meth:`compute`.

        :param float t0: The time at which to center the plot \
               (same units as light curve)
        :param float period: The period of the folding operation
        :param float dur: The transit duration in days. Default 0.2

        '''

        # Mask the planet
        self.mask_planet(t0, period, dur)

        # Whiten
        gp = GP(self.kernel, self.kernel_params, white=False)
        gp.compute(self.apply_mask(self.time), self.apply_mask(self.fraw_err))
        med = np.nanmedian(self.apply_mask(self.flux))
        y, _ = gp.predict(self.apply_mask(self.flux) - med, self.time)
        fwhite = (self.flux - y)
        fwhite /= np.nanmedian(fwhite)

        # Fold
        tfold = (self.time - t0 - period / 2.) % period - period / 2.

        # Crop
        inds = np.where(np.abs(tfold) < 2 * dur)[0]
        x = tfold[inds]
        y = fwhite[inds]

        # Plot
        fig, ax = pl.subplots(1, figsize=(9, 5))
        fig.subplots_adjust(bottom=0.125)
        ax.plot(x, y, 'k.', alpha=0.5)

        # Get ylims
        yfin = np.delete(y, np.where(np.isnan(y)))
        lo, hi = yfin[np.argsort(yfin)][[3, -3]]
        pad = (hi - lo) * 0.1
        ylim = (lo - pad, hi + pad)
        ax.set_ylim(*ylim)

        # Appearance
        ax.set_xlabel(r'Time (days)', fontsize=18)
        ax.set_ylabel(r'Normalized Flux', fontsize=18)
        fig.canvas.set_window_title(
            '%s %d' % (self._mission.IDSTRING, self.ID))

        pl.show()

    def plot_transit_model(self, show=True, fold=None, ax=None):
        '''
        Plot the light curve de-trended with a join instrumental + transit
        model with the best fit transit model overlaid. The transit model
        should be specified using the :py:obj:`transit_model` attribute
        and should be an instance or list of instances of
        :py:class:`everest.transit.TransitModel`.

        :param bool show: Show the plot, or return the `fig, ax` instances? \
               Default `True`
        :param str fold: The name of the planet/transit model on which to \
               fold. If only one model is present, can be set to \
               :py:obj:`True`. Default :py:obj:`False` \
               (does not fold the data).
        :param ax: A `matplotlib` axis instance to use for plotting. \
               Default :py:obj:`None`

        '''

        if self.transit_model is None:
            raise ValueError("No transit model provided!")
        if self.transit_depth is None:
            self.compute()
        if fold is not None:
            if (fold is True and len(self.transit_model) > 1) or \
               (type(fold) is not str):
                raise Exception(
                    "Kwarg `fold` should be the name of the transit " +
                    "model on which to fold the data.")
            if fold is True:
                # We are folding on the first index of `self.transit_model`
                fold = 0
            elif type(fold) is str:
                # Figure out the index of the transit model on which to fold
                fold = np.argmax(
                    [fold == tm.name for tm in self.transit_model])
            log.info('Plotting the transit model folded ' +
                     'on transit model index %d...' % fold)
        else:
            log.info('Plotting the transit model...')

        # Set up axes
        if ax is None:
            if fold is not None:
                fig, ax = pl.subplots(1, figsize=(8, 5))
            else:
                fig, ax = pl.subplots(1, figsize=(13, 6))
            fig.canvas.set_window_title('EVEREST Light curve')
        else:
            fig = pl.gcf()

        # Set up some stuff
        if self.cadence == 'sc':
            ms = 2
        else:
            ms = 4

        # Fold?
        if fold is not None:
            times = self.transit_model[fold].params.get('times', None)
            if times is not None:
                time = self.time - \
                    [times[np.argmin(np.abs(ti - times))] for ti in self.time]
                t0 = times[0]
            else:
                t0 = self.transit_model[fold].params.get('t0', 0.)
                period = self.transit_model[fold].params.get('per', 10.)
                time = (self.time - t0 - period / 2.) % period - period / 2.
            dur = 0.01 * \
                len(np.where(self.transit_model[fold](
                    np.linspace(t0 - 0.5, t0 + 0.5, 100)) < 0)[0])
        else:
            time = self.time
            ax.plot(self.apply_mask(time), self.apply_mask(self.flux),
                    ls='none', marker='.', color='k', markersize=ms, alpha=0.5)
            ax.plot(time[self.outmask], self.flux[self.outmask],
                    ls='none', marker='.', color='k', markersize=ms, alpha=0.5)
            ax.plot(time[self.transitmask], self.flux[self.transitmask],
                    ls='none', marker='.', color='k', markersize=ms, alpha=0.5)

        # Plot the transit + GP model
        med = np.nanmedian(self.apply_mask(self.flux))
        transit_model = \
            med * np.sum([depth * tm(self.time)
                          for tm, depth in zip(self.transit_model,
                                               self.transit_depth)], axis=0)
        gp = GP(self.kernel, self.kernel_params, white=False)
        gp.compute(self.apply_mask(self.time), self.apply_mask(self.fraw_err))
        y, _ = gp.predict(self.apply_mask(
            self.flux - transit_model) - med, self.time)
        if fold is not None:
            flux = (self.flux - y) / med
            ax.plot(self.apply_mask(time), self.apply_mask(flux),
                    ls='none', marker='.', color='k', markersize=ms, alpha=0.5)
            ax.plot(time[self.outmask], flux[self.outmask], ls='none',
                    marker='.', color='k', markersize=ms, alpha=0.5)
            ax.plot(time[self.transitmask], flux[self.transitmask],
                    ls='none', marker='.', color='k', markersize=ms, alpha=0.5)
            hires_time = np.linspace(-5 * dur, 5 * dur, 1000)
            hires_transit_model = 1 + \
                self.transit_depth[fold] * \
                self.transit_model[fold](hires_time + t0)
            ax.plot(hires_time, hires_transit_model, 'r-', lw=1, alpha=1)
        else:
            flux = self.flux
            y += med
            y += transit_model
            ax.plot(time, y, 'r-', lw=1, alpha=1)

        # Plot the bad data points
        bnmask = np.array(
            list(set(np.concatenate([self.badmask, self.nanmask]))), dtype=int)
        bmask = [i for i in self.badmask if i not in self.nanmask]
        ax.plot(time[bmask], flux[bmask], 'r.', markersize=ms, alpha=0.25)

        # Appearance
        ax.set_ylabel('EVEREST Flux', fontsize=18)
        ax.margins(0.01, 0.1)
        if fold is not None:
            ax.set_xlabel('Time From Transit Center (days)', fontsize=18)
            ax.set_xlim(-3 * dur, 3 * dur)
        else:
            ax.set_xlabel('Time (%s)' % self._mission.TIMEUNITS, fontsize=18)
            for brkpt in self.breakpoints[:-1]:
                ax.axvline(time[brkpt], color='r', ls='--', alpha=0.25)
            ax.get_yaxis().set_major_formatter(Formatter.Flux)

        # Get y lims that bound most of the flux
        if fold is not None:
            lo = np.min(hires_transit_model)
            pad = 1.5 * (1 - lo)
            ylim = (lo - pad, 1 + pad)
        else:
            f = np.delete(flux, bnmask)
            N = int(0.995 * len(f))
            hi, lo = f[np.argsort(f)][[N, -N]]
            pad = (hi - lo) * 0.1
            ylim = (lo - pad, hi + pad)
        ax.set_ylim(ylim)

        # Indicate off-axis outliers
        for i in np.where(flux < ylim[0])[0]:
            if i in bmask:
                color = "#ffcccc"
            else:
                color = "#ccccff"
            ax.annotate('', xy=(time[i], ylim[0]), xycoords='data',
                        xytext=(0, 15), textcoords='offset points',
                        arrowprops=dict(arrowstyle="-|>", color=color,
                        alpha=0.5))
        for i in np.where(flux > ylim[1])[0]:
            if i in bmask:
                color = "#ffcccc"
            else:
                color = "#ccccff"
            ax.annotate('', xy=(time[i], ylim[1]), xycoords='data',
                        xytext=(0, -15), textcoords='offset points',
                        arrowprops=dict(arrowstyle="-|>", color=color,
                        alpha=0.5))

        if show:
            pl.show()
            pl.close()
        else:
            return fig, ax
