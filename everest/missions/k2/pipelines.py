#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`pipelines.py` - Pipeline comparison
--------------------------------------------

Routines for retrieving, processing and plotting the light curves
generated by different pipelines.

'''

from __future__ import division, print_function, absolute_import, \
     unicode_literals
from ...config import EVEREST_SRC, EVEREST_DAT, EVEREST_DEV
from ...mathutils import SavGol
import os
import sys
import shutil
import k2plr
from k2plr.config import KPLR_ROOT
from six.moves import urllib
import matplotlib.pyplot as pl
import numpy as np
import warnings
try:
    import pyfits
except ImportError:
    try:
        import astropy.io.fits as pyfits
    except ImportError:
        raise Exception('Please install the `pyfits` package.')
import logging
log = logging.getLogger(__name__)

#: The supported pipelines
Pipelines = ['everest2', 'everest1', 'k2sff', 'k2sc', 'raw']


def get(ID, pipeline='everest2', campaign=None):
    '''
    Returns the `time` and `flux` for a given EPIC `ID` and
    a given `pipeline` name.

    '''

    log.info('Downloading %s light curve for %d...' % (pipeline, ID))

    # Dev version hack
    if EVEREST_DEV:
        if pipeline.lower() == 'everest2' or pipeline.lower() == 'k2sff':
            from . import Season, TargetDirectory, FITSFile
            if campaign is None:
                campaign = Season(ID)
            fits = os.path.join(TargetDirectory(
                ID, campaign), FITSFile(ID, campaign))
            newdir = os.path.join(KPLR_ROOT, "data", "everest", str(ID))
            if not os.path.exists(newdir):
                os.makedirs(newdir)
            if os.path.exists(fits):
                shutil.copy(fits, newdir)

    if pipeline.lower() == 'everest2':
        s = k2plr.EVEREST(ID, version=2, sci_campaign=campaign)
        time = s.time
        flux = s.flux
    elif pipeline.lower() == 'everest1':
        s = k2plr.EVEREST(ID, version=1, sci_campaign=campaign)
        time = s.time
        flux = s.flux
    elif pipeline.lower() == 'k2sff':
        s = k2plr.K2SFF(ID, sci_campaign=campaign)
        time = s.time
        flux = s.fcor
        # Normalize to the median flux
        s = k2plr.EVEREST(ID, version=2, sci_campaign=campaign)
        flux *= np.nanmedian(s.flux)
    elif pipeline.lower() == 'k2sc':
        s = k2plr.K2SC(ID, sci_campaign=campaign)
        time = s.time
        flux = s.pdcflux
    elif pipeline.lower() == 'raw':
        s = k2plr.EVEREST(ID, version=2, raw=True, sci_campaign=campaign)
        time = s.time
        flux = s.flux
    else:
        raise ValueError('Invalid pipeline: `%s`.' % pipeline)

    return time, flux


def plot(ID, pipeline='everest2', show=True, campaign=None):
    '''
    Plots the de-trended flux for the given EPIC `ID` and for
    the specified `pipeline`.

    '''

    # Get the data
    time, flux = get(ID, pipeline=pipeline, campaign=campaign)

    # Remove nans
    mask = np.where(np.isnan(flux))[0]
    time = np.delete(time, mask)
    flux = np.delete(flux, mask)

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

    # Show the CDPP
    from .k2 import CDPP
    ax.annotate('%.2f ppm' % CDPP(flux),
                xy=(0.98, 0.975), xycoords='axes fraction',
                ha='right', va='top', fontsize=12, color='r', zorder=99)

    # Appearance
    ax.margins(0, None)
    ax.set_xlabel("Time (BJD - 2454833)", fontsize=16)
    ax.set_ylabel("%s Flux" % pipeline.upper(), fontsize=16)
    fig.canvas.set_window_title("%s: EPIC %d" % (pipeline.upper(), ID))

    if show:
        pl.show()
        pl.close()
    else:
        return fig, ax


def get_cdpp(campaign, pipeline='everest2'):
    '''
    Computes the CDPP for a given `campaign` and a given `pipeline`.
    Stores the results in a file under "/missions/k2/tables/".

    '''

    # Imports
    from .k2 import CDPP
    from .utils import GetK2Campaign

    # Check pipeline
    assert pipeline.lower() in Pipelines, 'Invalid pipeline: `%s`.' % pipeline

    # Create file if it doesn't exist
    file = os.path.join(EVEREST_SRC, 'missions', 'k2',
                        'tables', 'c%02d_%s.cdpp' % (int(campaign), pipeline))
    if not os.path.exists(file):
        open(file, 'a').close()

    # Get all EPIC stars
    stars = GetK2Campaign(campaign, epics_only=True)
    nstars = len(stars)

    # Remove ones we've done
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        done = np.loadtxt(file, dtype=float)
    if len(done):
        done = [int(s) for s in done[:, 0]]
    stars = list(set(stars) - set(done))
    n = len(done) + 1

    # Open the output file
    with open(file, 'a', 1) as outfile:

        # Loop over all to get the CDPP
        for EPIC in stars:

            # Progress
            sys.stdout.write('\rRunning target %d/%d...' % (n, nstars))
            sys.stdout.flush()
            n += 1

            # Get the CDPP
            try:
                _, flux = get(EPIC, pipeline=pipeline, campaign=campaign)
                mask = np.where(np.isnan(flux))[0]
                flux = np.delete(flux, mask)
                cdpp = CDPP(flux)
            except (urllib.error.HTTPError, urllib.error.URLError,
                    TypeError, ValueError, IndexError):
                print("{:>09d} {:>15.3f}".format(EPIC, 0), file=outfile)
                continue

            # Log to file
            print("{:>09d} {:>15.3f}".format(EPIC, cdpp), file=outfile)


def get_outliers(campaign, pipeline='everest2', sigma=5):
    '''
    Computes the number of outliers for a given `campaign`
    and a given `pipeline`.
    Stores the results in a file under "/missions/k2/tables/".

    :param int sigma: The sigma level at which to clip outliers. Default 5

    '''

    # Imports
    from .utils import GetK2Campaign
    client = k2plr.API()

    # Check pipeline
    assert pipeline.lower() in Pipelines, 'Invalid pipeline: `%s`.' % pipeline

    # Create file if it doesn't exist
    file = os.path.join(EVEREST_SRC, 'missions', 'k2',
                        'tables', 'c%02d_%s.out' % (int(campaign), pipeline))
    if not os.path.exists(file):
        open(file, 'a').close()

    # Get all EPIC stars
    stars = GetK2Campaign(campaign, epics_only=True)
    nstars = len(stars)

    # Remove ones we've done
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        done = np.loadtxt(file, dtype=float)
    if len(done):
        done = [int(s) for s in done[:, 0]]
    stars = list(set(stars) - set(done))
    n = len(done) + 1

    # Open the output file
    with open(file, 'a', 1) as outfile:

        # Loop over all to get the CDPP
        for EPIC in stars:

            # Progress
            sys.stdout.write('\rRunning target %d/%d...' % (n, nstars))
            sys.stdout.flush()
            n += 1

            # Get the number of outliers
            try:
                time, flux = get(EPIC, pipeline=pipeline, campaign=campaign)

                # Get the raw K2 data
                tpf = os.path.join(KPLR_ROOT, "data", "k2",
                                   "target_pixel_files",
                                   "%09d" % EPIC,
                                   "ktwo%09d-c%02d_lpd-targ.fits.gz"
                                   % (EPIC, campaign))
                if not os.path.exists(tpf):
                    client.k2_star(EPIC).get_target_pixel_files(fetch=True)
                with pyfits.open(tpf) as f:
                    k2_qual = np.array(f[1].data.field('QUALITY'), dtype=int)
                    k2_time = np.array(
                        f[1].data.field('TIME'), dtype='float64')
                    mask = []
                    for b in [1, 2, 3, 4, 5, 6, 7, 8, 9,
                              11, 12, 13, 14, 16, 17]:
                        mask += list(np.where(k2_qual & 2 ** (b - 1))[0])
                    mask = np.array(sorted(list(set(mask))))

                # Fill in missing cadences, if any
                tol = 0.005
                if not ((len(time) == len(k2_time)) and (np.abs(time[0]
                        - k2_time[0]) < tol) and (np.abs(time[-1]
                                                  - k2_time[-1]) < tol)):
                    ftmp = np.zeros_like(k2_time) * np.nan
                    j = 0
                    for i, t in enumerate(k2_time):
                        if np.abs(time[j] - t) < tol:
                            ftmp[i] = flux[j]
                            j += 1
                            if j == len(time) - 1:
                                break
                    flux = ftmp

                # Remove flagged cadences
                flux = np.delete(flux, mask)

                # Remove nans
                nanmask = np.where(np.isnan(flux))[0]
                flux = np.delete(flux, nanmask)

                # Iterative sigma clipping
                inds = np.array([], dtype=int)
                m = 1
                while len(inds) < m:
                    m = len(inds)
                    f = SavGol(np.delete(flux, inds))
                    med = np.nanmedian(f)
                    MAD = 1.4826 * np.nanmedian(np.abs(f - med))
                    inds = np.append(inds, np.where(
                        (f > med + sigma * MAD) | (f < med - sigma * MAD))[0])
                nout = len(inds)
                ntot = len(flux)

            except (urllib.error.HTTPError, urllib.error.URLError,
                    TypeError, ValueError, IndexError):
                print("{:>09d} {:>5d} {:>5d}".format(
                    EPIC, -1, -1), file=outfile)
                continue

            # Log to file
            print("{:>09d} {:>5d} {:>5d}".format(
                EPIC, nout, ntot), file=outfile)
