#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`fits.py` - FITS conversion
-----------------------------------

Converts the :py:mod:`everest` de-trended light curves to the
*FITS* format using
`PyFits <http://stsdas.stsci.edu/download/docs/The_PyFITS_Handbook.pdf>`_.
These *FITS* files make up the public :py:mod:`everest` catalog.

'''

from __future__ import division, print_function, absolute_import, \
     unicode_literals
from . import __version__ as EVEREST_VERSION
from .config import EVEREST_DAT, EVEREST_SRC, QUALITY_BAD, QUALITY_NAN, \
     QUALITY_OUT, QUALITY_REC, QUALITY_TRN, EVEREST_MAJOR_MINOR
try:
    import pyfits
except ImportError:
    try:
        import astropy.io.fits as pyfits
    except ImportError:
        raise Exception('Please install the `pyfits` package.')
import os
import numpy as np
from time import strftime
import logging
log = logging.getLogger(__name__)

__all__ = ['MakeFITS']


def PrimaryHDU(model):
    '''
    Construct the primary HDU file containing basic header info.

    '''

    # Get mission cards
    cards = model._mission.HDUCards(model.meta, hdu=0)
    if 'KEPMAG' not in [c[0] for c in cards]:
        cards.append(('KEPMAG', model.mag, 'Kepler magnitude'))

    # Add EVEREST info
    cards.append(('COMMENT', '************************'))
    cards.append(('COMMENT', '*     EVEREST INFO     *'))
    cards.append(('COMMENT', '************************'))
    cards.append(('MISSION', model.mission, 'Mission name'))
    cards.append(('VERSION', EVEREST_MAJOR_MINOR, 'EVEREST pipeline version'))
    cards.append(('SUBVER', EVEREST_VERSION, 'EVEREST pipeline subversion'))
    cards.append(('DATE', strftime('%Y-%m-%d'),
                  'EVEREST file creation date (YYYY-MM-DD)'))

    # Create the HDU
    header = pyfits.Header(cards=cards)
    hdu = pyfits.PrimaryHDU(header=header)

    return hdu


def LightcurveHDU(model):
    '''
    Construct the data HDU file containing the arrays and the observing info.

    '''

    # Get mission cards
    cards = model._mission.HDUCards(model.meta, hdu=1)

    # Add EVEREST info
    cards.append(('COMMENT', '************************'))
    cards.append(('COMMENT', '*     EVEREST INFO     *'))
    cards.append(('COMMENT', '************************'))
    cards.append(('MISSION', model.mission, 'Mission name'))
    cards.append(('VERSION', EVEREST_MAJOR_MINOR, 'EVEREST pipeline version'))
    cards.append(('SUBVER', EVEREST_VERSION, 'EVEREST pipeline subversion'))
    cards.append(('DATE', strftime('%Y-%m-%d'),
                  'EVEREST file creation date (YYYY-MM-DD)'))
    cards.append(('MODEL', model.name, 'Name of EVEREST model used'))
    cards.append(('APNAME', model.aperture_name, 'Name of aperture used'))
    cards.append(('BPAD', model.bpad, 'Chunk overlap in cadences'))
    for c in range(len(model.breakpoints)):
        cards.append(
            ('BRKPT%02d' % (c + 1), model.breakpoints[c],
             'Light curve breakpoint'))
    cards.append(('CBVNUM', model.cbv_num, 'Number of CBV signals to recover'))
    cards.append(('CBVNITER', model.cbv_niter,
                  'Number of CBV SysRem iterations'))
    cards.append(('CBVWIN', model.cbv_win, 'Window size for smoothing CBVs'))
    cards.append(('CBVORD', model.cbv_order, 'Order when smoothing CBVs'))
    cards.append(('CDIVS', model.cdivs, 'Cross-validation subdivisions'))
    cards.append(('CDPP', model.cdpp, 'Average de-trended CDPP'))
    cards.append(('CDPPR', model.cdppr, 'Raw CDPP'))
    cards.append(('CDPPV', model.cdppv, 'Average validation CDPP'))
    cards.append(('CDPPG', model.cdppg, 'Average GP-de-trended CDPP'))
    for i in range(99):
        try:
            cards.append(('CDPP%02d' % (i + 1),
                          model.cdpp_arr[i] if not np.isnan(
                model.cdpp_arr[i]) else 0, 'Chunk de-trended CDPP'))
            cards.append(('CDPPR%02d' % (
                i + 1), model.cdppr_arr[i] if not np.isnan(
                model.cdppr_arr[i]) else 0, 'Chunk raw CDPP'))
            cards.append(('CDPPV%02d' % (i + 1),
                          model.cdppv_arr[i] if not np.isnan(
                model.cdppv_arr[i]) else 0, 'Chunk validation CDPP'))
        except:
            break
    cards.append(
        ('CVMIN', model.cv_min, 'Cross-validation objective function'))
    cards.append(
        ('GITER', model.giter, 'Number of GP optimiziation iterations'))
    cards.append(
        ('GMAXF', model.giter, 'Max number of GP function evaluations'))
    cards.append(('GPFACTOR', model.gp_factor,
                  'GP amplitude initialization factor'))
    cards.append(('KERNEL', model.kernel, 'GP kernel name'))
    if model.kernel == 'Basic':
        cards.append(
            ('GPWHITE', model.kernel_params[0],
             'GP white noise amplitude (e-/s)'))
        cards.append(
            ('GPRED', model.kernel_params[1],
             'GP red noise amplitude (e-/s)'))
        cards.append(
            ('GPTAU', model.kernel_params[2],
             'GP red noise timescale (days)'))
    elif model.kernel == 'QuasiPeriodic':
        cards.append(
            ('GPWHITE', model.kernel_params[0],
             'GP white noise amplitude (e-/s)'))
        cards.append(
            ('GPRED', model.kernel_params[1], 'GP red noise amplitude (e-/s)'))
        cards.append(('GPGAMMA', model.kernel_params[2], 'GP scale factor'))
        cards.append(('GPPER', model.kernel_params[3], 'GP period (days)'))
    for c in range(len(model.breakpoints)):
        for o in range(model.pld_order):
            cards.append(('LAMB%02d%02d' % (c + 1, o + 1),
                          model.lam[c][o], 'Cross-validation parameter'))
            if model.name == 'iPLD':
                cards.append(('RECL%02d%02d' % (c + 1, o + 1),
                              model.reclam[c][o],
                              'Cross-validation parameter'))
    cards.append(('LEPS', model.leps, 'Cross-validation tolerance'))
    cards.append(('MAXPIX', model.max_pixels, 'Maximum size of TPF aperture'))
    for i, source in enumerate(model.nearby[:99]):
        cards.append(('NRBY%02dID' %
                      (i + 1), source['ID'], 'Nearby source ID'))
        cards.append(
            ('NRBY%02dX' % (i + 1), source['x'], 'Nearby source X position'))
        cards.append(
            ('NRBY%02dY' % (i + 1), source['y'], 'Nearby source Y position'))
        cards.append(
            ('NRBY%02dM' % (i + 1), source['mag'], 'Nearby source magnitude'))
        cards.append(('NRBY%02dX0' %
                      (i + 1), source['x0'], 'Nearby source reference X'))
        cards.append(('NRBY%02dY0' %
                      (i + 1), source['y0'], 'Nearby source reference Y'))
    for i, n in enumerate(model.neighbors):
        cards.append(
            ('NEIGH%02d' % i, model.neighbors[i],
             'Neighboring star used to de-trend'))
    cards.append(('OITER', model.oiter, 'Number of outlier search iterations'))
    cards.append(('OPTGP', model.optimize_gp, 'GP optimization performed?'))
    cards.append(
        ('OSIGMA', model.osigma, 'Outlier tolerance (standard deviations)'))
    for i, planet in enumerate(model.planets):
        cards.append(
            ('P%02dT0' % (i + 1), planet[0], 'Planet transit time (days)'))
        cards.append(
            ('P%02dPER' % (i + 1), planet[1], 'Planet transit period (days)'))
        cards.append(
            ('P%02dDUR' % (i + 1), planet[2],
             'Planet transit duration (days)'))
    cards.append(('PLDORDER', model.pld_order, 'PLD de-trending order'))
    cards.append(('SATUR', model.saturated, 'Is target saturated?'))
    cards.append(('SATTOL', model.saturation_tolerance,
                  'Fractional saturation tolerance'))

    # Add the EVEREST quality flags to the QUALITY array
    quality = np.array(model.quality)
    quality[np.array(model.badmask, dtype=int)] += 2 ** (QUALITY_BAD - 1)
    quality[np.array(model.nanmask, dtype=int)] += 2 ** (QUALITY_NAN - 1)
    quality[np.array(model.outmask, dtype=int)] += 2 ** (QUALITY_OUT - 1)
    quality[np.array(model.recmask, dtype=int)] += 2 ** (QUALITY_REC - 1)
    quality[np.array(model.transitmask, dtype=int)] += 2 ** (QUALITY_TRN - 1)

    # When de-trending, we interpolated to fill in NaN fluxes. Here
    # we insert the NaNs back in, since there's no actual physical
    # information at those cadences.
    flux = np.array(model.flux)
    flux[model.nanmask] = np.nan

    # Create the arrays list
    arrays = [pyfits.Column(name='CADN', format='D', array=model.cadn),
              pyfits.Column(name='FLUX', format='D', array=flux, unit='e-/s'),
              pyfits.Column(name='FRAW', format='D',
                            array=model.fraw, unit='e-/s'),
              pyfits.Column(name='FRAW_ERR', format='D',
                            array=model.fraw_err, unit='e-/s'),
              pyfits.Column(name='QUALITY', format='J', array=quality),
              pyfits.Column(name='TIME', format='D',
                            array=model.time, unit='BJD - 2454833')]

    # Add the CBVs
    if model.fcor is not None:
        arrays += [pyfits.Column(name='FCOR', format='D',
                                 array=model.fcor, unit='e-/s')]
        for n in range(model.XCBV.shape[1]):
            arrays += [pyfits.Column(name='CBV%02d' %
                                     (n + 1), format='D',
                                     array=model.XCBV[:, n])]

    # Did we subtract a background term?
    if hasattr(model.bkg, '__len__'):
        arrays.append(pyfits.Column(name='BKG', format='D',
                                    array=model.bkg, unit='e-/s'))

    # Create the HDU
    header = pyfits.Header(cards=cards)
    cols = pyfits.ColDefs(arrays)
    hdu = pyfits.BinTableHDU.from_columns(cols, header=header, name='ARRAYS')

    return hdu


def PixelsHDU(model):
    '''
    Construct the HDU containing the pixel-level light curve.

    '''

    # Get mission cards
    cards = model._mission.HDUCards(model.meta, hdu=2)

    # Add EVEREST info
    cards = []
    cards.append(('COMMENT', '************************'))
    cards.append(('COMMENT', '*     EVEREST INFO     *'))
    cards.append(('COMMENT', '************************'))
    cards.append(('MISSION', model.mission, 'Mission name'))
    cards.append(('VERSION', EVEREST_MAJOR_MINOR, 'EVEREST pipeline version'))
    cards.append(('SUBVER', EVEREST_VERSION, 'EVEREST pipeline subversion'))
    cards.append(('DATE', strftime('%Y-%m-%d'),
                  'EVEREST file creation date (YYYY-MM-DD)'))

    # Create the HDU
    header = pyfits.Header(cards=cards)

    # The pixel timeseries
    arrays = [pyfits.Column(name='FPIX', format='%dD' %
                            model.fpix.shape[1], array=model.fpix)]

    # The first order PLD vectors for all the neighbors (npixels, ncadences)
    X1N = model.X1N
    if X1N is not None:
        arrays.append(pyfits.Column(name='X1N', format='%dD' %
                                    X1N.shape[1], array=X1N))

    cols = pyfits.ColDefs(arrays)
    hdu = pyfits.BinTableHDU.from_columns(cols, header=header, name='PIXELS')

    return hdu


def ApertureHDU(model):
    '''
    Construct the HDU containing the aperture used to de-trend.

    '''

    # Get mission cards
    cards = model._mission.HDUCards(model.meta, hdu=3)

    # Add EVEREST info
    cards.append(('COMMENT', '************************'))
    cards.append(('COMMENT', '*     EVEREST INFO     *'))
    cards.append(('COMMENT', '************************'))
    cards.append(('MISSION', model.mission, 'Mission name'))
    cards.append(('VERSION', EVEREST_MAJOR_MINOR, 'EVEREST pipeline version'))
    cards.append(('SUBVER', EVEREST_VERSION, 'EVEREST pipeline subversion'))
    cards.append(('DATE', strftime('%Y-%m-%d'),
                  'EVEREST file creation date (YYYY-MM-DD)'))

    # Create the HDU
    header = pyfits.Header(cards=cards)
    hdu = pyfits.ImageHDU(data=model.aperture,
                          header=header, name='APERTURE MASK')

    return hdu


def ImagesHDU(model):
    '''
    Construct the HDU containing sample postage stamp images of the target.

    '''

    # Get mission cards
    cards = model._mission.HDUCards(model.meta, hdu=4)

    # Add EVEREST info
    cards.append(('COMMENT', '************************'))
    cards.append(('COMMENT', '*     EVEREST INFO     *'))
    cards.append(('COMMENT', '************************'))
    cards.append(('MISSION', model.mission, 'Mission name'))
    cards.append(('VERSION', EVEREST_MAJOR_MINOR, 'EVEREST pipeline version'))
    cards.append(('SUBVER', EVEREST_VERSION, 'EVEREST pipeline subversion'))
    cards.append(('DATE', strftime('%Y-%m-%d'),
                  'EVEREST file creation date (YYYY-MM-DD)'))

    # The images
    format = '%dD' % model.pixel_images[0].shape[1]
    arrays = [pyfits.Column(name='STAMP1', format=format,
                            array=model.pixel_images[0]),
              pyfits.Column(name='STAMP2', format=format,
                            array=model.pixel_images[1]),
              pyfits.Column(name='STAMP3', format=format,
                            array=model.pixel_images[2])]

    # Create the HDU
    header = pyfits.Header(cards=cards)
    cols = pyfits.ColDefs(arrays)
    hdu = pyfits.BinTableHDU.from_columns(
        cols, header=header, name='POSTAGE STAMPS')

    return hdu


def HiResHDU(model):
    '''
    Construct the HDU containing the hi res image of the target.

    '''

    # Get mission cards
    cards = model._mission.HDUCards(model.meta, hdu=5)

    # Add EVEREST info
    cards.append(('COMMENT', '************************'))
    cards.append(('COMMENT', '*     EVEREST INFO     *'))
    cards.append(('COMMENT', '************************'))
    cards.append(('MISSION', model.mission, 'Mission name'))
    cards.append(('VERSION', EVEREST_MAJOR_MINOR, 'EVEREST pipeline version'))
    cards.append(('SUBVER', EVEREST_VERSION, 'EVEREST pipeline subversion'))
    cards.append(('DATE', strftime('%Y-%m-%d'),
                  'EVEREST file creation date (YYYY-MM-DD)'))

    # Create the HDU
    header = pyfits.Header(cards=cards)
    if model.hires is not None:
        hdu = pyfits.ImageHDU(
            data=model.hires, header=header, name='HI RES IMAGE')
    else:
        hdu = pyfits.ImageHDU(data=np.empty(
            (0, 0), dtype=float), header=header, name='HI RES IMAGE')
    return hdu


def MakeFITS(model, fitsfile=None):
    '''
    Generate a FITS file for a given :py:mod:`everest` run.

    :param model: An :py:mod:`everest` model instance

    '''

    # Get the fits file name
    if fitsfile is None:
        outfile = os.path.join(model.dir, model._mission.FITSFile(
            model.ID, model.season, model.cadence))
    else:
        outfile = os.path.join(model.dir, fitsfile)
    if os.path.exists(outfile) and not model.clobber:
        return
    elif os.path.exists(outfile):
        os.remove(outfile)

    log.info('Generating FITS file...')

    # Create the HDUs
    primary = PrimaryHDU(model)
    lightcurve = LightcurveHDU(model)
    pixels = PixelsHDU(model)
    aperture = ApertureHDU(model)
    images = ImagesHDU(model)
    hires = HiResHDU(model)

    # Combine to get the HDUList
    hdulist = pyfits.HDUList(
        [primary, lightcurve, pixels, aperture, images, hires])

    # Output to the FITS file
    hdulist.writeto(outfile)

    return
