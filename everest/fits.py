#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`fits.py` - FITS conversion
-----------------------------------

Converts the :py:mod:`everest` de-trended light curves to the
*FITS* format using `PyFits <http://stsdas.stsci.edu/download/docs/The_PyFITS_Handbook.pdf>`_.
These *FITS* files make up the public :py:mod:`everest` catalog.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from . import __version__ as EVEREST_VERSION
from .config import EVEREST_DAT, EVEREST_SRC, QUALITY_BAD, QUALITY_NAN, QUALITY_OUT
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
  cards = model._mission.HDUCards(model.meta, hdu = 0)
  
  # Add EVEREST info
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*     EVEREST INFO     *'))
  cards.append(('COMMENT', '************************'))
  cards.append(('MISSION', model.mission, 'Mission name'))
  cards.append(('VERSION', EVEREST_VERSION, 'EVEREST pipeline version'))
  cards.append(('DATE', strftime('%Y-%m-%d'), 'EVEREST file creation date (YYYY-MM-DD)'))
  
  # Create the HDU
  header = pyfits.Header(cards = cards)
  hdu = pyfits.PrimaryHDU(header = header)

  return hdu

def LongCadenceHDU(model):
  '''
  Construct the data HDU file containing the arrays and the observing info.
  
  '''
  
  # Get mission cards
  cards = model._mission.HDUCards(model.meta, hdu = 1)
      
  # Add EVEREST info
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*     EVEREST INFO     *'))
  cards.append(('COMMENT', '************************'))
  cards.append(('MISSION', model.mission, 'Mission name'))
  cards.append(('VERSION', EVEREST_VERSION, 'EVEREST pipeline version'))
  cards.append(('DATE', strftime('%Y-%m-%d'), 'EVEREST file creation date (YYYY-MM-DD)'))
  cards.append(('MODEL', model.name, 'Name of EVEREST model used'))
  cards.append(('APNAME', model.aperture_name, 'Name of aperture used'))
  cards.append(('BPAD', model.bpad, 'Chunk overlap in cadences'))
  for c in range(len(model.breakpoints)):
    cards.append(('BRKPT%d' % (c + 1), model.breakpoints[c], 'Light curve breakpoint'))
  cards.append(('CDIVS', model.cdivs, 'Cross-validation subdivisions'))
  cards.append(('CDPP6', model.cdpp6, 'Average de-trended 6-hr CDPP'))
  cards.append(('CDPPR', model.cdppr, 'Raw 6-hr CDPP'))
  cards.append(('CDPPV', model.cdppv, 'Average validation 6-hr CDPP'))
  cards.append(('CDPPG', model.gppp, 'Average GP-de-trended 6-hr CDPP'))
  for i in range(99):
    try:
      cards.append(('CDPP6%02d' % (i + 1), model.cdpp6_arr[i], 'Chunk de-trended 6-hr CDPP'))
      cards.append(('CDPPR%02d' % (i + 1), model.cdppr_arr[i], 'Chunk raw 6-hr CDPP'))
      cards.append(('CDPPV%02d' % (i + 1), model.cdppv_arr[i], 'Chunk validation 6-hr CDPP'))
    except:
      break
  cards.append(('GITER', model.giter, 'Number of GP optimiziation iterations'))
  cards.append(('GPFACTOR', model.gp_factor, 'GP amplitude initialization factor'))
  cards.append(('GPWHITE', model.kernel_params[0], 'GP white noise amplitude (e-/s)'))
  cards.append(('GPRED', model.kernel_params[1], 'GP red noise amplitude (e-/s)'))
  cards.append(('GPTAU', model.kernel_params[2], 'GP red noise timescale (days)'))
  for c in range(len(model.breakpoints)):
    for o in range(model.pld_order):
      cards.append(('LAMBDA%d%d' % (c + 1, o + 1), model.lam[c][o], 'Cross-validation parameter'))
  cards.append(('LEPS', model.leps, 'Cross-validation tolerance'))
  cards.append(('MAXPIX', model.max_pixels, 'Maximum size of TPF aperture'))
  for i, source in enumerate(model.nearby[:99]):
    cards.append(('NRBY%02dID' % (i + 1), source['ID'], 'Nearby source ID'))
    cards.append(('NRBY%02dX' % (i + 1), source['x'], 'Nearby source X position'))
    cards.append(('NRBY%02dY' % (i + 1), source['y'], 'Nearby source Y position'))
    cards.append(('NRBY%02dM' % (i + 1), source['mag'], 'Nearby source magnitude'))
    cards.append(('NRBY%02dX0' % (i + 1), source['x0'], 'Nearby source reference X'))
    cards.append(('NRBY%02dY0' % (i + 1), source['y0'], 'Nearby source reference Y'))
  for i, n in enumerate(model.neighbors):
    cards.append(('NEIGH%02d' % i, model.neighbors[i], 'Neighboring star used to de-trend'))
  cards.append(('OITER', model.oiter, 'Number of outlier search iterations'))
  cards.append(('OPTGP', model.optimize_gp, 'GP optimization performed?'))
  cards.append(('OSIGMA', model.osigma, 'Outlier tolerance (standard deviations)'))
  cards.append(('PLDORDER', model.pld_order, 'PLD de-trending order'))
  cards.append(('RECRSVE', model.recursive, 'Recursive PLD?'))
  cards.append(('SATUR', model.saturated, 'Is target saturated?'))
  cards.append(('SATTOL', model.saturation_tolerance, 'Fractional saturation tolerance'))
  
  # Add the EVEREST quality flags to the QUALITY array
  quality = np.array(model.quality)
  quality[model.badmask] += 2 ** (QUALITY_BAD - 1)
  quality[model.nanmask] += 2 ** (QUALITY_NAN - 1)
  quality[model.outmask] += 2 ** (QUALITY_OUT - 1)

  # Create the arrays list
  arrays = [pyfits.Column(name = 'CADN', format = 'D', array = model.cadn),
            pyfits.Column(name = 'FLUX', format = 'D', array = model.flux, unit = 'e-/s'),
            pyfits.Column(name = 'FRAW', format = 'D', array = model.fraw, unit = 'e-/s'),
            pyfits.Column(name = 'FRAW_ERR', format = 'D', array = model.fraw_err, unit = 'e-/s'),
            pyfits.Column(name = 'QUALITY', format = 'J', array = quality),
            pyfits.Column(name = 'TIME', format = 'D', array = model.time, unit = 'BJD - 2454833')]
  
  # Did we subtract a background term?
  if hasattr(model.bkg, '__len__'):
    arrays.append(pyfits.Column(name = 'BKG', format = 'D', array = model.bkg, unit = 'BJD - 2454833'))     
  
  # Create the HDU
  header = pyfits.Header(cards = cards)
  cols = pyfits.ColDefs(arrays)
  hdu = pyfits.BinTableHDU.from_columns(cols, header = header, name = 'LC ARRAYS')

  return hdu

def ShortCadenceHDU(model):
  '''
  Construct the data HDU file containing the arrays and the observing info.
  
  '''
  
  # Get mission cards
  cards = model._mission.HDUCards(model.meta, hdu = 6)
   
  # Add EVEREST info
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*     EVEREST INFO     *'))
  cards.append(('COMMENT', '************************'))
  cards.append(('MISSION', model.mission, 'Mission name'))
  cards.append(('VERSION', EVEREST_VERSION, 'EVEREST pipeline version'))
  cards.append(('DATE', strftime('%Y-%m-%d'), 'EVEREST file creation date (YYYY-MM-DD)'))
  cards.append(('MODEL', model.name, 'Name of EVEREST model used'))
  cards.append(('SC_BPAD', model.sc_bpad, 'Chunk overlap in cadences'))

  # Arrays
  sc_quality = np.array(model.sc_quality)
  sc_quality[model.sc_badmask] += 2 ** (QUALITY_BAD - 1)
  sc_quality[model.sc_nanmask] += 2 ** (QUALITY_NAN - 1)
  sc_quality[model.sc_outmask] += 2 ** (QUALITY_OUT - 1)
  arrays = [pyfits.Column(name = 'SC_CADN', format = 'D', array = model.sc_cadn),
            pyfits.Column(name = 'SC_FLUX', format = 'D', array = model.sc_flux, unit = 'e-/s'),
            pyfits.Column(name = 'SC_FRAW', format = 'D', array = model.sc_fraw, unit = 'e-/s'),
            pyfits.Column(name = 'SC_FERR', format = 'D', array = model.sc_fraw_err, unit = 'e-/s'),
            pyfits.Column(name = 'SC_QUAL', format = 'J', array = sc_quality),
            pyfits.Column(name = 'SC_TIME', format = 'D', array = model.sc_time, unit = 'BJD - 2454833')]
  
  # Did we subtract a background term?
  if hasattr(model.sc_bkg, '__len__'):
    arrays.append(pyfits.Column(name = 'SC_BKG', format = 'D', array = model.sc_bkg, unit = 'BJD - 2454833'))
  
  # Create the HDU
  header = pyfits.Header(cards = cards)
  cols = pyfits.ColDefs(arrays)
  hdu = pyfits.BinTableHDU.from_columns(cols, header = header, name = 'SC ARRAYS')

  return hdu

def LongCadencePixelHDU(model):
  '''
  Construct the HDU containing the pixel-level light curve.
  
  '''
  
  # Get mission cards
  cards = model._mission.HDUCards(model.meta, hdu = 2)
  
  # Add EVEREST info
  cards = []
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*     EVEREST INFO     *'))
  cards.append(('COMMENT', '************************'))
  cards.append(('MISSION', model.mission, 'Mission name'))
  cards.append(('VERSION', EVEREST_VERSION, 'EVEREST pipeline version'))
  cards.append(('DATE', strftime('%Y-%m-%d'), 'EVEREST file creation date (YYYY-MM-DD)'))
  
  # Create the HDU
  header = pyfits.Header(cards = cards)
  
  # The pixel timeseries
  arrays = [pyfits.Column(name = 'FPIX', format = '%dD' % model.fpix.shape[1], array = model.fpix)]
  
  # The first order PLD vectors for all the neighbors (npixels, ncadences)
  X1N = model.X1N
  if X1N is not None:
    arrays.append(pyfits.Column(name = 'X1N', format = '%dD' % X1N.shape[1], array = X1N))
    
  cols = pyfits.ColDefs(arrays)
  hdu = pyfits.BinTableHDU.from_columns(cols, header = header, name = 'LC PIXELS')

  return hdu

def ShortCadencePixelHDU(model):
  '''
  Construct the HDU containing the pixel-level light curve.
  
  '''

  # Get mission cards
  cards = model._mission.HDUCards(model.meta, hdu = 7)

  # Add EVEREST info
  cards = []
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*     EVEREST INFO     *'))
  cards.append(('COMMENT', '************************'))
  cards.append(('MISSION', model.mission, 'Mission name'))
  cards.append(('VERSION', EVEREST_VERSION, 'EVEREST pipeline version'))
  cards.append(('DATE', strftime('%Y-%m-%d'), 'EVEREST file creation date (YYYY-MM-DD)'))
  
  # Create the HDU
  header = pyfits.Header(cards = cards)
  
  # The pixel timeseries
  arrays = [pyfits.Column(name = 'SC_FPIX', format = '%dD' % model.sc_fpix.shape[1], array = model.sc_fpix)]
  
  # The first order PLD vectors for all the neighbors (npixels, ncadences)
  if model.sc_X1N is not None:
    arrays.append(pyfits.Column(name = 'SC_X1N', format = '%dD' % model.sc_X1N.shape[1], array = model.sc_X1N))
    
  cols = pyfits.ColDefs(arrays)
  hdu = pyfits.BinTableHDU.from_columns(cols, header = header, name = 'SC PIXELS')

  return hdu

def ApertureHDU(model):
  '''
  Construct the HDU containing the aperture used to de-trend.
  
  '''
  
  # Get mission cards
  cards = model._mission.HDUCards(model.meta, hdu = 3)
  
  # Add EVEREST info
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*     EVEREST INFO     *'))
  cards.append(('COMMENT', '************************'))
  cards.append(('MISSION', model.mission, 'Mission name'))
  cards.append(('VERSION', EVEREST_VERSION, 'EVEREST pipeline version'))
  cards.append(('DATE', strftime('%Y-%m-%d'), 'EVEREST file creation date (YYYY-MM-DD)'))
  
  # Create the HDU
  header = pyfits.Header(cards = cards)
  hdu = pyfits.ImageHDU(data = model.aperture, header = header, name = 'APERTURE MASK')

  return hdu

def ImagesHDU(model):
  '''
  Construct the HDU containing sample postage stamp images of the target.
  
  '''
  
  # Get mission cards
  cards = model._mission.HDUCards(model.meta, hdu = 4)
  
  # Add EVEREST info
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*     EVEREST INFO     *'))
  cards.append(('COMMENT', '************************'))
  cards.append(('MISSION', model.mission, 'Mission name'))
  cards.append(('VERSION', EVEREST_VERSION, 'EVEREST pipeline version'))
  cards.append(('DATE', strftime('%Y-%m-%d'), 'EVEREST file creation date (YYYY-MM-DD)'))
  
  # The images
  format = '%dD' % model.pixel_images[0].shape[1]
  arrays = [pyfits.Column(name = 'STAMP1', format = format, array = model.pixel_images[0]),
            pyfits.Column(name = 'STAMP2', format = format, array = model.pixel_images[1]),
            pyfits.Column(name = 'STAMP3', format = format, array = model.pixel_images[2])]
  
  # Create the HDU
  header = pyfits.Header(cards = cards)
  cols = pyfits.ColDefs(arrays)
  hdu = pyfits.BinTableHDU.from_columns(cols, header = header, name = 'POSTAGE STAMPS')

  return hdu

def HiResHDU(model):
  '''
  Construct the HDU containing the hi res image of the target.
  
  '''
  
  # Get mission cards
  cards = model._mission.HDUCards(model.meta, hdu = 5)
  
  # Add EVEREST info
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*     EVEREST INFO     *'))
  cards.append(('COMMENT', '************************'))
  cards.append(('MISSION', model.mission, 'Mission name'))
  cards.append(('VERSION', EVEREST_VERSION, 'EVEREST pipeline version'))
  cards.append(('DATE', strftime('%Y-%m-%d'), 'EVEREST file creation date (YYYY-MM-DD)'))
  
  # Create the HDU
  header = pyfits.Header(cards = cards)
  hdu = pyfits.ImageHDU(data = model.hires, header = header, name = 'HI RES IMAGE')

  return hdu

def MakeFITS(model):
  '''
  Generate a FITS file for a given :py:mod:`everest` run.
  
  :param model: An :py:mod:`everest` model instance
  
  '''
  
  log.info('Generating FITS file...')
  
  # Get the fits file name
  outfile = os.path.join(model.dir, model._mission.FITSFile(model.ID, model.season))
  if os.path.exists(outfile) and not model.clobber:
    return
  elif os.path.exists(outfile):
    os.remove(outfile)
      
  # Create the HDUs
  primary = PrimaryHDU(model)
  long_cadence = LongCadenceHDU(model)
  lc_pixels = LongCadencePixelHDU(model)
  aperture = ApertureHDU(model)
  images = ImagesHDU(model)
  hires = HiResHDU(model)
  
  # Combine to get the HDUList
  if model.has_sc:
    short_cadence = ShortCadenceHDU(model)
    sc_pixels = ShortCadencePixelHDU(model)
    hdulist = pyfits.HDUList([primary, long_cadence, lc_pixels, aperture, images, hires, short_cadence, sc_pixels])
  else:
    hdulist = pyfits.HDUList([primary, long_cadence, lc_pixels, aperture, images, hires])
  
  # Output to the FITS file
  hdulist.writeto(outfile)
  
  return