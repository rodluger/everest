#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`fits.py` - *FITS* conversion
-------------------------------------

Converts the :py:mod:`everest` de-trended light curves to the
*FITS* format using `PyFits <http://stsdas.stsci.edu/download/docs/The_PyFITS_Handbook.pdf>`_.
These *FITS* files make up the public :py:mod:`everest` catalog.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from ... import __version__ as EVEREST_VERSION
from ...config import EVEREST_DAT, EVEREST_SRC, QUALITY_BAD, QUALITY_NAN, QUALITY_OUT
import k2plr as kplr
from k2plr.config import KPLR_ROOT
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
  
  # Get info from the TPF Primary HDU Header
  tpf_header = model.meta[0]
  entries = ['TELESCOP', 'INSTRUME', 'OBJECT', 'KEPLERID', 'CHANNEL', 'MODULE', 
             'OUTPUT', 'CAMPAIGN', 'DATA_REL', 'OBSMODE', 'MISSION', 'TTABLEID',
             'RADESYS', 'RA_OBJ', 'DEC_OBJ',  'EQUINOX', 'KEPMAG']  
  cards = []
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*      KEPLER INFO     *'))
  cards.append(('COMMENT', '************************'))
  for entry in entries:
    cards.append(tuple(tpf_header[entry]))
  
  # Add EVEREST info
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*     EVEREST INFO     *'))
  cards.append(('COMMENT', '************************'))
  cards.append(('VERSION', EVEREST_VERSION, 'EVEREST pipeline version'))
  cards.append(('DATE', strftime('%Y-%m-%d'), 'EVEREST file creation date (YYYY-MM-DD)'))
  
  # Create the HDU
  header = pyfits.Header(cards = cards)
  hdu = pyfits.PrimaryHDU(header = header)

  return hdu

def LongCadenceHDU(model):
  '''
  Construct the data HDU file containing the arrays and Kepler observing info.
  
  '''
  
  # Get info from the TPF BinTable HDU Header
  tpf_header = model.meta[1]
  entries = ['WCSN4P', 'WCAX4P', '1CTY4P', '2CTY4P', '1CUN4P', '2CUN4P', '1CRV4P', 
             '2CRV4P', '1CDL4P', '2CDL4P', '1CRP4P', '2CRP4P', 'WCAX4', '1CTYP4', 
             '2CTYP4', '1CRPX4', '2CRPX4', '1CRVL4', '2CRVL4', '1CUNI4', '2CUNI4', 
             '1CDLT4', '2CDLT4', '11PC4', '12PC4', '21PC4', '22PC4', 'WCSN5P', 
             'WCAX5P', '1CTY5P', '2CTY5P', '1CUN5P', '2CUN5P', '1CRV5P', '2CRV5P', 
             '1CDL5P', '2CDL5P', '1CRP5P', '2CRP5P', 'WCAX5', '1CTYP5', '2CTYP5', 
             '1CRPX5', '2CRPX5', '1CRVL5', '2CRVL5', '1CUNI5', '2CUNI5', '1CDLT5', 
             '2CDLT5', '11PC5', '12PC5', '21PC5', '22PC5', 'WCSN6P', 'WCAX6P', 
             '1CTY6P', '2CTY6P', '1CUN6P', '2CUN6P', '1CRV6P', '2CRV6P', '1CDL6P', 
             '2CDL6P', '1CRP6P', '2CRP6P', 'WCAX6', '1CTYP6', '2CTYP6', '1CRPX6', 
             '2CRPX6', '1CRVL6', '2CRVL6', '1CUNI6', '2CUNI6', '1CDLT6', '2CDLT6', 
             '11PC6', '12PC6', '21PC6', '22PC6', 'WCSN7P', 'WCAX7P', '1CTY7P', 
             '2CTY7P', '1CUN7P', '2CUN7P', '1CRV7P', '2CRV7P', '1CDL7P', '2CDL7P', 
             '1CRP7P', '2CRP7P', 'WCAX7', '1CTYP7', '2CTYP7', '1CRPX7', '2CRPX7', 
             '1CRVL7', '2CRVL7', '1CUNI7', '2CUNI7', '1CDLT7', '2CDLT7', '11PC7', 
             '12PC7', '21PC7', '22PC7', 'WCSN8P', 'WCAX8P', '1CTY8P', '2CTY8P', 
             '1CUN8P', '2CUN8P', '1CRV8P', '2CRV8P', '1CDL8P', '2CDL8P', '1CRP8P', 
             '2CRP8P', 'WCAX8', '1CTYP8', '2CTYP8', '1CRPX8', '2CRPX8', '1CRVL8', 
             '2CRVL8', '1CUNI8', '2CUNI8', '1CDLT8', '2CDLT8', '11PC8', '12PC8', 
             '21PC8', '22PC8', 'WCSN9P', 'WCAX9P', '1CTY9P', '2CTY9P', '1CUN9P', 
             '2CUN9P', '1CRV9P', '2CRV9P', '1CDL9P', '2CDL9P', '1CRP9P', '2CRP9P', 
             'WCAX9', '1CTYP9', '2CTYP9', '1CRPX9', '2CRPX9', '1CRVL9', '2CRVL9', 
             '1CUNI9', '2CUNI9', '1CDLT9', '2CDLT9', '11PC9', '12PC9', '21PC9', 
             '22PC9', 'INHERIT', 'EXTNAME', 'EXTVER', 'TELESCOP', 'INSTRUME', 
             'OBJECT', 'KEPLERID', 'RADESYS', 'RA_OBJ', 'DEC_OBJ', 'EQUINOX', 
             'EXPOSURE', 'TIMEREF', 'TASSIGN', 'TIMESYS', 'BJDREFI', 'BJDREFF', 
             'TIMEUNIT', 'TELAPSE', 'LIVETIME', 'TSTART', 'TSTOP', 'LC_START', 
             'LC_END', 'DEADC', 'TIMEPIXR', 'TIERRELA', 'INT_TIME', 'READTIME', 
             'FRAMETIM', 'NUM_FRM', 'TIMEDEL', 'DATE-OBS', 'DATE-END', 'BACKAPP', 
             'DEADAPP', 'VIGNAPP', 'GAIN', 'READNOIS', 'NREADOUT', 'TIMSLICE', 
             'MEANBLCK', 'LCFXDOFF', 'SCFXDOFF']
  cards = []
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*      KEPLER INFO     *'))
  cards.append(('COMMENT', '************************'))
  for entry in entries:
    cards.append(tuple(tpf_header[entry]))
      
  # Add EVEREST info
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*     EVEREST INFO     *'))
  cards.append(('COMMENT', '************************'))
  cards.append(('VERSION', EVEREST_VERSION, 'EVEREST pipeline version'))
  cards.append(('DATE', strftime('%Y-%m-%d'), 'EVEREST file creation date (YYYY-MM-DD)'))
  cards.append(('MODEL', model.name, 'Name of EVEREST model used'))
  cards.append(('APNAME', model.aperture_name, 'Name of aperture used'))
  cards.append(('BPAD', model.bpad, 'Chunk overlap in cadences'))
  cards.append(('BRKPT', model.breakpoints[0], 'Light curve breakpoint'))
  cards.append(('CDIVS', model.cdivs, 'Cross-validation subdivisions'))
  cards.append(('CDPP6', model.cdpp6, 'Average de-trended 6-hr CDPP'))
  cards.append(('CDPP6A', model.cdpp6_arr[0], 'First chunk de-trended 6-hr CDPP'))
  try:
    cards.append(('CDPP6B', model.cdpp6_arr[1], 'Second chunk de-trended 6-hr CDPP'))
  except IndexError:
    cards.append(('CDPP6B', np.nan, 'Second chunk de-trended 6-hr CDPP'))
  cards.append(('CDPP6R', model.cdppr, 'Raw 6-hr CDPP'))
  cards.append(('CDPP6RA', model.cdppr_arr[0], 'First chunk raw 6-hr CDPP'))
  try:
    cards.append(('CDPP6RB', model.cdppr_arr[1], 'Second chunk raw 6-hr CDPP'))
  except IndexError:
    cards.append(('CDPP6RB', np.nan, 'Second chunk raw 6-hr CDPP'))
  cards.append(('CDPP6V', model.cdpp6, 'Average validation 6-hr CDPP'))
  cards.append(('CDPP6VA', model.cdpp6_arr[0], 'First chunk validation 6-hr CDPP'))
  try:
    cards.append(('CDPP6VB', model.cdpp6_arr[1], 'Second chunk validation 6-hr CDPP'))
  except IndexError:
    cards.append(('CDPP6VB', np.nan, 'Second chunk validation 6-hr CDPP'))
  cards.append(('CDPP6G', model.gppp, 'Average GP-de-trended 6-hr CDPP'))
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
  for i, n in enumerate(model.neighbors):
    cards.append(('NEIGH%02d' % i, model.neighbors[i], 'Neighboring star used to de-trend'))
  cards.append(('OITER', model.oiter, 'Number of outlier search iterations'))
  cards.append(('OPTGP', model.optimize_gp, 'GP optimization performed?'))
  cards.append(('OSIGMA', model.osigma, 'Outlier tolerance (standard deviations)'))
  cards.append(('PLDORDER', model.pld_order, 'PLD de-trending order'))
  cards.append(('RECRSVE', model.recursive, 'Recursive PLD?'))
  cards.append(('SATUR', model.saturated, 'Is target saturated?'))
  cards.append(('SATAP', model.saturated_aperture_name, 'Saturated aperture name'))
  cards.append(('SATTOL', model.saturation_tolerance, 'Fractional saturation tolerance'))
  
  # Add the EVEREST quality flags to the QUALITY array
  quality = np.array(model.quality)
  quality[model.badmask] += QUALITY_BAD
  quality[model.nanmask] += QUALITY_NAN
  quality[model.outmask] += QUALITY_OUT
  
  # Create the arrays list
  arrays = [pyfits.Column(name = 'CADN', format = 'D', array = model.cadn),
            pyfits.Column(name = 'FLUX', format = 'D', array = model.flux, unit = 'e-/s'),
            pyfits.Column(name = 'FRAW', format = 'D', array = model.fraw, unit = 'e-/s'),
            pyfits.Column(name = 'FRAW_ERR', format = 'D', array = model.fraw_err, unit = 'e-/s'),
            pyfits.Column(name = 'QUALITY', format = 'J', array = quality),
            pyfits.Column(name = 'TIME', format = 'D', array = model.time, unit = 'BJD - 2454833')]
  
  # Did we subtract a background term?
  if hasattr(model.bkg, '__len__'):
    arrays.append(pyfits.Column(name = 'BKG', format = 'D', array = time, unit = 'BJD - 2454833'))     
  
  # Create the HDU
  header = pyfits.Header(cards = cards)
  cols = pyfits.ColDefs(arrays)
  hdu = pyfits.BinTableHDU.from_columns(cols, header = header, name = 'LC ARRAYS')

  return hdu

def ShortCadenceHDU(model):
  '''
  Construct the data HDU file containing the arrays and Kepler observing info.
  
  '''
  
  # Get info from the TPF BinTable HDU Header
  tpf_header = model.meta[1]
  entries = ['WCSN4P', 'WCAX4P', '1CTY4P', '2CTY4P', '1CUN4P', '2CUN4P', '1CRV4P', 
             '2CRV4P', '1CDL4P', '2CDL4P', '1CRP4P', '2CRP4P', 'WCAX4', '1CTYP4', 
             '2CTYP4', '1CRPX4', '2CRPX4', '1CRVL4', '2CRVL4', '1CUNI4', '2CUNI4', 
             '1CDLT4', '2CDLT4', '11PC4', '12PC4', '21PC4', '22PC4', 'WCSN5P', 
             'WCAX5P', '1CTY5P', '2CTY5P', '1CUN5P', '2CUN5P', '1CRV5P', '2CRV5P', 
             '1CDL5P', '2CDL5P', '1CRP5P', '2CRP5P', 'WCAX5', '1CTYP5', '2CTYP5', 
             '1CRPX5', '2CRPX5', '1CRVL5', '2CRVL5', '1CUNI5', '2CUNI5', '1CDLT5', 
             '2CDLT5', '11PC5', '12PC5', '21PC5', '22PC5', 'WCSN6P', 'WCAX6P', 
             '1CTY6P', '2CTY6P', '1CUN6P', '2CUN6P', '1CRV6P', '2CRV6P', '1CDL6P', 
             '2CDL6P', '1CRP6P', '2CRP6P', 'WCAX6', '1CTYP6', '2CTYP6', '1CRPX6', 
             '2CRPX6', '1CRVL6', '2CRVL6', '1CUNI6', '2CUNI6', '1CDLT6', '2CDLT6', 
             '11PC6', '12PC6', '21PC6', '22PC6', 'WCSN7P', 'WCAX7P', '1CTY7P', 
             '2CTY7P', '1CUN7P', '2CUN7P', '1CRV7P', '2CRV7P', '1CDL7P', '2CDL7P', 
             '1CRP7P', '2CRP7P', 'WCAX7', '1CTYP7', '2CTYP7', '1CRPX7', '2CRPX7', 
             '1CRVL7', '2CRVL7', '1CUNI7', '2CUNI7', '1CDLT7', '2CDLT7', '11PC7', 
             '12PC7', '21PC7', '22PC7', 'WCSN8P', 'WCAX8P', '1CTY8P', '2CTY8P', 
             '1CUN8P', '2CUN8P', '1CRV8P', '2CRV8P', '1CDL8P', '2CDL8P', '1CRP8P', 
             '2CRP8P', 'WCAX8', '1CTYP8', '2CTYP8', '1CRPX8', '2CRPX8', '1CRVL8', 
             '2CRVL8', '1CUNI8', '2CUNI8', '1CDLT8', '2CDLT8', '11PC8', '12PC8', 
             '21PC8', '22PC8', 'WCSN9P', 'WCAX9P', '1CTY9P', '2CTY9P', '1CUN9P', 
             '2CUN9P', '1CRV9P', '2CRV9P', '1CDL9P', '2CDL9P', '1CRP9P', '2CRP9P', 
             'WCAX9', '1CTYP9', '2CTYP9', '1CRPX9', '2CRPX9', '1CRVL9', '2CRVL9', 
             '1CUNI9', '2CUNI9', '1CDLT9', '2CDLT9', '11PC9', '12PC9', '21PC9', 
             '22PC9', 'INHERIT', 'EXTNAME', 'EXTVER', 'TELESCOP', 'INSTRUME', 
             'OBJECT', 'KEPLERID', 'RADESYS', 'RA_OBJ', 'DEC_OBJ', 'EQUINOX', 
             'EXPOSURE', 'TIMEREF', 'TASSIGN', 'TIMESYS', 'BJDREFI', 'BJDREFF', 
             'TIMEUNIT', 'TELAPSE', 'LIVETIME', 'TSTART', 'TSTOP', 'LC_START', 
             'LC_END', 'DEADC', 'TIMEPIXR', 'TIERRELA', 'INT_TIME', 'READTIME', 
             'FRAMETIM', 'NUM_FRM', 'TIMEDEL', 'DATE-OBS', 'DATE-END', 'BACKAPP', 
             'DEADAPP', 'VIGNAPP', 'GAIN', 'READNOIS', 'NREADOUT', 'TIMSLICE', 
             'MEANBLCK', 'LCFXDOFF', 'SCFXDOFF']
  cards = []
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*      KEPLER INFO     *'))
  cards.append(('COMMENT', '************************'))
  for entry in entries:
    cards.append(tuple(tpf_header[entry]))
      
  # Add EVEREST info
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*     EVEREST INFO     *'))
  cards.append(('COMMENT', '************************'))
  cards.append(('VERSION', EVEREST_VERSION, 'EVEREST pipeline version'))
  cards.append(('DATE', strftime('%Y-%m-%d'), 'EVEREST file creation date (YYYY-MM-DD)'))
  cards.append(('MODEL', model.name, 'Name of EVEREST model used'))
  cards.append(('APNAME', model.aperture_name, 'Name of aperture used'))
  cards.append(('BPAD', model.bpad, 'Chunk overlap in cadences'))
  cards.append(('BRKPT', model.breakpoints[0], 'Light curve breakpoint'))
  cards.append(('CDIVS', model.cdivs, 'Cross-validation subdivisions'))
  cards.append(('CDPP6', model.cdpp6, 'Average de-trended 6-hr CDPP'))
  cards.append(('CDPP6A', model.cdpp6_arr[0], 'First chunk de-trended 6-hr CDPP'))
  try:
    cards.append(('CDPP6B', model.cdpp6_arr[1], 'Second chunk de-trended 6-hr CDPP'))
  except IndexError:
    cards.append(('CDPP6B', np.nan, 'Second chunk de-trended 6-hr CDPP'))
  cards.append(('CDPP6R', model.cdppr, 'Raw 6-hr CDPP'))
  cards.append(('CDPP6RA', model.cdppr_arr[0], 'First chunk raw 6-hr CDPP'))
  try:
    cards.append(('CDPP6RB', model.cdppr_arr[1], 'Second chunk raw 6-hr CDPP'))
  except IndexError:
    cards.append(('CDPP6RB', np.nan, 'Second chunk raw 6-hr CDPP'))
  cards.append(('CDPP6V', model.cdpp6, 'Average validation 6-hr CDPP'))
  cards.append(('CDPP6VA', model.cdpp6_arr[0], 'First chunk validation 6-hr CDPP'))
  try:
    cards.append(('CDPP6VB', model.cdpp6_arr[1], 'Second chunk validation 6-hr CDPP'))
  except IndexError:
    cards.append(('CDPP6VB', np.nan, 'Second chunk validation 6-hr CDPP'))
  cards.append(('CDPP6G', model.gppp, 'Average GP-de-trended 6-hr CDPP'))
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
  cards.append(('OITER', model.oiter, 'Number of outlier search iterations'))
  cards.append(('OPTGP', model.optimize_gp, 'GP optimization performed?'))
  cards.append(('OSIGMA', model.osigma, 'Outlier tolerance (standard deviations)'))
  cards.append(('PLDORDER', model.pld_order, 'PLD de-trending order'))
  cards.append(('RECRSVE', model.recursive, 'Recursive PLD?'))
  cards.append(('SATUR', model.saturated, 'Is target saturated?'))
  cards.append(('SATAP', model.saturated_aperture_name, 'Saturated aperture name'))
  cards.append(('SATTOL', model.saturation_tolerance, 'Fractional saturation tolerance'))
  
  # Arrays
  sc_quality = np.array(model.sc_quality)
  sc_quality[model.sc_badmask] += QUALITY_BAD
  sc_quality[model.sc_nanmask] += QUALITY_NAN
  sc_quality[model.sc_outmask] += QUALITY_OUT
  arrays = [pyfits.Column(name = 'SC_CADN', format = 'D', array = model.sc_cadn),
            pyfits.Column(name = 'SC_FLUX', format = 'D', array = model.sc_flux, unit = 'e-/s'),
            pyfits.Column(name = 'SC_FRAW', format = 'D', array = model.sc_fraw, unit = 'e-/s'),
            pyfits.Column(name = 'SC_FRAW_ERR', format = 'D', array = model.sc_fraw_err, unit = 'e-/s'),
            pyfits.Column(name = 'SC_QUALITY', format = 'J', array = sc_quality),
            pyfits.Column(name = 'SC_TIME', format = 'D', array = model.sc_time, unit = 'BJD - 2454833')]
  
  # Create the HDU
  header = pyfits.Header(cards = cards)
  cols = pyfits.ColDefs(arrays)
  hdu = pyfits.BinTableHDU.from_columns(cols, header = header, name = 'SC ARRAYS')

  return hdu

def LongCadencePixelHDU(model):
  '''
  Construct the HDU containing the pixel-level light curve.
  
  '''

  # Add EVEREST info
  cards = []
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*     EVEREST INFO     *'))
  cards.append(('COMMENT', '************************'))
  cards.append(('VERSION', EVEREST_VERSION, 'EVEREST pipeline version'))
  cards.append(('DATE', strftime('%Y-%m-%d'), 'EVEREST file creation date (YYYY-MM-DD)'))
  
  # Create the HDU
  header = pyfits.Header(cards = cards)
  
  # The pixel timeseries
  arrays = [pyfits.Column(name = 'FPIX', format = '%dD' % model.fpix.shape[1], array = model.fpix)]
  
  # The first order PLD vectors for all the neighbors (npixels, ncadences)
  X1N = model.XNeighbors[0]
  if X1N is not None:
    arrays.append(pyfits.Column(name = 'X1N', format = '%dD' % X1N.shape[1], array = X1N))
    
  cols = pyfits.ColDefs(arrays)
  hdu = pyfits.BinTableHDU.from_columns(cols, header = header, name = 'LC PIXELS')

  return hdu

def ShortCadencePixelHDU(model):
  '''
  Construct the HDU containing the pixel-level light curve.
  
  '''

  # Add EVEREST info
  cards = []
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*     EVEREST INFO     *'))
  cards.append(('COMMENT', '************************'))
  cards.append(('VERSION', EVEREST_VERSION, 'EVEREST pipeline version'))
  cards.append(('DATE', strftime('%Y-%m-%d'), 'EVEREST file creation date (YYYY-MM-DD)'))
  
  # Create the HDU
  header = pyfits.Header(cards = cards)
  
  # The pixel timeseries
  arrays = [pyfits.Column(name = 'SC_FPIX', format = '%dD' % model.sc_fpix.shape[1], array = model.sc_fpix)]
  
  # The first order PLD vectors for all the neighbors (npixels, ncadences)
  if model.sc_X1N is not None:
    arrays.append(pyfits.Column(name = 'sc_X1N', format = '%dD' % model.sc_X1N.shape[1], array = model.sc_X1N))
    
  cols = pyfits.ColDefs(arrays)
  hdu = pyfits.BinTableHDU.from_columns(cols, header = header, name = 'SC PIXELS')

  return hdu

def ApertureHDU(model):
  '''
  Construct the HDU containing the aperture used to de-trend.
  
  '''
  
  # Get info from the TPF BinTable HDU Header
  tpf_header = model.meta[2]
  entries = ['TELESCOP', 'INSTRUME', 'OBJECT', 'KEPLERID', 'RADESYS', 'RA_OBJ', 
             'DEC_OBJ', 'EQUINOX', 'WCSAXES', 'CTYPE1', 'CTYPE2', 'CRPIX1', 
             'CRPIX2', 'CRVAL1', 'CRVAL2', 'CUNIT1', 'CUNIT2', 'CDELT1', 
             'CDELT2', 'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2', 'WCSNAMEP', 'WCSAXESP', 
             'CTYPE1P', 'CUNIT1P', 'CRPIX1P', 'CRVAL1P', 'CDELT1P', 'CTYPE2P', 
             'CUNIT2P', 'CRPIX2P', 'CRVAL2P', 'CDELT2P', 'NPIXSAP', 'NPIXMISS']
  cards = []
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*      KEPLER INFO     *'))
  cards.append(('COMMENT', '************************'))
  for entry in entries:
    cards.append(tuple(tpf_header[entry]))
  
  # Add EVEREST info
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*     EVEREST INFO     *'))
  cards.append(('COMMENT', '************************'))
  cards.append(('VERSION', EVEREST_VERSION, 'EVEREST pipeline version'))
  cards.append(('DATE', strftime('%Y-%m-%d'), 'EVEREST file creation date (YYYY-MM-DD)'))
  
  # Create the HDU
  header = pyfits.Header(cards = cards)
  hdu = pyfits.ImageHDU(data = model.aperture, header = header, name = 'APERTURE MASK')

  return hdu

def MakeFITS(model):
  '''
  Generate a FITS file for a given :py:mod:`everest` run.
  
  :param model: An :py:mod:`everest` model instance
  
  '''
  
  # Get the fits file name
  outfile = os.path.join(model.dir, 'hlsp_everest_k2_llc_%d-c%02d_kepler_v%s_lc.fits' % 
                        (model.ID, model.season, EVEREST_VERSION))
  if os.path.exists(outfile) and not model.clobber:
    return
  elif os.path.exists(outfile):
    os.remove(outfile)
      
  # Create the HDUs
  primary = PrimaryHDU(model)
  long_cadence = LongCadenceHDU(model)
  lc_pixels = LongCadencePixelHDU(model)
  aperture = ApertureHDU(model)
  
  # Combine to get the HDUList
  if model.has_sc:
    short_cadence = ShortCadenceHDU(model)
    sc_pixels = ShortCadencePixelHDU(model)
    hdulist = pyfits.HDUList([primary, long_cadence, lc_pixels, aperture, short_cadence, sc_pixels])
  else:
    hdulist = pyfits.HDUList([primary, long_cadence, lc_pixels, aperture])
  
  # Output to the FITS file
  hdulist.writeto(outfile)
  
  return