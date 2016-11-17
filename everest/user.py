#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`user.py` - User-facing routines
----------------------------------------

.. todo:: 
   - Plot raw flux
   - Plot SC flux
   - Add easy masking
   
'''

from __future__ import division, print_function, absolute_import, unicode_literals
from . import __version__ as EVEREST_VERSION
from .basecamp import Basecamp
from .pld import *
from .data import HasShortCadence, Season
from .missions import Missions
from .gp import GetCovariance
from .config import QUALITY_BAD, QUALITY_NAN, QUALITY_OUT
from .utils import InitLog
import os, sys
import numpy as np
import matplotlib.pyplot as pl
try:
  import pyfits
except ImportError:
  try:
    import astropy.io.fits as pyfits
  except ImportError:
    raise Exception('Please install the `pyfits` package.')
import logging
log = logging.getLogger(__name__)

def Everest(ID, mission = 'k2', quiet = False, **kwargs):
  '''
  
  '''

  # Grab some info
  season = Season(ID, mission)
  target_dir = Missions[mission].TargetDirectory(ID, season)
  fitsfile = os.path.join(target_dir, 'hlsp_everest_k2_llc_%d-c%02d_kepler_v%s_lc.fits' % 
                         (ID, season, EVEREST_VERSION))
  mission = pyfits.getheader(fitsfile, 0)['MISSION']
  model_name = pyfits.getheader(fitsfile, 1)['MODEL']

  class Star(eval(model_name + 'Base'), Basecamp):
    '''
  
    '''
    
    def __repr__(self):
      '''
      
      '''
      
      return "<everest.Star(%d)>" % self.ID
    
    def __init__(self, ID, mission, model_name, fitsfile, quiet, **kwargs):
      '''
      
      '''
      
      self.ID = ID
      self.mission = mission
      self.model_name = model_name
      self.fitsfile = fitsfile
      if not quiet:
        screen_level = logging.DEBUG
      else:
        screen_level = logging.CRITICAL
      log_level = kwargs.get('log_level', logging.DEBUG)
      InitLog(self.logfile, logging.DEBUG, screen_level, False)
      log.info("Loading FITS file for %d." % (self.ID))
      self.load_fits()
      log.info("Initializing %s model for %d." % (self.name, self.ID))
      self.init_model()
    
    @property
    def name(self):
      '''
      
      '''
      
      return self.model_name
    
    def init_model(self):
      '''
      
      '''
      
      self._A = [[None for i in range(self.pld_order)] for b in self.breakpoints]
      self._B = [[None for i in range(self.pld_order)] for b in self.breakpoints]
      self._mK = [None for b in self.breakpoints]
      self._f = [None for b in self.breakpoints]
      self._X = [None for i in range(self.pld_order)]
      self._weights = None
      self.K = GetCovariance(self.kernel_params, self.time, self.fraw_err)
      self.get_X()
      
    def load_fits(self):
      '''
      
      '''
      
      with pyfits.open(self.fitsfile) as f:
        
        # Params and long cadence data
        self.has_sc = HasShortCadence(self.ID, season = self.season)
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
        self.cadn = f[1].data['CADN']
        self.cdivs = f[1].header['CDIVS']
        self.cdpp6 = f[1].header['CDPP6']
        self.cdppr = f[1].header['CDPPR']
        self.cdppv = f[1].header['CDPPV']
        self.gppp = f[1].header['CDPPG']
        self.fpix = f[2].data['FPIX']
        self.pixel_images = [f[4].data['STAMP1'], f[4].data['STAMP2'], f[4].data['STAMP3']]
        self.fraw = f[1].data['FRAW']
        self.fraw_err = f[1].data['FRAW_ERR']
        self.giter = f[1].header['GITER']
        self.gp_factor = f[1].header['GPFACTOR']
        self.hires = f[5].data
        self.kernel_params = np.array([f[1].header['GPWHITE'], 
                                       f[1].header['GPRED'], 
                                       f[1].header['GPTAU']])
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
            self.nearby.append({'ID': id, 'x': x, 'y': y, 'mag': mag, 'x0': x0, 'y0': y0})
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
        self.quality = f[1].data['QUALITY']
        self.recursive = f[1].header['RECRSVE']
        self.saturated = f[1].header['SATUR']
        self.saturated_aperture_name = f[1].header['SATAP']
        self.saturation_tolerance = f[1].header['SATTOL']
        self.time = f[1].data['TIME']
        
        # Short cadence
        if self.has_sc:
          try:
            self.sc_X1N = f[7].data['SC_X1N']
          except KeyError:
            self.sc_X1N = None
          try:
            self.sc_bkg = f[6].data['SC_BKG']
          except KeyError:
            self.sc_bkg = 0.
          self.sc_bpad = f[6].header['SC_BPAD']
          self.sc_cadn = f[6].data['SC_CADN']
          self.sc_fpix = f[7].data['SC_FPIX']
          self.sc_fraw = f[6].data['SC_FRAW']
          self.sc_fraw_err = f[6].data['SC_FERR']
          self.sc_model = self.sc_fraw - f[6].data['SC_FLUX']
          self.sc_quality = f[6].data['SC_QUAL']
          self.sc_time = f[6].data['SC_TIME']
        else:
          self.sc_X1N = None
          self.sc_bkg = None
          self.sc_bpad = None
          self.sc_cadn = None
          self.sc_time = None
          self.sc_fpix = None
          self.sc_fraw = None
          self.sc_fraw_err = None
          self.sc_model = None
          self.sc_quality = None
          self.sc_bkg = None
          
        # Chunk arrays
        self.breakpoints = []
        self.cdpp6_arr = []
        self.cdppv_arr = []
        self.cdppr_arr = []
        for c in range(99):
          try:
            self.breakpoints.append(f[1].header['BRKPT%d' % (c + 1)])
            self.cdpp6_arr.append(f[1].header['CDPP6%02d' % (c + 1)])
            self.cdppr_arr.append(f[1].header['CDPPR%02d' % (c + 1)])
            self.cdppv_arr.append(f[1].header['CDPPV%02d' % (c + 1)])
          except KeyError:
            break
        self.lam = [[f[1].header['LAMBDA%d%d' % (c + 1, o + 1)] for o in range(self.pld_order)] 
                     for c in range(len(self.breakpoints))]
        
      # Masks
      self.badmask = np.where(self.quality & 2 ** (QUALITY_BAD - 1))[0]
      self.nanmask = np.where(self.quality & 2 ** (QUALITY_NAN - 1))[0]
      self.outmask = np.where(self.quality & 2 ** (QUALITY_OUT - 1))[0]
      if self.has_sc:
        self.sc_badmask = np.where(self.sc_quality & 2 ** (QUALITY_BAD - 1))[0]
        self.sc_nanmask = np.where(self.sc_quality & 2 ** (QUALITY_NAN - 1))[0]
        self.sc_outmask = np.where(self.sc_quality & 2 ** (QUALITY_OUT - 1))[0]
      else:
        self.sc_badmask = None
        self.sc_nanmask = None
        self.sc_outmask = None
        
      # These are not stored in the fits file; we don't need them
      self.apertures = None
      self.Xpos = None
      self.Ypos = None
      self.fpix_err = None
      self.parent_model = None
      self.lambda_arr = None
      self.sc_Xpos = None
      self.sc_Ypos = None
      self.sc_fpix_err = None
      self.meta = None
      self.sc_meta = None
      self.transitmask = np.array([], dtype = int)
      self.sc_transitmask = np.array([], dtype = int)

  return Star(ID, mission, model_name, fitsfile, quiet, **kwargs)