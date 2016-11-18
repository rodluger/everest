#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`user.py` - User-facing routines
----------------------------------------

.. todo:: 
   - Add easy masking
   
'''

from __future__ import division, print_function, absolute_import, unicode_literals
from . import __version__ as EVEREST_VERSION
from . import missions
from .basecamp import Basecamp
from .pld import *
from .gp import GetCovariance
from .config import QUALITY_BAD, QUALITY_NAN, QUALITY_OUT, EVEREST_DEV, EVEREST_FITS
from .utils import InitLog, Formatter
import george
import os, sys, platform
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
import logging
log = logging.getLogger(__name__)

def DownloadFile(ID, mission = 'k2', filename = None, clobber = False):
  '''
  Download a given :py:mod:`everest` file from MAST.
  
  :param bool clobber: If `True`, download and overwrite existing files. Default `False`
  
  '''
  
  # Grab some info
  season = getattr(missions, mission).Season(ID)
  path = getattr(missions, mission).TargetDirectory(ID, season)
  relpath = getattr(missions, mission).TargetDirectory(ID, season, relative = True)
  if filename is None:
    filename = getattr(missions, mission).FITSFile(ID, season)
  
  # Check if file exists
  if not os.path.exists(path):
    os.makedirs(path)
  elif os.path.exists(os.path.join(path, filename)) and not clobber:
    log.info('Found cached file.')
    return os.path.join(path, filename)
  
  # Get file URL
  log.info('Downloading the file...')
  try:
    fitsurl = getattr(missions, mission).FITSUrl(ID, season)
    if not fitsurl.endswith('/'):
      fitsurl += '/'
  except AssertionError:
    fitsurl = None
   
  if (not EVEREST_DEV) and (url is not None):
  
    # Download the data
    r = urllib.request.Request(url + filename)
    handler = urllib.request.urlopen(r)
    code = handler.getcode()
    if int(code) != 200:
      raise Exception("Error code {0} for URL '{1}'".format(code, url + filename))
    data = handler.read()
    
    # Atomically save to disk
    f = NamedTemporaryFile("wb", delete=False)
    f.write(data)
    f.flush()
    os.fsync(f.fileno())
    f.close()
    shutil.move(f.name, os.path.join(path, filename))
    
  else:
    
    # This section is for pre-publication/development use only!
    if EVEREST_FITS is None:
      raise Exception("Unable to locate the file.")
    
    # Get the url
    inpath = os.path.join(EVEREST_FITS, relpath, filename)
    outpath = os.path.join(path, filename)

    # Download the data
    subprocess.call(['scp', inpath, outpath])
  
  # Success?
  if os.path.exists(os.path.join(path, filename)):
    return os.path.join(path, filename)
  else:
    raise Exception("Unable to download the file.")

def ShowDVS(ID, mission = 'k2', model = 'nPLD', clobber = False):
  '''
  
  '''
  
  file = DownloadFile(ID, mission = mission, 
                      filename = model + '.pdf', 
                      clobber = clobber)  
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
    log.info("Unable to open the pdf. The full path is")
    log.info(file)

def Everest(ID, mission = 'k2', quiet = False, clobber = False, **kwargs):
  '''
  
  '''
  
  # Initialize preliminary logging
  if not quiet:
    screen_level = logging.DEBUG
  else:
    screen_level = logging.CRITICAL
  InitLog(None, logging.DEBUG, screen_level, False)

  # Download the FITS file if necessary
  fitsfile = DownloadFile(ID, mission = mission, clobber = clobber)
  model_name = pyfits.getheader(fitsfile, 1)['MODEL']

  class Star(eval(model_name + 'Base'), Basecamp):
    '''
  
    '''
    
    def __repr__(self):
      '''
      
      '''
      
      return "<everest.Star(%d)>" % self.ID
    
    def __init__(self, ID, mission, model_name, fitsfile, quiet, clobber, **kwargs):
      '''
      
      '''
      
      self.ID = ID
      self.mission = mission
      self.model_name = model_name
      self.fitsfile = fitsfile
      self.clobber = clobber
      if not quiet:
        screen_level = logging.DEBUG
      else:
        screen_level = logging.CRITICAL
      log_level = kwargs.get('log_level', logging.DEBUG)
      InitLog(self.logfile, logging.DEBUG, screen_level, False)
      self.download_fits()
      self.load_fits()
      self.init_model()
    
    @property
    def name(self):
      '''
      
      '''
      
      return self.model_name
    
    @property
    def X(self):
      '''
      
      '''
      
      if self._X is None:
        self._X = [None for i in range(self.pld_order)]
        self.get_X()
      return self._X
    
    def download_fits(self):
      '''
      
      '''
      
      pass
    
    def init_model(self):
      '''
      
      '''
      
      log.info("Initializing %s model for %d." % (self.name, self.ID))
      self._A = [[None for i in range(self.pld_order)] for b in self.breakpoints]
      self._B = [[None for i in range(self.pld_order)] for b in self.breakpoints]
      self._mK = [None for b in self.breakpoints]
      self._f = [None for b in self.breakpoints]
      self._X = None
      self._weights = None
      self.K = GetCovariance(self.kernel_params, self.time, self.fraw_err)
      
    def load_fits(self):
      '''
      
      '''
      
      log.info("Loading FITS file for %d." % (self.ID))
      with pyfits.open(self.fitsfile) as f:
        
        # Params and long cadence data
        self.has_sc = self._mission.HasShortCadence(self.ID, season = self.season)
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
      self.saturated_aperture_name = None
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
    
    def plot_aperture(self, show = True):
      '''
      
      '''
      
      # Set up the axes
      fig, ax = pl.subplots(2,2, figsize = (6, 8))
      fig.subplots_adjust(top = 0.975, bottom = 0.025, left = 0.05, 
                          right = 0.95, hspace = 0.05, wspace = 0.05)
      ax = ax.flatten()
      fig.canvas.set_window_title('%s %d' % (self._mission.IDSTRING, self.ID))
      super(Star, self).plot_aperture(ax, labelsize = 12) 
      
      if show:
        pl.show()
        pl.close()
      else:
        return fig, ax
    
    def plot_weights(self, show = True):
      '''
      
      '''
      
      # Set up the axes
      fig = pl.figure(figsize = (12, 12))
      fig.subplots_adjust(top = 0.95, bottom = 0.025, left = 0.1, right = 0.92)
      fig.canvas.set_window_title('%s %d' % (self._mission.IDSTRING, self.ID))
      ax = [pl.subplot2grid((80, 130), (20 * j, 25 * i), colspan = 23, rowspan = 18) 
            for j in range(len(self.breakpoints) * 2) for i in range(1 + 2 * (self.pld_order - 1))]
      cax = [pl.subplot2grid((80, 130), (20 * j, 25 * (1 + 2 * (self.pld_order - 1))), 
             colspan = 4, rowspan = 18) for j in range(len(self.breakpoints) * 2)]
      ax = np.array(ax).reshape(2 * len(self.breakpoints), -1)
      cax = np.array(cax)
      super(Star, self).plot_weights(ax, cax)
      
      if show:
        pl.show()
        pl.close()
      else:
        return fig, ax, cax

    def plot(self, show = True, plot_raw = True, plot_gp = True, 
             plot_bad = True, plot_out = True, plot_sc = False):
      '''
      Plots the final de-trended light curve.
    
      '''

      log.info('Plotting the light curve...')
    
      # Set up axes
      if plot_raw:
        fig, axes = pl.subplots(2, figsize = (13, 9), sharex = True)
        fig.subplots_adjust(hspace = 0.1)
        axes = [axes[1], axes[0]]
        if plot_sc:
          fluxes = [self.sc_flux, self.sc_fraw]
        else:
          fluxes = [self.flux, self.fraw]
        labels = ['EVEREST Flux', 'Raw Flux']
      else:
        fig, axes = pl.subplots(1, figsize = (13, 6))
        axes = [axes]
        if plot_sc:
          fluxes = [self.sc_flux]
        else:
          fluxes = [self.flux]
        labels = ['EVEREST Flux']
      fig.canvas.set_window_title('EVEREST Light curve')
      
      # Set up some stuff
      if plot_sc:
        time = self.sc_time
        badmask = self.sc_badmask
        nanmask = self.sc_nanmask
        outmask = self.sc_outmask
        fraw_err = self.sc_fraw_err
        breakpoints = []
        plot_gp = False
        ms = 2
      else:
        time = self.time
        badmask = self.badmask
        nanmask = self.nanmask
        outmask = self.outmask
        fraw_err = self.fraw_err
        breakpoints = self.breakpoints
        ms = 4
      
      # Get the cdpps
      cdpps = [[self.get_cdpp(self.flux), self.get_cdpp_arr(self.flux)],
               [self.get_cdpp(self.fraw), self.get_cdpp_arr(self.fraw)]]
      self.cdpp6 = cdpps[0][0]
      self.cdpp6_arr = cdpps[0][1]
      
      for n, ax, flux, label, cdpp in zip([0,1], axes, fluxes, labels, cdpps):
        
        # Initialize CDPP
        cdpp6 = cdpp[0]
        cdpp6_arr = cdpp[1]
          
        # Plot the good data points
        ax.plot(self.apply_mask(time, sc = plot_sc), self.apply_mask(flux, sc = plot_sc), ls = 'none', marker = '.', color = 'k', markersize = ms, alpha = 0.5)
    
        # Plot the outliers
        bnmask = np.array(list(set(np.concatenate([badmask, nanmask]))), dtype = int)
        O1 = lambda x: x[outmask]
        O2 = lambda x: x[bnmask]
        if plot_out:
          ax.plot(O1(time), O1(flux), ls = 'none', color = "#777777", marker = '.', markersize = ms, alpha = 0.5)
        if plot_bad:
          ax.plot(O2(time), O2(flux), 'r.', markersize = ms, alpha = 0.25)

        # Plot the GP
        if n == 0 and plot_gp:
          M = lambda x: np.delete(x, bnmask)
          _, amp, tau = self.kernel_params
          gp = george.GP(amp ** 2 * george.kernels.Matern32Kernel(tau ** 2))
          gp.compute(self.apply_mask(time, sc = plot_sc), self.apply_mask(fraw_err, sc = plot_sc))
          med = np.nanmedian(self.apply_mask(flux, sc = plot_sc))
          y, _ = gp.predict(self.apply_mask(flux, sc = plot_sc) - med, time)
          y += med
          ax.plot(M(time), M(y), 'r-', lw = 0.5, alpha = 0.5)

        # Appearance
        if n == 0: 
          ax.set_xlabel('Time (%s)' % self._mission.TIMEUNITS, fontsize = 18)
        ax.set_ylabel(label, fontsize = 18)
        for brkpt in breakpoints[:-1]:
          ax.axvline(time[brkpt], color = 'r', ls = '--', alpha = 0.25)
        if len(cdpp6_arr) == 2:
          ax.annotate('%.2f ppm' % cdpp6_arr[0], xy = (0.02, 0.975), xycoords = 'axes fraction', 
                      ha = 'left', va = 'top', fontsize = 12, color = 'r', zorder = 99)
          ax.annotate('%.2f ppm' % cdpp6_arr[1], xy = (0.98, 0.975), xycoords = 'axes fraction', 
                      ha = 'right', va = 'top', fontsize = 12, color = 'r', zorder = 99)
        else:
          ax.annotate('%.2f ppm' % cdpp6, xy = (0.02, 0.975), xycoords = 'axes fraction', 
                      ha = 'left', va = 'top', fontsize = 12, color = 'r', zorder = 99)
        ax.margins(0.01, 0.1)          
    
        # Get y lims that bound 99% of the flux
        f = np.concatenate([np.delete(f, bnmask) for f in fluxes])
        N = int(0.995 * len(f))
        hi, lo = f[np.argsort(f)][[N,-N]]
        fsort = f[np.argsort(f)]
        pad = (hi - lo) * 0.1
        ylim = (lo - pad, hi + pad)
        ax.set_ylim(ylim)   
        ax.get_yaxis().set_major_formatter(Formatter.Flux)
    
        # Indicate off-axis outliers
        for i in np.where(flux < ylim[0])[0]:
          if i in bnmask:
            color = "#ffcccc"
            if not plot_bad: 
              continue
          elif i in outmask:
            color = "#cccccc"
            if not plot_out:
              continue
          else:
            color = "#ccccff"
          ax.annotate('', xy=(time[i], ylim[0]), xycoords = 'data',
                      xytext = (0, 15), textcoords = 'offset points',
                      arrowprops=dict(arrowstyle = "-|>", color = color))
        for i in np.where(flux > ylim[1])[0]:
          if i in bnmask:
            color = "#ffcccc"
            if not plot_bad:
              continue
          elif i in outmask:
            color = "#cccccc"
            if not plot_out:
              continue
          else:
            color = "#ccccff"
          ax.annotate('', xy=(time[i], ylim[1]), xycoords = 'data',
                      xytext = (0, -15), textcoords = 'offset points',
                      arrowprops=dict(arrowstyle = "-|>", color = color))
      
      # Show total CDPP improvement
      pl.figtext(0.5, 0.94, '%s %d' % (self._mission.IDSTRING, self.ID), fontsize = 18, ha = 'center', va = 'bottom')
      pl.figtext(0.5, 0.905, r'$%.2f\ \mathrm{ppm} \rightarrow %.2f\ \mathrm{ppm}$' % (self.cdppr, self.cdpp6), fontsize = 14, ha = 'center', va = 'bottom')
      
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
      
      '''
      
      ShowDVS(self.ID, mission = self.mission, model = self.model_name, clobber = self.clobber)
    
    def plot_pipeline(self, **kwargs):
      '''
      
      '''
      
      return getattr(missions, mission).pipelines.plot(self.ID, **kwargs)
    
  return Star(ID, mission, model_name, fitsfile, quiet, clobber, **kwargs)