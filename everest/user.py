#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`user.py` - User-facing routines
----------------------------------------
   
'''

from __future__ import division, print_function, absolute_import, unicode_literals
from . import __version__ as EVEREST_VERSION
from . import missions
from .basecamp import Basecamp
from .gp import GetCovariance
from .config import QUALITY_BAD, QUALITY_NAN, QUALITY_OUT, QUALITY_REC, EVEREST_DEV, EVEREST_FITS
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

def DownloadFile(ID, mission = 'k2', cadence = 'lc', filename = None, clobber = False):
  '''
  Download a given :py:mod:`everest` file from MAST.
  
  :param bool clobber: If `True`, download and overwrite existing files. Default `False`
  
  '''
  
  # Grab some info
  season = getattr(missions, mission).Season(ID)
  path = getattr(missions, mission).TargetDirectory(ID, season)
  relpath = getattr(missions, mission).TargetDirectory(ID, season, relative = True)
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

class Everest(Basecamp):
  '''
  
  '''
  
  def __init__(self, ID, mission = 'k2', quiet = False, clobber = False, cadence = 'lc', **kwargs):
    '''
    
    '''
    
    # Read kwargs
    self.ID = ID
    self.mission = mission
    self.clobber = clobber
    
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
    self.fitsfile = DownloadFile(ID, mission = mission, clobber = clobber, cadence = cadence)
    self.model_name = pyfits.getheader(self.fitsfile, 1)['MODEL']
    self.download_fits()
    self.load_fits()
    self._weights = None

  def __repr__(self):
    '''
    
    '''
    
    return "<everest.Everest(%d)>" % self.ID

  @property
  def name(self):
    '''
    
    '''
    
    return self.model_name
  
  def reset(self):
    '''
    
    '''
    
    self.load_fits()
    self._weights = None
  
  def compute(self):
    '''
    
    '''
    
    if self.model_name == 'rPLD':
      self.get_norm()
    super(Everest, self).compute()
  
  def get_norm(self):
    '''
    
    '''
    
    log.info('Computing the PLD normalization...')
    
    # Loop over all chunks
    mod = [None for b in self.breakpoints]
    for b, brkpt in enumerate(self.breakpoints):
      
      # Unmasked chunk
      c = self.get_chunk(b)
      
      # Masked chunk (original mask plus user transit mask)
      inds = np.array(list(set(np.concatenate([self.transitmask, self.recmask]))), dtype = int)
      M = np.delete(np.arange(len(self.time)), inds, axis = 0)
      if b > 0:
        m = M[(M > self.breakpoints[b - 1] - self.bpad) & (M <= self.breakpoints[b] + self.bpad)]
      else:
        m = M[M <= self.breakpoints[b] + self.bpad]

      # This block of the masked covariance matrix
      mK = GetCovariance(self.kernel_params, self.time[m], self.fraw_err[m])
      
      # Get median
      med = np.nanmedian(self.fraw[m])
      
      # Normalize the flux
      f = self.fraw[m] - med
      
      # The X^2 matrices
      A = np.zeros((len(m), len(m)))
      B = np.zeros((len(c), len(m)))
      
      # Loop over all orders
      for n in range(self.pld_order):
        XM = self.X(n,m)
        XC = self.X(n,c)
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
        model = np.concatenate([model, m[self.bpad:-self.bpad] + offset])
  
      # Last chunk
      offset = model[-1] - mod[-1][self.bpad - 1]
      model = np.concatenate([model, mod[-1][self.bpad:] + offset])      
  
    else:

      model = mod[0]
  
    # Subtract the global median
    model -= np.nanmedian(model)

    # Save the norm
    self._norm = self.fraw - model
    
  def download_fits(self):
    '''
    TODO.
    
    '''
    
    pass
        
  def load_fits(self):
    '''
    
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
      self.saturated = f[1].header['SATUR']
      self.saturation_tolerance = f[1].header['SATTOL']
      self.time = f[1].data['TIME']
      self._norm = np.array(self.fraw)
      
      # Chunk arrays
      self.breakpoints = []
      self.cdpp6_arr = []
      self.cdppv_arr = []
      self.cdppr_arr = []
      for c in range(99):
        try:
          self.breakpoints.append(f[1].header['BRKPT%02d' % (c + 1)])
          self.cdpp6_arr.append(f[1].header['CDPP6%02d' % (c + 1)])
          self.cdppr_arr.append(f[1].header['CDPPR%02d' % (c + 1)])
          self.cdppv_arr.append(f[1].header['CDPPV%02d' % (c + 1)])
        except KeyError:
          break
      self.lam = [[f[1].header['LAMB%02d%02d' % (c + 1, o + 1)] for o in range(self.pld_order)] 
                   for c in range(len(self.breakpoints))]
      if self.model_name == 'rPLD':
        self.reclam = [[f[1].header['RECL%02d%02d' % (c + 1, o + 1)] for o in range(self.pld_order)] 
                        for c in range(len(self.breakpoints))]
    # Masks
    self.badmask = np.where(self.quality & 2 ** (QUALITY_BAD - 1))[0]
    self.nanmask = np.where(self.quality & 2 ** (QUALITY_NAN - 1))[0]
    self.outmask = np.where(self.quality & 2 ** (QUALITY_OUT - 1))[0]
    self.recmask = np.where(self.quality & 2 ** (QUALITY_REC - 1))[0]  
      
    # These are not stored in the fits file; we don't need them
    self.saturated_aperture_name = None
    self.apertures = None
    self.Xpos = None
    self.Ypos = None
    self.fpix_err = None
    self.parent_model = None
    self.lambda_arr = None
    self.meta = None
    self.transitmask = np.array([], dtype = int)
  
  def plot_aperture(self, show = True):
    '''
    
    '''
    
    # Set up the axes
    fig, ax = pl.subplots(2,2, figsize = (6, 8))
    fig.subplots_adjust(top = 0.975, bottom = 0.025, left = 0.05, 
                        right = 0.95, hspace = 0.05, wspace = 0.05)
    ax = ax.flatten()
    fig.canvas.set_window_title('%s %d' % (self._mission.IDSTRING, self.ID))
    super(Everest, self).plot_aperture(ax, labelsize = 12) 
    
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
    super(Everest, self).plot_weights(ax, cax)
    
    if show:
      pl.show()
      pl.close()
    else:
      return fig, ax, cax

  def plot(self, show = True, plot_raw = True, plot_gp = True, 
           plot_bad = True, plot_out = True):
    '''
    Plots the final de-trended light curve.
  
    '''

    log.info('Plotting the light curve...')
  
    # Set up axes
    if plot_raw:
      fig, axes = pl.subplots(2, figsize = (13, 9), sharex = True)
      fig.subplots_adjust(hspace = 0.1)
      axes = [axes[1], axes[0]]
      fluxes = [self.flux, self.fraw]
      labels = ['EVEREST Flux', 'Raw Flux']
    else:
      fig, axes = pl.subplots(1, figsize = (13, 6))
      axes = [axes]
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
    self.cdpp6 = cdpps[0][0]
    self.cdpp6_arr = cdpps[0][1]
    
    for n, ax, flux, label, cdpp in zip([0,1], axes, fluxes, labels, cdpps):
      
      # Initialize CDPP
      cdpp6 = cdpp[0]
      cdpp6_arr = cdpp[1]
        
      # Plot the good data points
      ax.plot(self.apply_mask(time), self.apply_mask(flux), ls = 'none', marker = '.', color = 'k', markersize = ms, alpha = 0.5)
  
      # Plot the outliers
      bnmask = np.array(list(set(np.concatenate([badmask, nanmask]))), dtype = int)
      O1 = lambda x: x[outmask]
      O2 = lambda x: x[bnmask]
      O3 = lambda x: x[transitmask]
      if plot_out:
        ax.plot(O1(time), O1(flux), ls = 'none', color = "#777777", marker = '.', markersize = ms, alpha = 0.5)
      if plot_bad:
        ax.plot(O2(time), O2(flux), 'r.', markersize = ms, alpha = 0.25)
      ax.plot(O3(time), O3(flux), 'b.', markersize = ms, alpha = 0.25)
      
      # Plot the GP
      if n == 0 and plot_gp and self.cadence != 'sc':
        M = lambda x: np.delete(x, bnmask)
        _, amp, tau = self.kernel_params
        gp = george.GP(amp ** 2 * george.kernels.Matern32Kernel(tau ** 2))
        gp.compute(self.apply_mask(time), self.apply_mask(fraw_err))
        med = np.nanmedian(self.apply_mask(flux))
        y, _ = gp.predict(self.apply_mask(flux) - med, time)
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
  
  def plot_pipeline(self, *args, **kwargs):
    '''
    
    '''
    
    return getattr(missions, mission).pipelines.plot(self.ID, *args, **kwargs)
  
  def mask_planet(self, t0, period, dur = 0.2):
    '''
    
    '''
    
    mask = []
    t0 += np.ceil((self.time[0] - dur - t0) / period) * period
    for t in np.arange(t0, self.time[-1] + dur, period):
      mask.extend(np.where(np.abs(self.time - t) < dur / 2.)[0])
    self.transitmask = np.array(list(set(np.concatenate([self.transitmask, mask]))))