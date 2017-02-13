#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`user.py` - User Python routines
----------------------------------------

This is the gateway to the :py:obj:`everest` catalog, containing
all of the user-facing code.

- :py:class:`Everest` is the main user-facing class for interfacing with the catalog
- :py:func:`DVS` downloads and plots the data validation summary for a given target

Instantiating an :py:class:`Everest` class automatically downloads the light curve
from the online MAST catalog. So, to get started, all you need to do is run

.. code-block :: python

   import everest
   star = everest.Everest(201367065)

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from . import __version__ as EVEREST_VERSION
from . import __subversion__ as EVEREST_SUBVERSION
from . import missions
from .basecamp import Basecamp
from .detrender import pPLD
from .gp import GetCovariance
from .config import QUALITY_BAD, QUALITY_NAN, QUALITY_OUT, QUALITY_REC, QUALITY_TRN, EVEREST_DEV, EVEREST_FITS
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
import six
from six.moves import urllib
from tempfile import NamedTemporaryFile
import shutil
from distutils.version import LooseVersion
import logging
log = logging.getLogger(__name__)

def DownloadFile(ID, mission = 'k2', cadence = 'lc', filename = None, clobber = False, **kwargs):
  '''
  Download a given :py:mod:`everest` file from MAST.
  
  :param str mission: The mission name. Default `k2`
  :param str cadence: The light curve cadence. Default `lc`
  :param str filename: The name of the file to download. Default :py:obj:`None`, in which case the default \
                       FITS file is retrieved.
  :param bool clobber: If :py:obj:`True`, download and overwrite existing files. Default :py:obj:`False`
  :param int season: The observing season (campaign) number
  
  '''
  
  # Grab some info
  season = kwargs.get('season', getattr(missions, mission).Season(ID))
  if season is None:
    if getattr(missions, mission).ISTARGET(ID):
      raise ValueError("Target not found in local database. Consider upgrading EVEREST by running `pip install everest --upgrade`.")
    else:
      raise ValueError('Invalid target ID.')
  elif hasattr(season, '__len__'):
    raise ValueError("Multiple seasons available for this target. Please specify the desired one with the `season` kwarg.")
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
    log.error("Error code {0} for URL '{1}'".format(code, fitsurl + filename))
    
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
    raise Exception("Unable to download the file.")

def DVS(ID, mission = 'k2', clobber = False, cadence = 'lc', model = 'nPLD', **kwargs):
  '''
  Show the data validation summary (DVS) for a given target.
  
  :param str mission: The mission name. Default `k2`
  :param str cadence: The light curve cadence. Default `lc`
  :param bool clobber: If :py:obj:`True`, download and overwrite existing files. Default :py:obj:`False`
  :param int season: The observing season (campaign) number
  
  '''
  
  # Get season
  season = kwargs.get('season', getattr(missions, mission).Season(ID))
  if hasattr(season, '__len__'):
    raise ValueError("Multiple seasons available for this target. Please specify the desired one with the `season` kwarg.")
  
  # Get file name
  if model == 'nPLD':
    filename = getattr(missions, mission).DVSFile(ID, season, cadence)
  else:
    if cadence == 'sc':
      filename = model + '.sc.pdf'
    else:
      filename = model + '.pdf'
    
  file = DownloadFile(ID, mission = mission, 
                      filename = filename, 
                      clobber = clobber, season = season)  
  
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
  The main user-accessible :py:mod:`everest` class for interfacing with the light curves
  stored on MAST. Instantiating this class downloads the current :py:mod:`everest` FITS 
  file for the requested target and populates the class instance with the light curve
  data and attributes. Many of the methods are inherited from :py:class:`everest.Basecamp`.
  
  :param int ID: The target ID. For `k2`, this is the `EPIC` number of the star.
  :param str mission: The mission name. Default `k2`
  :param bool quiet: Suppress :py:obj:`stdout` messages? Default :py:obj:`False`
  :param str cadence: The light curve cadence. Default `lc`
  :param bool clobber: If :py:obj:`True`, download and overwrite existing files. Default :py:obj:`False`
  :param int season: The observing season (campaign)
  
  '''
  
  def __init__(self, ID, mission = 'k2', quiet = False, clobber = False, cadence = 'lc', **kwargs):
    '''
    
    '''
    
    # Read kwargs
    self.ID = ID
    self.mission = mission
    self.clobber = clobber
    self._season = kwargs.get('season', getattr(missions, mission).Season(ID))
    if hasattr(self._season, '__len__'):
      raise ValueError("Multiple seasons available for this target. Please specify the desired one with the `season` kwarg.")
    
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
    self.fitsfile = DownloadFile(ID, mission = mission, clobber = clobber, cadence = cadence, season = self.season)
    self.model_name = pyfits.getheader(self.fitsfile, 1)['MODEL']
    self._weights = None
    
    # Check the pipeline version. Do we need to upgrade?
    subversion = pyfits.getheader(self.fitsfile, 1).get('SUBVER', None)
    if subversion is not None:
      if LooseVersion(subversion) > LooseVersion(EVEREST_SUBVERSION):
        raise Exception("Desired light curve was generated with EVEREST version %s, but current version is %s.\n" % (subversion, EVEREST_SUBVERSION) +
                        "Please upgrade EVEREST by running `pip install everest --upgrade`.")
    
    # Load the FITS file
    self.load_fits()

  def __repr__(self):
    '''
    
    '''
    
    return "<everest.Everest(%d)>" % self.ID

  @property
  def name(self):
    '''
    Returns the name of the :py:mod:`everest` model used to generate this light curve.
    
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
    Re-compute the :py:mod:`everest` model for the given value of :py:obj:`lambda`.
    For long cadence `k2` light curves, this should take several seconds. For short
    cadence `k2` light curves, it may take a few minutes.
    Note that this is a simple wrapper around :py:func:`everest.Basecamp.compute`.
    
    '''
    
    # If we're doing iterative PLD, get the normalization
    if self.model_name == 'iPLD':
      self._get_norm()
    
    # Compute as usual
    super(Everest, self).compute()
    
    # Make NaN cadences NaNs
    if self.nsub == 1:
      self.flux[self.nanmask] = np.nan
    else:
      for k in range(self.nsub):
        self.flux[k][self.nanmask[k]] = np.nan
    
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
        
  def load_fits(self):
    '''
    Load the FITS file from disk and populate the class instance with its data.
    
    '''
    
    log.info("Loading FITS file for %d." % (self.ID))
    with pyfits.open(self.fitsfile) as f:
      
      # Params and long cadence data
      self.loaded = True
      self.is_parent = False
      
      # Get the number of subseasons
      n = 0
      for c in range(99):
        n = f[1].header.get('SUBSN%02d' % (c + 1), n)
      self.nsub = n + 1
      
      # Get the HDU numbers
      PRIMARY = 0
      LIGHTCURVE = 1
      if self.nsub == 1:
        PIXELS = 2
        APERTURE = 3
        IMAGES = 4
        HIRES = 5
      else:
        PIXELS = list(range(2, 2 + self.nsub))
        APERTURE = list(range(2 + self.nsub, 2 + 2 * self.nsub))
        IMAGES = 2 + 2 * self.nsub
        HIRES = 2 + 2 * self.nsub + 1

      # Get the model params
      self.bpad = f[LIGHTCURVE].header['BPAD']
      self.cbv_minstars = []
      self.cbv_num = f[LIGHTCURVE].header['CBVNUM']
      self.cbv_niter = f[LIGHTCURVE].header['CBVNITER']
      self.cbv_win = f[LIGHTCURVE].header['CBVWIN']
      self.cbv_order = f[LIGHTCURVE].header['CBVORD']
      self.cdivs = f[LIGHTCURVE].header['CDIVS']
      self.cdpp = f[LIGHTCURVE].header['CDPP']
      self.cdppr = f[LIGHTCURVE].header['CDPPR']
      self.cdppv = f[LIGHTCURVE].header['CDPPV']
      self.cdppg = f[LIGHTCURVE].header['CDPPG']
      self.cv_min = f[LIGHTCURVE].header['CVMIN']
      self.giter = f[LIGHTCURVE].header['GITER']
      self.gmaxf = f[LIGHTCURVE].header.get('GMAXF', 200)
      self.gp_factor = f[LIGHTCURVE].header['GPFACTOR']
      self.hires = f[5].data
      self.kernel_params = np.array([f[LIGHTCURVE].header['GPWHITE'], 
                                     f[LIGHTCURVE].header['GPRED'], 
                                     f[LIGHTCURVE].header['GPTAU']])
      self.pld_order = f[LIGHTCURVE].header['PLDORDER']
      self.lam_idx = self.pld_order
      self.leps = f[LIGHTCURVE].header['LEPS']
      self.mag = f[PRIMARY].header['KEPMAG']
      self.max_pixels = f[LIGHTCURVE].header['MAXPIX']
      self.oiter = f[LIGHTCURVE].header['OITER']
      self.optimize_gp = f[LIGHTCURVE].header['OPTGP']
      self.osigma = f[LIGHTCURVE].header['OSIGMA']
      self.saturated = f[LIGHTCURVE].header['SATUR']
      self.saturation_tolerance = f[LIGHTCURVE].header['SATTOL']
      
      # Get the sub-season/segment-specific stuff
      self.breakpoints = []
      self.cdpp_arr = []
      self.cdppv_arr = []
      self.cdppr_arr = []
      subseason = []
      for c in range(99):
        try:
          self.breakpoints.append(int(f[LIGHTCURVE].header['BRKPT%02d' % (c + 1)]))
          self.cdpp_arr.append(f[LIGHTCURVE].header['CDPP%02d' % (c + 1)])
          self.cdppr_arr.append(f[LIGHTCURVE].header['CDPPR%02d' % (c + 1)])
          self.cdppv_arr.append(f[LIGHTCURVE].header['CDPPV%02d' % (c + 1)])
          if f[LIGHTCURVE].header.get('SUBSN%02d' % (c + 1), None) is not None:
            subseason.append(f[LIGHTCURVE].header['SUBSN%02d' % (c + 1)])
        except KeyError:
          break
      if self.nsub == 1:
        self.nseg = len(self.breakpoints)
      else:
        self.breakpoints = np.split(self.breakpoints, np.where(np.diff(subseason) != 0)[0] + 1)
        self.nseg = len([item for sublist in self.breakpoints for item in sublist])
      self.lam = [[f[LIGHTCURVE].header['LAMB%02d%02d' % (c + 1, o + 1)] for o in range(self.pld_order)] 
                   for c in range(len(self._breakpoints))]
      if self.model_name == 'iPLD':
        self.reclam = [[f[LIGHTCURVE].header['RECL%02d%02d' % (c + 1, o + 1)] for o in range(self.pld_order)] 
                        for c in range(len(self._breakpoints))]
      
      # Get the aperture(s)
      if self.nsub == 1:
        self.aperture = f[APERTURE].data
        self.aperture_name = f[LIGHTCURVE].header['APNAME']
      else:
        self.aperture = [f[n].data for n in APERTURE]
        self.aperture_name = [f[LIGHTCURVE].header['APNAME%02d' % (k + 1)] for k in range(self.nsub)]
      
      # Get the pixel images
      self.pixel_images = [f[IMAGES].data['STAMP1'], f[IMAGES].data['STAMP2'], f[IMAGES].data['STAMP3']]
      
      # Get the neighboring PLD vectors
      try:
        if self.nsub == 1:
          self.X1N = f[PIXELS].data['X1N']
        else:
          self.X1N = [f[n].data['X1N'] for n in PIXELS]
      except KeyError:
        self.X1N = None
      
      # Get the masked planets
      self.planets = []
      for i in range(99):
        try:
          t0 = f[LIGHTCURVE].header['P%02dT0' % (i + 1)]
          per = f[LIGHTCURVE].header['P%02dPER' % (i + 1)]
          dur = f[LIGHTCURVE].header['P%02dDUR' % (i + 1)]
          self.planets.append((t0, per, dur))
        except KeyError:
          break
      
      # Get the neighboring targets
      self.nearby = []
      for i in range(99):
        try:
          ID = f[LIGHTCURVE].header['NRBY%02dID' % (i + 1)]
          x = f[LIGHTCURVE].header['NRBY%02dX' % (i + 1)]
          y = f[LIGHTCURVE].header['NRBY%02dY' % (i + 1)]
          mag = f[LIGHTCURVE].header['NRBY%02dM' % (i + 1)]
          x0 = f[LIGHTCURVE].header['NRBY%02dX0' % (i + 1)]
          y0 = f[LIGHTCURVE].header['NRBY%02dY0' % (i + 1)]
          self.nearby.append({'ID': ID, 'x': x, 'y': y, 'mag': mag, 'x0': x0, 'y0': y0})
        except KeyError:
          break
      self.neighbors = []
      for c in range(99):
        try:
          self.neighbors.append(f[LIGHTCURVE].header['NEIGH%02d' % (c + 1)])
        except KeyError:
          break

      # 1D Arrays
      try:
        self.bkg = f[LIGHTCURVE].data['BKG']
      except KeyError:
        self.bkg = 0.
      self.cadn = f[LIGHTCURVE].data['CADN']
      self.fraw = f[LIGHTCURVE].data['FRAW']
      self.fraw_err = f[LIGHTCURVE].data['FRAW_ERR']
      self.model = self.fraw - f[LIGHTCURVE].data['FLUX']
      self.quality = f[LIGHTCURVE].data['QUALITY']
      self.time = f[LIGHTCURVE].data['TIME']
      self._norm = np.array(self.fraw)
      
      # CBVs
      self.XCBV = np.empty((len(self.time), 0))
      for i in range(99):
        try:
          self.XCBV = np.hstack([self.XCBV, f[LIGHTCURVE].data['CBV%02d' % (i + 1)].reshape(-1, 1)])
        except KeyError:
          break
          
      # Split into sub-seasons
      if self.nsub > 1:
        split = np.cumsum([b[-1] + 1 for b in self.breakpoints[:-1]])
        if self.bkg == 0:
          self.bkg = [0. for k in range(self.nsub)]
        else:
          self.bkg = np.split(self.bkg, split)
        self.cadn = np.split(self.cadn, split)
        self.fraw = np.split(self.fraw, split)
        self.fraw_err = np.split(self.fraw_err, split)
        self.model = np.split(self.model, split)
        self.quality = np.split(self.quality, split)
        self.time = np.split(self.time, split)
        self._norm = np.split(self._norm, split)
        self.XCBV = np.split(self.XCBV, split, axis = 0)

      # 2D Arrays
      if self.nsub == 1:
        self.fpix = f[PIXELS].data['FPIX']
      else:
        self.fpix = [f[n].data['FPIX'] for n in PIXELS]

      # Masks
      if self.nsub == 1:
        self.badmask = np.where(self.quality & 2 ** (QUALITY_BAD - 1))[0]
        self.nanmask = np.where(self.quality & 2 ** (QUALITY_NAN - 1))[0]
        self.outmask = np.where(self.quality & 2 ** (QUALITY_OUT - 1))[0]
        self.recmask = np.where(self.quality & 2 ** (QUALITY_REC - 1))[0]  
        self.transitmask = np.where(self.quality & 2 ** (QUALITY_TRN - 1))[0]
      else:
        self.badmask = [np.where(q & 2 ** (QUALITY_BAD - 1))[0] for q in self.quality]
        self.nanmask = [np.where(q & 2 ** (QUALITY_NAN - 1))[0] for q in self.quality]
        self.outmask = [np.where(q & 2 ** (QUALITY_OUT - 1))[0] for q in self.quality]
        self.recmask = [np.where(q & 2 ** (QUALITY_REC - 1))[0] for q in self.quality]  
        self.transitmask = [np.where(q & 2 ** (QUALITY_TRN - 1))[0] for q in self.quality]

    # These are not stored in the fits file; we don't need them
    self.saturated_aperture_name = None
    self.apertures = None
    self.Xpos = None
    self.Ypos = None
    self.fpix_err = None
    self.parent_model = None
    self.lambda_arr = None
    self.meta = None
  
  def plot_aperture(self, show = True):
    '''
    Plot sample postage stamps for the target with the aperture outline marked,
    as well as a high-res target image (if available).
    
    :param bool show: Show the plot or return the `(fig, ax)` instance? Default :py:obj:`True`
    
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

  def plot(self, show = True, plot_raw = True, plot_gp = True, 
           plot_bad = True, plot_out = True, plot_cbv = True,
           simple = False):
    '''
    Plots the final de-trended light curve.
    
    :param bool show: Show the plot or return the `(fig, ax)` instance? Default :py:obj:`True`
    :param bool plot_raw: Show the raw light curve? Default :py:obj:`True`
    :param bool plot_gp: Show the GP model prediction? Default :py:obj:`True`
    :param bool plot_bad: Show and indicate the bad data points? Default :py:obj:`True`
    :param bool plot_out: Show and indicate the outliers? Default :py:obj:`True`
    :param bool plot_cbv: Plot the CBV-corrected light curve? Default :py:obj:`True`. If \
                          :py:obj:`False`, plots the de-trended but uncorrected light curve.
  
    '''

    log.info('Plotting the light curve...')
    
    # Set up axes
    if plot_raw:
      fig, axes = pl.subplots(2, figsize = (11, 7), sharex = True)
      fig.subplots_adjust(hspace = 0.1)
      axes = [axes[1], axes[0]]
      if plot_cbv:
        fluxes = [self.fcor, self.fraw]
      else:
        fluxes = [self.flux, self.fraw]
      labels = ['EVEREST Flux', 'Raw Flux']
    else:
      fig, axes = pl.subplots(1, figsize = (13, 6))
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
    f = np.array(self.flux)
    fr = np.array(self.fraw)
    if self.nsub > 1:
      for k in range(self.nsub):
        f[k] = self.apply_mask(f[k], k = k)
        fr[k] = self.apply_mask(fr[k], k = k)
      f = np.concatenate(f)
      fr = np.concatenate(fr)
    else:
      f = self.apply_mask(f)
      fr = self.apply_mask(fr)
    cdpps = [[self.get_cdpp(f), self.get_cdpp_arr(self.flux)],
             [self.get_cdpp(fr), self.get_cdpp_arr(self.fraw)]]
    self.cdpp = cdpps[0][0]
    self.cdpp_arr = cdpps[0][1]
    
    # Get y lims that bound 99% of the flux
    if self.nsub == 1:
      f = np.concatenate([np.delete(f, np.array(list(set(np.concatenate([badmask, nanmask]))), dtype = int)) for f in fluxes])
    else:
      f = np.concatenate([np.delete(f[k], np.array(list(set(np.concatenate([badmask[k], nanmask[k]]))), dtype = int)) for k in range(self.nsub) for f in fluxes])
    N = int(0.995 * len(f))
    hi, lo = f[np.argsort(f)][[N,-N]]
    pad = (hi - lo) * 0.1
    ylim = (lo - pad, hi + pad)

    # Loop over axes
    for n, ax, flux, label, c in zip([0,1], axes, fluxes, labels, cdpps):
    
      # Loop over all sub-seasons
      for k in range(self.nsub):
      
        # No sub-seasons?
        if self.nsub == 1:
          k = slice(None, None, None)
      
        # Initialize CDPP
        cdpp = c[0]
        cdpp_arr = c[1]
        
        # Plot the good data points
        ax.plot(self.apply_mask(time[k], k = k), self.apply_mask(flux[k], k = k), ls = 'none', marker = '.', color = 'k', markersize = ms, alpha = 0.5)
  
        # Plot the outliers
        bnmask = np.array(list(set(np.concatenate([badmask[k], nanmask[k]]))), dtype = int)
        bmask = [i for i in self.badmask[k] if i not in self.nanmask[k]]
        O1 = lambda x: x[outmask[k]]
        O2 = lambda x: x[bmask]
        O3 = lambda x: x[transitmask[k]]

        if plot_out:
          ax.plot(O1(time[k]), O1(flux[k]), ls = 'none', color = "#777777", marker = '.', markersize = ms, alpha = 0.5)
        if plot_bad:
          ax.plot(O2(time[k]), O2(flux[k]), 'r.', markersize = ms, alpha = 0.25)
        ax.plot(O3(time[k]), O3(flux[k]), 'b.', markersize = ms, alpha = 0.25)
      
        # Plot the GP
        if n == 0 and plot_gp and self.cadence != 'sc' and (len(self.mask[k]) < len(self.time[k])):
          _, amp, tau = self.kernel_params
          gp = george.GP(amp ** 2 * george.kernels.Matern32Kernel(tau ** 2))
          gp.compute(self.apply_mask(time[k], k = k), self.apply_mask(fraw_err[k], k = k))
          med = np.nanmedian(self.apply_mask(flux[k], k = k))
          y, _ = gp.predict(self.apply_mask(flux[k], k = k) - med, time[k])
          y += med
          ax.plot(self.apply_mask(time[k], k = k), self.apply_mask(y, k = k), 'r-', lw = 0.5, alpha = 0.5)

        # Appearance
        if (self.nsub == 1 or k == 0):
          if n == 0: 
            ax.set_xlabel('Time (%s)' % self._mission.TIMEUNITS, fontsize = 18)
          ax.set_ylabel(label, fontsize = 18)
          for brkpt in breakpoints[k][:-1]:
            ax.axvline(time[brkpt], color = 'r', ls = '--', alpha = 0.25)
          if len(cdpp_arr) == 2:
            if not np.isnan(cdpp_arr[0]):
              ax.annotate('%.2f ppm' % cdpp_arr[0], xy = (0.02, 0.975), xycoords = 'axes fraction', 
                          ha = 'left', va = 'top', fontsize = 12, color = 'r', zorder = 99)
            if not np.isnan(cdpp_arr[1]):
              ax.annotate('%.2f ppm' % cdpp_arr[1], xy = (0.98, 0.975), xycoords = 'axes fraction', 
                          ha = 'right', va = 'top', fontsize = 12, color = 'r', zorder = 99)
          elif len(cdpp_arr) < 6:
            for n in range(len(cdpp_arr)):
              if n > 0:
                x = (self.time[self.breakpoints[k][n - 1]] - self.time[k][0]) / (self.time[k][-1] - self.time[k][0]) + 0.02
              else:
                x = 0.02
              if not np.isnan(cdpp_arr[n]):
                ax.annotate('%.2f ppm' % cdpp_arr[n], xy = (x, 0.975), xycoords = 'axes fraction', 
                            ha = 'left', va = 'top', fontsize = 10, zorder = 99, color = 'r')
          else:
            ax.annotate('%.2f ppm' % cdpp, xy = (0.02, 0.975), xycoords = 'axes fraction', 
                        ha = 'left', va = 'top', fontsize = 12, color = 'r', zorder = 99)
          ax.margins(0.01, 0.1)          

        # Indicate off-axis outliers
        for i in np.where(flux[k] < ylim[0])[0]:
          if i in bmask:
            color = "#ffcccc"
            if not plot_bad: 
              continue
          elif i in outmask[k]:
            color = "#cccccc"
            if not plot_out:
              continue
          elif i in nanmask[k]:
            continue
          else:
            color = "#ccccff"
          ax.annotate('', xy=(time[k][i], ylim[0]), xycoords = 'data',
                      xytext = (0, 15), textcoords = 'offset points',
                      arrowprops=dict(arrowstyle = "-|>", color = color))
        for i in np.where(flux[k] > ylim[1])[0]:
          if i in bmask:
            color = "#ffcccc"
            if not plot_bad:
              continue
          elif i in outmask[k]:
            color = "#cccccc"
            if not plot_out:
              continue
          elif i in nanmask[k]:
            continue
          else:
            color = "#ccccff"
          ax.annotate('', xy=(time[k][i], ylim[1]), xycoords = 'data',
                      xytext = (0, -15), textcoords = 'offset points',
                      arrowprops=dict(arrowstyle = "-|>", color = color))
      
      ax.set_ylim(ylim)   
      ax.get_yaxis().set_major_formatter(Formatter.Flux)
      
    # Show total CDPP improvement
    pl.figtext(0.5, 0.94, '%s %d' % (self._mission.IDSTRING, self.ID), fontsize = 18, ha = 'center', va = 'bottom')
    pl.figtext(0.5, 0.905, r'$%.2f\ \mathrm{ppm} \rightarrow %.2f\ \mathrm{ppm}$' % (self.cdppr, self.cdpp), fontsize = 14, ha = 'center', va = 'bottom')
    
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
    
    DVS(self.ID, mission = self.mission, model = self.model_name, clobber = self.clobber, season = self.season)
  
  def plot_pipeline(self, pipeline, *args, **kwargs):
    '''
    Plots the light curve for the target de-trended with a given pipeline.
    
    :param str pipeline: The name of the pipeline (lowercase). Options are \
           'everest2', 'everest1', and other mission-specific pipelines. For \
           `K2`, the available pipelines are 'k2sff' and 'k2sc'.
    
    Additional :py:obj:`args` and :py:obj:`kwargs` are passed directly to
    the :py:func:`pipelines.plot` function of the mission.
    
    '''
    
    if pipeline != 'everest2':
      kwargs.update({'season', self.season})
      return getattr(missions, self.mission).pipelines.plot(self.ID, pipeline, *args, **kwargs)
    
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
      if self.nsub == 1:
        bnmask = np.concatenate([self.nanmask, self.badmask])
        time = np.delete(self.time, bnmask)
        flux = np.delete(y, bnmask)
      else:
        bnmask = [np.concatenate([self.nanmask[k], self.badmask[k]]) for k in range(self.nsub)]
        time = [np.delete(self.time[k], bnmask[k]) for k in range(self.nsub)]
        flux = [np.delete(y[k], bnmask[k]) for k in range(self.nsub)]
  
      # Plot it
      fig, ax = pl.subplots(1, figsize = (10, 4))
      fig.subplots_adjust(bottom = 0.15)
      for k in range(self.nsub):
        if self.nsub == 1: k = slice(None, None, None)
        ax.plot(time[k], flux[k], "k.", markersize = 3, alpha = 0.5)
  
      # Axis limits
      if self.nsub == 1:
        N = int(0.995 * len(flux))
        hi, lo = flux[np.argsort(flux)][[N,-N]]
      else:
        N = int(0.995 * len(np.concatenate(flux)))
        hi, lo = np.concatenate(flux)[np.argsort(np.concatenate(flux))][[N,-N]]
      pad = (hi - lo) * 0.1
      ylim = (lo - pad, hi + pad)
      ax.set_ylim(ylim)
      
      # Plot bad data points
      for k in range(self.nsub):
        if self.nsub == 1: k = slice(None, None, None)
        ax.plot(self.time[k][self.badmask[k]], y[k][self.badmask[k]], "r.", markersize = 3, alpha = 0.2)
      
      # Show the CDPP
      ax.annotate('%.2f ppm' % self.get_cdpp(), 
                  xy = (0.98, 0.975), xycoords = 'axes fraction', 
                  ha = 'right', va = 'top', fontsize = 12, color = 'r', zorder = 99)
  
      # Appearance
      ax.margins(0, None)
      ax.set_xlabel("Time (%s)" % self._mission.TIMEUNITS, fontsize = 16)
      ax.set_ylabel(ylabel, fontsize = 16)
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
    
    kwargs.update({'season', self.season})
    return getattr(missions, self.mission).pipelines.get(self.ID, *args, **kwargs)
  
  def mask_planet(self, t0, period, dur = 0.2):
    '''
    Mask all of the transits/eclipses of a given planet/EB. After calling
    this method, you must re-compute the model by calling :py:meth:`compute`
    in order for the mask to take effect.
    
    :param float t0: The time of first transit (same units as light curve)
    :param float period: The period of the planet in days
    :param foat dur: The transit duration in days. Default 0.2
    
    '''
    
    for k in range(self.nsub):      
      if self.nsub == 1:
        k = slice(None, None, None)
      mask = []
      t0 += np.ceil((self.time[k][0] - dur - t0) / period) * period
      for t in np.arange(t0, self.time[k][-1] + dur, period):
        mask.extend(np.where(np.abs(self.time[k] - t) < dur / 2.)[0])
      self.transitmask[k] = np.array(list(set(np.concatenate([self.transitmask[k], mask]))))
    
  def _plot_weights(self, show = True):
    '''
    .. warning:: Untested!
    
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
    
    # Check number of segments
    if len(self.breakpoints) > 3:
      log.error('Cannot currently plot weights for light curves with more than 3 segments.')
      return
    
    # Loop over all PLD orders and over all chunks
    npix = len(self.fpix[1])
    ap = self.aperture.flatten()
    ncol = 1 + 2 * (len(self.weights[0]) - 1)
    raw_weights = np.zeros((len(self.breakpoints), ncol, self.aperture.shape[0], self.aperture.shape[1]), dtype = float)
    scaled_weights = np.zeros((len(self.breakpoints), ncol, self.aperture.shape[0], self.aperture.shape[1]), dtype = float)
    
    # Loop over orders
    for o in range(len(self.weights[0])):
      if o == 0:
        oi = 0
      else:
        oi = 1 + 2 * (o - 1)
        
      # Loop over chunks
      for b in range(len(self.weights)):
      
        c = self.get_chunk(b)
        rw_ii = np.zeros(npix); rw_ij = np.zeros(npix)
        sw_ii = np.zeros(npix); sw_ij = np.zeros(npix)
        X = np.nanmedian(self.X(o,c), axis = 0)
      
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
        rw = np.zeros_like(ap, dtype = float)
        sw = np.zeros_like(ap, dtype = float)
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
          rw = np.zeros_like(ap, dtype = float)
          sw = np.zeros_like(ap, dtype = float)
          n = 0
          for i, a in enumerate(ap):
            if (a & 1):
              rw[i] = rw_ij[n]
              sw[i] = sw_ij[n]
              n += 1
          raw_weights[b][oi + 1] = rw.reshape(*self.aperture.shape)
          scaled_weights[b][oi + 1] = sw.reshape(*self.aperture.shape)

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
        imr = ax[2 * b, o].imshow(raw_weights[b][o], aspect = 'auto', interpolation = 'nearest', cmap = rdbu, origin = 'lower', vmin = -rmax, vmax = rmax)
        ims = ax[2 * b + 1, o].imshow(scaled_weights[b][o], aspect = 'auto', interpolation = 'nearest', cmap = rdbu, origin = 'lower', vmin=-smax, vmax=smax)
    
      # Colorbars
      def fmt(x, pos):
        a, b = '{:.0e}'.format(x).split('e')
        b = int(b)
        if float(a) > 0:
          a = r'+' + a
        elif float(a) == 0:
          return ''
        return r'${} \times 10^{{{}}}$'.format(a, b) 
      cbr = pl.colorbar(imr, cax = cax[2 * b], format = FuncFormatter(fmt))
      cbr.ax.tick_params(labelsize = 8) 
      cbs = pl.colorbar(ims, cax = cax[2 * b + 1], format = FuncFormatter(fmt))
      cbs.ax.tick_params(labelsize = 8) 
  
    # Plot aperture contours
    def PadWithZeros(vector, pad_width, iaxis, kwargs):
      vector[:pad_width[0]] = 0
      vector[-pad_width[1]:] = 0
      return vector
    ny, nx = self.aperture.shape
    contour = np.zeros((ny,nx))
    contour[np.where(self.aperture)] = 1
    contour = np.lib.pad(contour, 1, PadWithZeros)
    highres = zoom(contour, 100, order = 0, mode='nearest') 
    extent = np.array([-1, nx, -1, ny])
    for axis in ax.flatten():
      axis.contour(highres, levels=[0.5], extent=extent, origin='lower', colors='r', linewidths=1)
      
      # Check for saturated columns
      for x in range(self.aperture.shape[0]):
        for y in range(self.aperture.shape[1]):
          if self.aperture[x][y] == AP_SATURATED_PIXEL:
            axis.fill([y - 0.5, y + 0.5, y + 0.5, y - 0.5], 
                      [x - 0.5, x - 0.5, x + 0.5, x + 0.5], fill = False, hatch='xxxxx', color = 'r', lw = 0)
      
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
      axis.set_title(titles[i], fontsize = 12)
    for j in range(len(self.weights)):
      ax[2 * j, 0].text(-0.55, -0.15, r'$%d$' % (j + 1), fontsize = 16, transform = ax[2 * j, 0].transAxes)
      ax[2 * j, 0].set_ylabel(r'$w_{ij}$', fontsize = 18)
      ax[2 * j + 1, 0].set_ylabel(r'$\bar{X}_{ij} \cdot w_{ij}$', fontsize = 18)
    
    if show:
      pl.show()
      pl.close()
    else:
      return fig, ax, cax
      
  def _plot_chunks(self, show = True, plot_bad = True, plot_out = True):
    '''
  
    '''

    log.info('Plotting the light curve...')
  
    # Set up axes    
    fig, ax = pl.subplots(len(self.breakpoints), figsize = (10, 8))
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
      mK = GetCovariance(self.kernel_params, self.time[m], self.fraw_err[m])
      med = np.nanmedian(self.fraw[m])
      f = self.fraw[m] - med
      A = np.zeros((len(m), len(m)))
      B = np.zeros((len(c), len(m)))
      for n in range(self.pld_order):
        if (self.lam_idx >= n) and (self.lam[b][n] is not None):
          XM = self.X(n,m)
          XC = self.X(n,c)
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
      ax[i].plot(self.apply_mask(self.time), self.apply_mask(flux), ls = 'none', marker = '.', color = 'k', markersize = ms, alpha = 0.5)
  
      # Plot the outliers
      bnmask = np.array(list(set(np.concatenate([self.badmask, self.nanmask]))), dtype = int)
      O1 = lambda x: x[self.outmask]
      O2 = lambda x: x[bnmask]
      O3 = lambda x: x[self.transitmask]
      if plot_out:
        ax[i].plot(O1(self.time), O1(flux), ls = 'none', color = "#777777", marker = '.', markersize = ms, alpha = 0.5)
      if plot_bad:
        ax[i].plot(O2(self.time), O2(flux), 'r.', markersize = ms, alpha = 0.25)
      ax[i].plot(O3(self.time), O3(flux), 'b.', markersize = ms, alpha = 0.25)
      
      # Appearance
      if i == len(self.breakpoints) - 1: 
        ax[i].set_xlabel('Time (%s)' % self._mission.TIMEUNITS, fontsize = 18)
      ax[i].set_ylabel('Flux %d' % (i + 1), fontsize = 18)
      for brkpt in self.breakpoints[:-1]:
        ax[i].axvline(self.time[brkpt], color = 'r', ls = '--', alpha = 0.25)
      if len(self.breakpoints) == 2:
        ax[i].annotate('%.2f ppm' % cdpp_arr[0], xy = (0.02, 0.975), xycoords = 'axes fraction', 
                       ha = 'left', va = 'top', fontsize = 12, color = 'r', zorder = 99)
        ax[i].annotate('%.2f ppm' % cdpp_arr[1], xy = (0.98, 0.975), xycoords = 'axes fraction', 
                       ha = 'right', va = 'top', fontsize = 12, color = 'r', zorder = 99)
      elif len(self.breakpoints) < 6:
        for n in range(len(self.breakpoints)):
          if n > 0:
            x = (self.time[self.breakpoints[n - 1]] - self.time[0]) / (self.time[-1] - self.time[0]) + 0.02
          else:
            x = 0.02
          ax[i].annotate('%.2f ppm' % cdpp_arr[n], xy = (x, 0.975), xycoords = 'axes fraction', 
                         ha = 'left', va = 'top', fontsize = 10, zorder = 99, color = 'r')
      else:
        ax[i].annotate('%.2f ppm' % cdpp_arr[0], xy = (0.02, 0.975), xycoords = 'axes fraction', 
                       ha = 'left', va = 'top', fontsize = 12, color = 'r', zorder = 99)
      ax[i].margins(0.01, 0.1)
      if i == 0:
        a = self.time[0]
      else:
        a = self.time[self.breakpoints[i-1]]
      b = self.time[self.breakpoints[i]]
      ax[i].axvspan(a, b, color = 'b', alpha = 0.1, zorder = -99)
      
      # Get y lims that bound 99% of the flux
      f = np.concatenate([np.delete(f, bnmask) for f in fluxes])
      N = int(0.995 * len(f))
      hi, lo = f[np.argsort(f)][[N,-N]]
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
        ax[i].annotate('', xy=(self.time[j], ylim[0]), xycoords = 'data',
                       xytext = (0, 15), textcoords = 'offset points',
                       arrowprops=dict(arrowstyle = "-|>", color = color))
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
        ax[i].annotate('', xy=(self.time[j], ylim[1]), xycoords = 'data',
                       xytext = (0, -15), textcoords = 'offset points',
                       arrowprops=dict(arrowstyle = "-|>", color = color))

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
  
  def optimize(self, piter = 3, pmaxf = 300, ppert = 0.1):
    '''
    Runs :py:obj:`pPLD` on the target in an attempt to further optimize the
    values of the PLD priors. See :py:class:`everest.detrender.pPLD`.
    
    '''
    
    self._save_npz()
    optimized = pPLD(self.ID, piter = piter, pmaxf = pmaxf, ppert = ppert, debug = True, clobber = True)
    optimized.publish()
    self.reset()
  
  def plot_folded(self, t0, period, dur = 0.2):
    '''
    Plot the light curve folded on a given `period` and centered at `t0`. 
    When plotting folded transits, please mask them using :py:meth:`mask_planet`
    and re-compute the model using :py:meth:`compute`.
    
    :param float t0: The time at which to center the plot (same units as light curve)
    :param float period: The period of the folding operation
    :param float dur: The transit duration in days. Default 0.2
    
    '''
    
    # Mask the planet
    self.mask_planet(t0, period, dur)
    
    # Whiten
    if self.nsub == 1:
      _, amp, tau = self.kernel_params
      gp = george.GP(amp ** 2 * george.kernels.Matern32Kernel(tau ** 2))
      gp.compute(self.apply_mask(self.time), self.apply_mask(self.fraw_err))
      med = np.nanmedian(self.apply_mask(self.flux))
      y, _ = gp.predict(self.apply_mask(self.flux) - med, self.time)
      fwhite = (self.flux - y)
      fwhite /= np.nanmedian(fwhite)
      time = np.array(self.time)
    else:
      fwhite = np.array([])
      _, amp, tau = self.kernel_params
      gp = george.GP(amp ** 2 * george.kernels.Matern32Kernel(tau ** 2))
      for k in range(self.nsub):
        gp.compute(self.apply_mask(self.time[k], k = k), self.apply_mask(self.fraw_err[k], k = k))
        med = np.nanmedian(self.apply_mask(self.flux[k], k = k))
        y, _ = gp.predict(self.apply_mask(self.flux[k], k = k) - med, self.time[k])
        fw = (self.flux[k] - y)
        fwhite = np.append(fwhite, fw / np.nanmedian(fw))
      time = np.concatenate(self.time)
    
    # Fold
    tfold = (time - t0 - period / 2.) % period - period / 2. 
    
    # Crop
    inds = np.where(np.abs(tfold) < 2 * dur)[0]
    x = tfold[inds]
    y = fwhite[inds]
    
    # Plot
    fig, ax = pl.subplots(1, figsize = (9, 5))
    fig.subplots_adjust(bottom = 0.125)
    ax.plot(x, y, 'k.', alpha = 0.5)
    
    # Get ylims
    yfin = np.delete(y, np.where(np.isnan(y)))
    lo, hi = yfin[np.argsort(yfin)][[3,-3]]
    pad = (hi - lo) * 0.1
    ylim = (lo - pad, hi + pad)
    ax.set_ylim(*ylim)
    
    # Appearance
    ax.set_xlabel(r'Time (days)', fontsize = 18)
    ax.set_ylabel(r'Normalized Flux', fontsize = 18)
    fig.canvas.set_window_title('%s %d' % (self._mission.IDSTRING, self.ID))
    
    pl.show()