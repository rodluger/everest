#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`k2.py` - Main mission routines
---------------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from ... import __version__ as EVEREST_VERSION
from . import cbv
from .aux import *
from ...config import EVEREST_SRC, EVEREST_DAT, EVEREST_DEV, MAST_ROOT, MAST_VERSION
from ...utils import DataContainer, sort_like, AP_COLLAPSED_PIXEL, AP_SATURATED_PIXEL
from ...math import SavGol, Interpolate, Scatter
try:
  import pyfits
except ImportError:
  try:
    import astropy.io.fits as pyfits
  except ImportError:
    raise Exception('Please install the `pyfits` package.')
import matplotlib.pyplot as pl
from matplotlib.ticker import ScalarFormatter, MaxNLocator
import k2plr as kplr; kplr_client = kplr.API()
from k2plr.config import KPLR_ROOT
import numpy as np
from tempfile import NamedTemporaryFile
import random
import os, sys, shutil
import logging
log = logging.getLogger(__name__)

__all__ = ['Setup', 'Season', 'Breakpoints', 'GetData', 'GetNeighbors', 
           'Statistics', 'TargetDirectory', 'HasShortCadence', 
           'InjectionStatistics', 'HDUCards', 'FITSFile', 'FITSUrl', 'CDPP']

def Setup():
  '''
  
  '''
  
  if not os.path.exists(os.path.join(EVEREST_DAT, 'k2', 'stats')):
    os.makedirs(os.path.join(EVEREST_DAT, 'k2', 'stats'))
  if not os.path.exists(os.path.join(EVEREST_DAT, 'k2', 'cbv')):
    os.makedirs(os.path.join(EVEREST_DAT, 'k2', 'cbv'))
  GetK2Stars(clobber = False)

def Season(EPIC, **kwargs):
  '''
  Returns the campaign number for a given EPIC target.
  
  '''
  
  return Campaign(EPIC, **kwargs)

def Breakpoints(EPIC, cadence = 'lc', **kwargs):  
  '''
  ..note :: The number corresponding to a given breakpoint is the number \
            of cadences *since the beginning of the campaign*.
  
  '''
  
  # Get the campaign number
  campaign = Season(EPIC)
  
  # Select LC or SC
  if cadence == 'lc':
    breakpoints = {
                   0: [665],          # OK
                   1: [2210],         # OK
                   2: [2042],         # OK
                   3: [2140],         # OK
                   4: [520, 2153],    # OK
                   5: [1774],         # OK
                   6: [2143],         # OK
                   7: [1192, 2319],   # OK
                   8: [1950],         # OK
                   9: [],
                  10: [],
                  11: [],
                  12: [],
                  13: [],
                  14: [],
                  15: [],
                  16: [],
                  17: []
                  }
  elif cadence == 'sc':
    breakpoints = {
                  0: np.array(np.linspace(0, 112590, 31)[1:-1], dtype = int),     # OK
                  1: np.array([ 8044,   12066,  16088,  20135,  24132,  28154,    # OK
                                32176,  36198,  40220,  44242,  48264,  52286,  
                                56308,  60330,  64352,  68374,  72396,  76418,  
                                80440,  84462,  88509,  92506,  96528,  100550,
                                104572, 108594, 112616, 116638]),
                  2: np.array(np.linspace(0, 115680, 31)[1:-1], dtype = int),     # OK
                  3: np.array([  3316,  6772,   10158,  13694,  16930,  20316,    # OK
                                23702,  27088,  30474,  33860,  37246,  40632,  
                                44018,  47404,  50790,  54176,  57562,  60948,  
                                64334,  67720,  71106,  74492,  77878,  81264,
                                84650,  88036,  91422,  94808,  98194]),
                  4: np.array(np.linspace(0, 101580, 31)[1:-1], dtype = int),     # OK
                  5: np.array([  3663,   7326,  10989,  14652,  18315,  21978,    # OK
                                25641,  29304,  32967,  36630,  40293,  43956,  
                                47619,  51282,  54945,  58608,  62271,  65934,  
                                69597,  73260,  76923,  80646,  84249,  87912,
                                91575,  95238,  98901, 102564, 106227]),
                  6: np.array(np.linspace(0, 115890, 31)[1:-1], dtype = int),     # OK
                  7: np.array(np.linspace(0, 121290, 31)[1:-1], dtype = int),     # OK
                  8: np.array(np.linspace(0, 115590, 31)[1:-1], dtype = int),     # Unclear
                  9: [],
                 10: [],
                 11: [],
                 12: [],
                 13: [],
                 14: [],
                 15: [],
                 16: [],
                 17: []
                  }
  else:
    raise ValueError("Invalid value for the cadence.")
  
  # Return
  if campaign in breakpoints:
    return breakpoints[campaign]
  else:
    return None

def CDPP(flux, mask = [], cadence = 'lc'):
  '''
  Compute the proxy 6-hr CDPP metric.
  
  '''
  
  if cadence == 'lc':
    rmswin = 13
    svgwin = 49
  else:
    rmswin = 13 * 30
    svgwin = 50 * 30 - 1
  flux_savgol = SavGol(np.delete(flux, mask), win = svgwin)
  if len(flux_savgol):
    return Scatter(flux_savgol / np.nanmedian(flux_savgol), remove_outliers = True, win = rmswin) 
  else:
    return np.nan
    
def GetData(EPIC, season = None, cadence = 'lc', clobber = False, delete_raw = False, 
            aperture_name = 'k2sff_15', saturated_aperture_name = 'k2sff_19',
            max_pixels = 75, download_only = False, saturation_tolerance = -0.1, 
            bad_bits = [1,2,3,4,5,6,7,8,9,11,12,13,14,16,17], **kwargs):
  '''

  '''
  
  # Campaign no.
  if season is None:
    campaign = Season(EPIC)
  else:
    campaign = season

  # Is there short cadence data available for this target?
  short_cadence = HasShortCadence(EPIC, season = campaign)
  if cadence == 'sc' and not short_cadence:
    raise ValueError("Short cadence data not available for this target.")
    
  # Local file name
  filename = os.path.join(EVEREST_DAT, 'k2', 'c%02d' % campaign, 
                         ('%09d' % EPIC)[:4] + '00000', ('%09d' % EPIC)[4:],
                         'data.npz')
  
  # Download?
  if clobber or not os.path.exists(filename):

    # Get the TPF
    tpf = os.path.join(KPLR_ROOT, 'data', 'k2', 'target_pixel_files', 
                       str(EPIC), 'ktwo%09d-c%02d_lpd-targ.fits.gz' % (EPIC, campaign))
    sc_tpf = os.path.join(KPLR_ROOT, 'data', 'k2', 'target_pixel_files', 
                          str(EPIC), 'ktwo%09d-c%02d_spd-targ.fits.gz' % (EPIC, campaign))
    if clobber or not os.path.exists(tpf):                 
      kplr_client.k2_star(EPIC).get_target_pixel_files(fetch = True)

    with pyfits.open(tpf) as f:
      qdata = f[1].data
      
      # Get the TPF aperture
      tpf_aperture = (f[2].data & 2) // 2
    
      # Get the enlarged TPF aperture
      tpf_big_aperture = np.array(tpf_aperture)
      for i in range(tpf_big_aperture.shape[0]):
        for j in range(tpf_big_aperture.shape[1]):
          if f[2].data[i][j] == 1:
            for n in [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]:
              if n[0] >= 0 and n[0] < tpf_big_aperture.shape[0]:
                if n[1] >= 0 and n[1] < tpf_big_aperture.shape[1]:
                  if tpf_aperture[n[0]][n[1]] == 1:
                    tpf_big_aperture[i][j] = 1
    
    # Is there short cadence data?
    if short_cadence:
      with pyfits.open(sc_tpf) as f:
        sc_qdata = f[1].data
    
    # Get K2SFF apertures
    try:
      k2sff = kplr.K2SFF(EPIC)
      k2sff_apertures = k2sff.apertures
      if delete_raw:
        os.remove(k2sff._file)
    except:
      k2sff_apertures = [None for i in range(20)]
    
    # Make a dict of all our apertures
    # We're not getting K2SFF apertures 0-9 any more
    apertures = {'tpf': tpf_aperture, 'tpf_big': tpf_big_aperture}
    for i in range(10,20):
      apertures.update({'k2sff_%02d' % i: k2sff_apertures[i]})
    
    # Get the header info
    fitsheader = [pyfits.getheader(tpf, 0).cards,
                  pyfits.getheader(tpf, 1).cards,
                  pyfits.getheader(tpf, 2).cards]
    if short_cadence:
      sc_fitsheader = [pyfits.getheader(sc_tpf, 0).cards,
                       pyfits.getheader(sc_tpf, 1).cards,
                       pyfits.getheader(sc_tpf, 2).cards]
    else:
      sc_fitsheader = None
      
    # Get a hi res image of the target
    hires = GetHiResImage(EPIC)
    
    # Get nearby sources
    nearby = GetSources(EPIC)
    
    # Delete?
    if delete_raw:
      os.remove(tpf)
      if short_cadence:
        os.remove(sc_tpf)
  
    # Get the arrays
    cadn = np.array(qdata.field('CADENCENO'), dtype='int32')
    time = np.array(qdata.field('TIME'), dtype='float64')
    fpix = np.array(qdata.field('FLUX'), dtype='float64')
    fpix_err = np.array(qdata.field('FLUX_ERR'), dtype='float64')
    qual = np.array(qdata.field('QUALITY'), dtype=int)
    
    # Get rid of NaNs in the time array by interpolating
    naninds = np.where(np.isnan(time))
    time = Interpolate(np.arange(0, len(time)), naninds, time)
    
    # Get the motion vectors (if available!)
    pc1 = np.array(qdata.field('POS_CORR1'), dtype='float64')
    pc2 = np.array(qdata.field('POS_CORR2'), dtype='float64')
    if not np.all(np.isnan(pc1)) and not np.all(np.isnan(pc2)):
      pc1 = Interpolate(time, np.where(np.isnan(pc1)), pc1)
      pc2 = Interpolate(time, np.where(np.isnan(pc2)), pc2)
    else:
      pc1 = None
      pc2 = None
    
    # Do the same for short cadence
    if short_cadence:
      sc_cadn = np.array(sc_qdata.field('CADENCENO'), dtype='int32')
      sc_time = np.array(sc_qdata.field('TIME'), dtype='float64')
      sc_fpix = np.array(sc_qdata.field('FLUX'), dtype='float64')
      sc_fpix_err = np.array(sc_qdata.field('FLUX_ERR'), dtype='float64')
      sc_qual = np.array(sc_qdata.field('QUALITY'), dtype=int)
      sc_naninds = np.where(np.isnan(sc_time))
      sc_time = Interpolate(np.arange(0, len(sc_time)), sc_naninds, sc_time)
      sc_pc1 = np.array(sc_qdata.field('POS_CORR1'), dtype='float64')
      sc_pc2 = np.array(sc_qdata.field('POS_CORR2'), dtype='float64')
      if not np.all(np.isnan(sc_pc1)) and not np.all(np.isnan(sc_pc2)):
        sc_pc1 = Interpolate(sc_time, np.where(np.isnan(sc_pc1)), sc_pc1)
        sc_pc2 = Interpolate(sc_time, np.where(np.isnan(sc_pc2)), sc_pc2)
      else:
        sc_pc1 = None
        sc_pc2 = None
    else:
      sc_cadn = None
      sc_time = None
      sc_fpix = None
      sc_fpix_err = None
      sc_qual = None
      sc_pc1 = None
      sc_pc2 = None
    
    # Static pixel images for plotting
    pixel_images = [fpix[0], fpix[len(fpix) // 2], fpix[len(fpix) - 1]]
    
    # Atomically write to disk.
    # http://stackoverflow.com/questions/2333872/atomic-writing-to-file-with-python
    if not os.path.exists(os.path.dirname(filename)):
      os.makedirs(os.path.dirname(filename))                   
    f = NamedTemporaryFile("wb", delete = False)
    np.savez_compressed(f, cadn = cadn, time = time, fpix = fpix, fpix_err = fpix_err, 
                        qual = qual, apertures = apertures,  
                        pc1 = pc1, pc2 = pc2, fitsheader = fitsheader,
                        pixel_images = pixel_images, nearby = nearby, hires = hires,
                        sc_cadn = sc_cadn, sc_time = sc_time, sc_fpix = sc_fpix,
                        sc_fpix_err = sc_fpix_err, sc_qual = sc_qual,
                        sc_pc1 = sc_pc1, sc_pc2 = sc_pc2, sc_fitsheader = sc_fitsheader)
    f.flush()
    os.fsync(f.fileno())
    f.close()
    shutil.move(f.name, filename)
    
    if download_only:
      return
  
  # Load
  data = np.load(filename)
  apertures = data['apertures'][()]
  pixel_images = data['pixel_images']
  nearby = data['nearby']
  hires = data['hires']
  
  if cadence == 'lc':
    fitsheader = data['fitsheader']
    cadn = data['cadn']
    time = data['time']
    fpix = data['fpix']
    fpix_err = data['fpix_err']
    qual = data['qual']
    pc1 = data['pc1']
    pc2 = data['pc2']
  elif cadence == 'sc':
    fitsheader = data['sc_fitsheader']
    cadn = data['sc_cadn']
    time = data['sc_time']
    fpix = data['sc_fpix']
    fpix_err = data['sc_fpix_err']
    qual = data['sc_qual']
    pc1 = data['sc_pc1']
    pc2 = data['sc_pc2']
  else:
    raise ValueError("Invalid value for the cadence.")
    
  # Select the "saturated aperture" to check if the star is saturated
  # If it is, we will use this aperture instead
  if saturated_aperture_name == 'custom':
    saturated_aperture = GetCustomAperture(data)
  else:
    if saturated_aperture_name is None:
      saturated_aperture_name = 'k2sff_19'
    saturated_aperture = apertures[saturated_aperture_name]
    if saturated_aperture is None:
      log.error("Invalid aperture selected. Defaulting to `tpf_big`.")
      saturated_aperture_name = 'tpf_big'
      saturated_aperture = apertures[saturated_aperture_name]
    
  # Compute the saturation flux and the 97.5th percentile 
  # flux in each pixel of the saturated aperture. We're going
  # to compare these to decide if the star is saturated.
  satflx = SaturationFlux(EPIC) * (1. + saturation_tolerance)
  f97 = np.zeros((fpix.shape[1], fpix.shape[2]))
  for i in range(fpix.shape[1]):
    for j in range(fpix.shape[2]):
      if saturated_aperture[i,j]:
        # Let's remove NaNs...
        tmp = np.delete(fpix[:,i,j], np.where(np.isnan(fpix[:,i,j])))
        # ... and really bad outliers...
        if len(tmp):
          f = SavGol(tmp)
          med = np.nanmedian(f)
          MAD = 1.4826 * np.nanmedian(np.abs(f - med))
          bad = np.where((f > med + 10. * MAD) | (f < med - 10. * MAD))[0]
          np.delete(tmp, bad)
          # ... so we can compute the 97.5th percentile flux
          i97 = int(0.975 * len(tmp))
          tmp = tmp[np.argsort(tmp)[i97]]
          f97[i,j] = tmp
  
  # Check if any of the pixels are actually saturated
  if np.nanmax(f97) <= satflx:
    log.info("No saturated columns detected.")
    saturated = False
  else:
    log.info("Saturated pixel(s) found. Switching to aperture `%s`." % saturated_aperture_name)
    aperture_name = saturated_aperture_name
    saturated = True
  
  # Now grab the aperture we'll actually use
  if aperture_name == 'custom':
    aperture = GetCustomAperture(data)
  else:
    if aperture_name is None:
      aperture_name = 'k2sff_15'
    aperture = apertures[aperture_name]
    if aperture is None:
      log.error("Invalid aperture selected. Defaulting to `tpf_big`.")
      aperture_name = 'tpf_big'
      aperture = apertures[aperture_name]

  # Now we check if the aperture is too big. Can lead to memory errors...
  # Treat saturated and unsaturated stars differently.
  if saturated:
    
    # We need to check if we have too many pixels *after* collapsing the columns.
    # Sort the apertures in decreasing order of pixels, but keep the aperture
    # chosen by the user first.
    aperture_names = np.array(list(apertures.keys()))
    npix_per_aperture = np.array([np.sum(apertures[k]) for k in aperture_names])
    aperture_names = aperture_names[np.argsort(npix_per_aperture)[::-1]]
    aperture_names = np.append([aperture_name], np.delete(aperture_names, np.argmax(aperture_names == aperture_name)))
    
    # Loop through them. Pick the first one that satisfies the `max_pixels` constraint
    for aperture_name in aperture_names:        
      aperture = apertures[aperture_name]
      aperture[np.isnan(fpix[0])] = 0
      ncol = 0
      apcopy = np.array(aperture)
      for j in range(apcopy.shape[1]):
        if np.any(f97[:,j] > satflx):
          apcopy[:,j] = 0
          ncol += 1
      if np.sum(apcopy) + ncol <= max_pixels:
        break
    assert np.sum(apcopy) + ncol <= max_pixels, "No apertures available with fewer than %d pixels. Aborting." % max_pixels
  
    # Now, finally, we collapse the saturated columns into single pixels
    # and make the pixel array 2D
    ncol = 0
    fpixnew = []
    ferrnew = []
    
    # HACK: K2SFF sometimes clips the heads/tails of saturated columns
    # That's really bad, since that's where all the information is. Let's
    # artificially extend the aperture by two pixels at the top and bottom
    # of each saturated column. This *could* increase contamination, but
    # it's unlikely since the saturated target is by definition really bright
    ext = 0
    for j in range(aperture.shape[1]):
      if np.any(f97[:,j] > satflx):
        for i in range(aperture.shape[0]):
          if (aperture[i,j] == 0) and (np.nanmedian(fpix[:,i,j]) > 0):
            if (i+2 < aperture.shape[0]) and aperture[i+2,j] == 1:
              aperture[i,j] = 2
              ext += 1
            elif (i+1 < aperture.shape[0]) and aperture[i+1,j] == 1:
              aperture[i,j] = 2
              ext += 1
            elif (i-1 >= 0) and aperture[i-1,j] == 1:
              aperture[i,j] = 2
              ext += 1
            elif (i-2 >= 0) and aperture[i-2,j] == 1:
              aperture[i,j] = 2
              ext += 1
    if ext: 
      log.info("Extended saturated columns by %d pixel(s)." % ext)        
    
    for j in range(aperture.shape[1]):
      if np.any(f97[:,j] > satflx):
        marked = False
        collapsed = np.zeros(len(fpix[:,0,0]))
        collapsed_err2 = np.zeros(len(fpix[:,0,0]))
        for i in range(aperture.shape[0]):
          if aperture[i,j]:
            if not marked:
              aperture[i,j] = AP_COLLAPSED_PIXEL
              marked = True
            else:
              aperture[i,j] = AP_SATURATED_PIXEL
            collapsed += fpix[:,i,j]
            collapsed_err2 += fpix_err[:,i,j] ** 2
        if np.any(collapsed):
          fpixnew.append(collapsed)
          ferrnew.append(np.sqrt(collapsed_err2))
          ncol += 1
      else:
        for i in range(aperture.shape[0]):
          if aperture[i,j]:
            fpixnew.append(fpix[:,i,j])
            ferrnew.append(fpix_err[:,i,j])
    fpix2D = np.array(fpixnew).T
    fpix_err2D = np.array(ferrnew).T
    log.info("Collapsed %d saturated column(s)." % ncol)

  else:
    
    # Check if there are too many pixels
    if np.sum(aperture) > max_pixels:
  
      # This case is simpler: we just pick the largest aperture that's less than or equal to `max_pixels`
      keys = list(apertures.keys())
      npix = np.array([np.sum(apertures[k]) for k in keys])
      aperture_name = keys[np.argmax(npix * (npix <= max_pixels))]
      aperture = apertures[aperture_name]
      aperture[np.isnan(fpix[0])] = 0
      assert np.sum(aperture) <= max_pixels, "No apertures available with fewer than %d pixels. Aborting." % max_pixels
      log.warn("Selected aperture is too big. Proceeding with aperture `%s` instead." % aperture_name)
    
    # Make the pixel flux array 2D
    aperture[np.isnan(fpix[0])] = 0
    ap = np.where(aperture & 1)
    fpix2D = np.array([f[ap] for f in fpix], dtype='float64')
    fpix_err2D = np.array([p[ap] for p in fpix_err], dtype='float64')
    
  # Compute the background
  binds = np.where(aperture ^ 1)
  if RemoveBackground(EPIC) and (len(binds[0]) > 0):
    bkg = np.nanmedian(np.array([f[binds] for f in fpix], dtype='float64'), axis = 1)
    # Uncertainty of the median: http://davidmlane.com/hyperstat/A106993.html
    bkg_err = 1.253 * np.nanmedian(np.array([e[binds] for e in fpix_err], 
                      dtype='float64'), axis = 1) / np.sqrt(len(binds[0]))
    bkg = bkg.reshape(-1, 1)
    bkg_err = bkg_err.reshape(-1, 1)
  else:
    bkg = 0.
    bkg_err = 0.
  
  # Make everything 2D and remove the background
  fpix = fpix2D - bkg
  fpix_err = np.sqrt(fpix_err2D ** 2 + bkg_err ** 2)
  flux = np.sum(fpix, axis = 1)
  ferr = np.sqrt(np.sum(fpix_err ** 2, axis = 1))
    
  # Get NaN data points
  nanmask = np.where(np.isnan(flux) | (flux == 0))[0]
  
  # Get flagged data points -- we won't train our model on them                         
  badmask = []
  for b in bad_bits:
    badmask += list(np.where(qual & 2 ** (b - 1))[0])
  
  # Flag >10 sigma outliers -- same thing.
  tmpmask = np.array(list(set(np.concatenate([badmask, nanmask]))))
  t = np.delete(time, tmpmask)
  f = np.delete(flux, tmpmask)
  f = SavGol(f)
  med = np.nanmedian(f)
  MAD = 1.4826 * np.nanmedian(np.abs(f - med))
  bad = np.where((f > med + 10. * MAD) | (f < med - 10. * MAD))[0]
  badmask.extend([np.argmax(time == t[i]) for i in bad])
  
  # Campaign 2 hack: the first day or two are screwed up
  if campaign == 2:
    badmask.extend(np.where(time < 2061.5)[0])
  
  # TODO: Fix time offsets in first half of 
  # Campaign 0. See note in everest 1.0 code
  
  # Finalize the mask
  badmask = np.array(sorted(list(set(badmask))))
  
  # Interpolate the nans
  fpix = Interpolate(time, nanmask, fpix)
  fpix_err = Interpolate(time, nanmask, fpix_err)

  # Return
  data = DataContainer()
  data.ID = EPIC
  data.campaign = campaign
  data.cadn = cadn
  data.time = time
  data.fpix = fpix
  data.fpix_err = fpix_err
  data.nanmask = nanmask
  data.badmask = badmask
  data.aperture = aperture
  data.aperture_name = aperture_name
  data.apertures = apertures
  data.quality = qual
  data.Xpos = pc1
  data.Ypos = pc2
  data.meta = fitsheader
  data.mag = fitsheader[0]['KEPMAG'][1]
  data.pixel_images = pixel_images
  data.nearby = nearby
  data.hires = hires
  data.saturated = saturated
  data.bkg = bkg
    
  return data

def GetNeighbors(EPIC, model = None, neighbors = 10, mag_range = (11., 13.), 
                 cdpp_range = None, aperture_name = 'k2sff_15', 
                 cadence = 'lc', **kwargs):
  '''
  Return `neighbors` random bright stars on the same module as `EPIC`.
  
  '''
  
  # Zero neighbors?
  if neighbors == 0:
    return []
  
  # Get the IDs
  campaign = Season(EPIC)
  epics, kepmags, channels, short_cadence = np.array(GetK2Stars()[campaign]).T
  short_cadence = np.array(short_cadence, dtype = bool)
  epics = np.array(epics, dtype = int)
  c = GetNeighboringChannels(Channel(EPIC))
  
  # Manage kwargs
  if aperture_name is None:
    aperture_name = 'k2sff_15'
  if mag_range is None:
    mag_lo = -np.inf
    mag_hi = np.inf
  else:
    mag_lo = mag_range[0]
    mag_hi = mag_range[1]
    # K2-specific tweak. The short cadence stars are preferentially really bright ones,
    # so we won't get many neighbors if we stick to the default magnitude range! I'm 
    # therefore enforcing a lower magnitude cut-off of 8.
    if cadence == 'sc':
      mag_lo = 8.
  if cdpp_range is None:
    cdpp_lo = -np.inf
    cdpp_hi = np.inf
  else:
    cdpp_lo = cdpp_range[0]
    cdpp_hi = cdpp_range[1]
  targets = []
  
  # First look for nearby targets, then relax the constraint
  # If still no targets, widen magnitude range
  for n in range(3):
    
    if n == 0:
      nearby = True
    elif n == 1:
      nearby = False
    elif n == 2:
      mag_lo -= 1
      mag_hi += 1
    
    # Loop over all stars
    for star, kp, channel, sc in zip(epics, kepmags, channels, short_cadence):
    
      # Preliminary vetting
      if not (((channel in c) if nearby else True) and (kp < mag_hi) and (kp > mag_lo) and (sc if cadence == 'sc' else True)):
        continue
      
      # Reject if self or if already in list
      if (star == EPIC) or (star in targets):
        continue
    
      # Ensure raw light curve file exists
      if not os.path.exists(os.path.join(TargetDirectory(star, campaign), 'data.npz')):
        continue
      
      # Ensure crowding is OK. This is quite conservative, as we
      # need to prevent potential astrophysical false positive contamination
      # from crowded planet-hosting neighbors when doing neighboring PLD.
      contam = False
      data = np.load(os.path.join(TargetDirectory(star, campaign), 'data.npz'))
      aperture = data['apertures'][()][aperture_name]
      for source in data['nearby'][()]:
        # Ignore self
        if source['ID'] == star:
          continue
        # Ignore really dim stars
        if source['mag'] < kp - 5:
          continue
        # Compute source position
        x = int(np.round(source['x'] - source['x0']))
        y = int(np.round(source['y'] - source['y0']))
        # If the source is within two pixels of the edge
        # of the target aperture, reject the target
        for j in [x - 2, x - 1, x, x + 1, x + 2]:
          if j < 0:
            # Outside the postage stamp
            continue
          for i in [y - 2, y - 1, y, y + 1, y + 2]:
            if i < 0:
              # Outside the postage stamp
              continue
            try:
              if aperture[i][j]:
                # Oh-oh!
                contam = True
            except IndexError:
              # Out of bounds... carry on!
              pass
      if contam:
        continue
    
      # Reject if the model is not present
      if model is not None:
        if not os.path.exists(os.path.join(TargetDirectory(star, campaign), model + '.npz')):
          continue
        
        # Reject if CDPP out of range
        if cdpp_range is not None:
          cdpp = np.load(os.path.join(TargetDirectory(star, campaign), model + '.npz'))['cdpp']
          if (cdpp > cdpp_hi) or (cdpp < cdpp_lo):
            continue
    
      # Passed all the tests!
      targets.append(star)
    
      # Do we have enough? If so, return
      if len(targets) == neighbors:
        random.shuffle(targets)
        return targets
  
  # If we get to this point, we didn't find enough neighbors... 
  # Return what we have anyway.
  return targets
    
def Statistics(campaign = 0, clobber = False, model = 'nPLD', injection = False, compare_to = 'everest1', plot = True, **kwargs):
  '''
  
  '''
  
  # Is this an injection run?
  if injection:
    return InjectionStatistics(campaign = campaign, clobber = clobber, model = model, plot = plot, **kwargs)
  
  # Compute the statistics
  sub = np.array([s[0] for s in GetK2Campaign(campaign)], dtype = int)
  outfile = os.path.join(EVEREST_DAT, 'k2', 'stats', '%s_c%02d.cdpp' % (model, int(campaign)))
  if clobber or not os.path.exists(outfile):
    with open(outfile, 'w') as f:
      print("EPIC               Kp           Raw CDPP     Everest CDPP      Validation        Outliers[1]     Outliers[2]     Datapoints     Saturated", file = f)
      print("---------          ------       ---------    ------------      ----------        -----------     -----------     ----------     ---------", file = f)
      all = GetK2Campaign(int(campaign))
      stars = np.array([s[0] for s in all], dtype = int)
      kpmgs = np.array([s[1] for s in all], dtype = float)
      for i, _ in enumerate(stars):
        sys.stdout.write('\rProcessing target %d/%d...' % (i + 1, len(stars)))
        sys.stdout.flush()
        nf = os.path.join(EVEREST_DAT, 'k2', 'c%02d' % campaign, 
                         ('%09d' % stars[i])[:4] + '00000', 
                         ('%09d' % stars[i])[4:], model + '.npz')
        try:
          data = np.load(nf)
          
          # Remove NaNs and flagged cadences
          flux = np.delete(data['fraw'] - data['model'], np.array(list(set(np.concatenate([data['nanmask'], data['badmask']])))))
          # Iterative sigma clipping to get 5 sigma outliers
          inds = np.array([], dtype = int)
          m = 1
          while len(inds) < m:
            m = len(inds)
            ff = SavGol(np.delete(flux, inds))
            med = np.nanmedian(ff)
            MAD = 1.4826 * np.nanmedian(np.abs(ff - med))
            inds = np.append(inds, np.where((ff > med + 5. * MAD) | (ff < med - 5. * MAD))[0])
          nout = len(inds)
          ntot = len(flux)

          print("{:>09d} {:>15.3f} {:>15.3f} {:>15.3f} {:>15.3f} {:>15d} {:>15d} {:>15d} {:>15d}".format(stars[i], kpmgs[i], data['cdppr'][()], data['cdpp'][()], data['cdppv'][()], len(data['outmask']), nout, ntot, int(data['saturated'])), file = f)
        except:
          print("{:>09d} {:>15.3f} {:>15.3f} {:>15.3f} {:>15.3f} {:>15d} {:>15d} {:>15d} {:>15d}".format(stars[i], kpmgs[i], np.nan, np.nan, np.nan, 0, 0, 0, 0), file = f)
      print("")
  
  if plot:
  
    # Load all stars
    epic, kp, cdpp6r, cdpp6, cdpp6v, _, out, tot, saturated = np.loadtxt(outfile, unpack = True, skiprows = 2)
    epic = np.array(epic, dtype = int)
    out = np.array(out, dtype = int)
    tot = np.array(tot, dtype = int)
    saturated = np.array(saturated, dtype = int)

    # Get only stars in this subcampaign
    inds = np.array([e in sub for e in epic])
    epic = epic[inds]
    kp = kp[inds]
    # HACK: Campaign 0 magnitudes are reported only to the nearest tenth,
    # so let's add a little noise to spread them out for nicer plotting
    kp = kp + 0.1 * (0.5 - np.random.random(len(kp)))
    cdpp6r = cdpp6r[inds]
    cdpp6 = cdpp6[inds]
    cdpp6v = cdpp6v[inds]
    out = out[inds]
    tot = tot[inds]
    saturated = saturated[inds]
    sat = np.where(saturated == 1)
    unsat = np.where(saturated == 0)
    if not np.any([not np.isnan(x) for x in cdpp6]):
      raise Exception("No targets to plot.")

    # Control transparency
    alpha_kepler = 0.03
    alpha_unsat = min(0.1, 2000. / (1 + len(unsat[0])))
    alpha_sat = min(1., 180. / (1 + len(sat[0])))

    # Get the comparison model stats
    if compare_to.lower() == 'everest1':
      epic_1, cdpp6_1 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 
                                   'tables', 'c%02d_everest1.cdpp' % int(campaign)), unpack = True)
      cdpp6_1 = sort_like(cdpp6_1, epic, epic_1)   
      # Outliers
      epic_1, out_1, tot_1 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 
                                        'tables', 'c%02d_everest1.out' % int(campaign)), unpack = True)
      out_1 = sort_like(out_1, epic, epic_1)
      tot_1 = sort_like(tot_1, epic, epic_1)
    elif compare_to.lower() == 'k2sc':
      epic_1, cdpp6_1 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 
                                   'tables', 'c%02d_k2sc.cdpp' % int(campaign)), unpack = True)
      cdpp6_1 = sort_like(cdpp6_1, epic, epic_1)
      # Outliers
      epic_1, out_1, tot_1 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 
                                        'tables', 'c%02d_k2sc.out' % int(campaign)), unpack = True)
      out_1 = sort_like(out_1, epic, epic_1)
      tot_1 = sort_like(tot_1, epic, epic_1)
    elif compare_to.lower() == 'k2sff':
      epic_1, cdpp6_1 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 
                                   'tables', 'c%02d_k2sff.cdpp' % int(campaign)), unpack = True)
      cdpp6_1 = sort_like(cdpp6_1, epic, epic_1)
      # Outliers
      epic_1, out_1, tot_1 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 
                                        'tables', 'c%02d_k2sff.out' % int(campaign)), unpack = True)
      out_1 = sort_like(out_1, epic, epic_1)
      tot_1 = sort_like(tot_1, epic, epic_1)
    elif compare_to.lower() == 'kepler':
      kic, kepler_kp, kepler_cdpp6 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 
                                         'tables', 'kepler.cdpp'), unpack = True)
    else:
      compfile = os.path.join(EVEREST_DAT, 'k2', 'stats', '%s_c%02d.cdpp' % (compare_to, int(campaign)))
      epic_1, _, _, cdpp6_1, _, _, _ = np.loadtxt(compfile, unpack = True, skiprows = 2)
      epic_1 = np.array(epic_1, dtype = int)
      inds = np.array([e in sub for e in epic_1])
      epic_1 = epic_1[inds]
      cdpp6_1 = cdpp6_1[inds]
      cdpp6_1 = sort_like(cdpp6_1, epic, epic_1) 
 
    # ------ 1. Plot cdpp vs. mag
    if compare_to.lower() != 'kepler':
      fig = pl.figure(figsize = (16, 5))
      ax = [pl.subplot2grid((120, 120), (0,  0), colspan=35, rowspan=120),
            pl.subplot2grid((120, 120), (0,  40), colspan=35, rowspan=120),
            pl.subplot2grid((120, 120), (0,  80), colspan=35, rowspan=55),
            pl.subplot2grid((120, 120), (65,  80), colspan=35, rowspan=55)]
    else:
      fig = pl.figure(figsize = (12, 5))
      ax = [pl.subplot2grid((120, 75), (0,  0), colspan=35, rowspan=120),
            None,
            pl.subplot2grid((120, 75), (0,  40), colspan=35, rowspan=55),
            pl.subplot2grid((120, 75), (65,  40), colspan=35, rowspan=55)]
    fig.canvas.set_window_title('K2 Campaign %s: %s versus %s' % (campaign, model, compare_to))
    fig.subplots_adjust(left = 0.05, right = 0.95, bottom = 0.125, top = 0.9)
    bins = np.arange(7.5,18.5,0.5)
    if compare_to.lower() != 'kepler':
      ax[0].scatter(kp[unsat], cdpp6_1[unsat], color = 'y', marker = '.', alpha = alpha_unsat)
      ax[0].scatter(kp[sat], cdpp6_1[sat], color = 'y', marker = 's', alpha = alpha_sat, s = 5)
      ax[0].scatter(kp[unsat], cdpp6[unsat], color = 'b', marker = '.', alpha = alpha_unsat, picker = True)
      ax[0].scatter(kp[sat], cdpp6[sat], color = 'b', marker = 's', alpha = alpha_sat, s = 5, picker = True)
      for y, style in zip([cdpp6_1, cdpp6], ['yo', 'bo']):
        by = np.zeros_like(bins) * np.nan
        for b, bin in enumerate(bins):
          i = np.where((y > -np.inf) & (y < np.inf) & (kp >= bin - 0.5) & (kp < bin + 0.5))[0]
          if len(i) > 10:
            by[b] = np.median(y[i])
        ax[0].plot(bins, by, style)
    else:  
      ax[0].scatter(kepler_kp, kepler_cdpp6, color = 'y', marker = '.', alpha = alpha_kepler)
      ax[0].scatter(kp, cdpp6, color = 'b', marker = '.', alpha = alpha_unsat, picker = True)
      for x, y, style in zip([kepler_kp, kp], [kepler_cdpp6, cdpp6], ['yo', 'bo']):
        by = np.zeros_like(bins) * np.nan
        for b, bin in enumerate(bins):
          i = np.where((y > -np.inf) & (y < np.inf) & (x >= bin - 0.5) & (x < bin + 0.5))[0]
          if len(i) > 10:
            by[b] = np.median(y[i])
        ax[0].plot(bins, by, style)
    ax[0].set_ylim(-10, 500)
    ax[0].set_xlim(8, 18)
    ax[0].set_xlabel('Kepler Magnitude', fontsize = 18)
    ax[0].set_title('CDPP6 (ppm)', fontsize = 18)

    # ------ 2. Plot the equivalent of the Aigrain+16 figure
    if compare_to.lower() != 'kepler':
      x = kp
      y = (cdpp6 - cdpp6_1) / cdpp6_1
      yv = (cdpp6v - cdpp6_1) / cdpp6_1
      ax[1].scatter(x[unsat], y[unsat], color = 'b', marker = '.', alpha = alpha_unsat, zorder = -1, picker = True)
      ax[1].scatter(x[sat], y[sat], color = 'r', marker = '.', alpha = alpha_sat, zorder = -1, picker = True)
      ax[1].set_ylim(-1,1)
      ax[1].set_xlim(8,18)
      ax[1].axhline(0, color = 'gray', lw = 2, zorder = -99, alpha = 0.5)
      ax[1].axhline(0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)
      ax[1].axhline(-0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)
      bins = np.arange(7.5,18.5,0.5)
      # Bin the CDPP
      by = np.zeros_like(bins) * np.nan
      for b, bin in enumerate(bins):
        i = np.where((y > -np.inf) & (y < np.inf) & (x >= bin - 0.5) & (x < bin + 0.5))[0]
        if len(i) > 10:
          by[b] = np.median(y[i])
      ax[1].plot(bins[:9], by[:9], 'k--', lw = 2)
      ax[1].plot(bins[8:], by[8:], 'k-', lw = 2)
      ax[1].set_title(r'Relative CDPP', fontsize = 18)
      ax[1].set_xlabel('Kepler Magnitude', fontsize = 18)
      
    # ------ 3. Plot the outliers
    i = np.argsort(out)
    a = int(0.95 * len(out))
    omax = out[i][a]
    if compare_to.lower() != 'kepler':
      j = np.argsort(out_1)
      b = int(0.95 * len(out_1))
      omax = max(omax, out_1[j][b])
    ax[2].hist(out, 25, range = (0, omax), histtype = 'step', color = 'b')
    if compare_to.lower() != 'kepler':
      ax[2].hist(out_1, 25, range = (0, omax), histtype = 'step', color = 'y')
    ax[2].margins(0, None)
    ax[2].set_title('Number of Outliers', fontsize = 18)
    
    # Plot the total number of data points
    i = np.argsort(tot)
    a = int(0.05 * len(tot))
    b = int(0.95 * len(tot))
    tmin = tot[i][a]
    tmax = tot[i][b]
    if compare_to.lower() != 'kepler':
      j = np.argsort(tot_1)
      c = int(0.05 * len(tot_1))
      d = int(0.95 * len(tot_1))
      tmin = min(tmin, tot_1[j][c])
      tmax = max(tmax, tot_1[j][d])
    ax[3].hist(tot, 25, range = (tmin, tmax), histtype = 'step', color = 'b')
    if compare_to.lower() != 'kepler':
      ax[3].hist(tot_1, 25, range = (tmin, tmax), histtype = 'step', color = 'y')
    ax[3].margins(0, None)
    ax[3].set_xlabel('Number of Data Points', fontsize = 18)
    
    # Pickable points
    Picker = StatsPicker([ax[0], ax[1]], [kp, kp], [cdpp6, y], epic, campaign, model = model, compare_to = compare_to)
    fig.canvas.mpl_connect('pick_event', Picker)
    
    # Show
    pl.show()

def HasShortCadence(EPIC, season = None):
  '''
  Returns `True` if short cadence data is available for this target.
  
  '''
  
  if season is None:
    season = Campaign(EPIC)
    if season is None:
      return None
  stars = GetK2Campaign(season)
  i = np.where([s[0] == EPIC for s in stars])[0]
  if len(i):
    return stars[i[0]][3]
  else:
    return None

def InjectionStatistics(campaign = 0, clobber = False, model = 'nPLD', plot = True, show = True, **kwargs):
  '''
  
  '''
  
  # Compute the statistics
  stars = GetK2Campaign(campaign, epics_only = True)
  if type(campaign) is int:
    outfile = os.path.join(EVEREST_DAT, 'k2', 'stats', '%s_c%02d.inj' % (model, campaign))
  else:
    outfile = os.path.join(EVEREST_DAT, 'k2', 'stats', '%s_c%04.1f.inj' % (model, campaign))
  if clobber or not os.path.exists(outfile):
    with open(outfile, 'w') as f:
      print("EPIC         Depth         UControl      URecovered    MControl      MRecovered", file = f)
      print("---------    ----------    ----------    ----------    ----------    ----------", file = f)
      for i, _ in enumerate(stars):
        sys.stdout.write('\rProcessing target %d/%d...' % (i + 1, len(stars)))
        sys.stdout.flush()
        path = os.path.join(EVEREST_DAT, 'k2', 'c%02d' % int(campaign), 
                           ('%09d' % stars[i])[:4] + '00000', 
                           ('%09d' % stars[i])[4:])

        # Loop over all depths
        for depth in [0.01, 0.001, 0.0001]:
          
          try:
            
            # Unmasked
            data = np.load(os.path.join(path, '%s_Inject_U%g.npz' % (model, depth)))
            assert depth == data['inject'][()]['depth'], ""
            ucontrol = data['inject'][()]['rec_depth_control']
            urecovered = data['inject'][()]['rec_depth']
        
            # Masked
            data = np.load(os.path.join(path, '%s_Inject_M%g.npz' % (model, depth)))
            assert depth == data['inject'][()]['depth'], ""
            mcontrol = data['inject'][()]['rec_depth_control']
            mrecovered = data['inject'][()]['rec_depth']
        
            # Log it
            print("{:>09d} {:>13.8f} {:>13.8f} {:>13.8f} {:>13.8f} {:>13.8f}".format(stars[i], depth, ucontrol, urecovered, mcontrol, mrecovered), file = f)
        
          except:
            pass
          
      print("")
  
  if plot:
  
    # Load the statistics
    try:
      epic, depth, ucontrol, urecovered, mcontrol, mrecovered = np.loadtxt(outfile, unpack = True, skiprows = 2)
    except ValueError:
      raise Exception("No targets to plot.")
    
    # Normalize to the injected depth
    ucontrol /= depth
    urecovered /= depth
    mcontrol /= depth
    mrecovered /= depth
    
    # Set up the plot
    fig, ax = pl.subplots(3,2, figsize = (9,12))
    fig.subplots_adjust(hspace = 0.29)
    ax[0,0].set_title(r'Default', fontsize = 18)
    ax[0,1].set_title(r'Masked', fontsize = 18)
    ax[0,0].set_ylabel(r'$D_0 = 10^{-2}$', rotation = 90, fontsize = 18, labelpad = 10)
    ax[1,0].set_ylabel(r'$D_0 = 10^{-3}$', rotation = 90, fontsize = 18, labelpad = 10)
    ax[2,0].set_ylabel(r'$D_0 = 10^{-4}$', rotation = 90, fontsize = 18, labelpad = 10)
    
    # Define some useful stuff for plotting
    depths = [1e-2, 1e-3, 1e-4]
    ranges = [(0.75, 1.25), (0.5, 1.5), (0., 2.)]
    nbins = [30, 30, 20]
    ymax = [0.4, 0.25, 0.16]
    xticks = [[0.75, 0.875, 1., 1.125, 1.25], [0.5, 0.75, 1., 1.25, 1.5], [0., 0.5, 1., 1.5, 2.0]]
    
    # Plot
    for i in range(3): 
      
      # Indices for this plot
      idx = np.where(depth == depths[i])
      
      for j, control, recovered in zip([0, 1], [ucontrol[idx], mcontrol[idx]], [urecovered[idx], mrecovered[idx]]):
      
        # Control
        ax[i,j].hist(control, bins = nbins[i], range = ranges[i], color = 'r', 
                     histtype = 'step', weights = np.ones_like(control) / len(control))
  
        # Recovered
        ax[i,j].hist(recovered, bins = nbins[i], range = ranges[i], color = 'b', 
                     histtype = 'step', weights = np.ones_like(recovered) / len(recovered))
      
        # Indicate center
        ax[i,j].axvline(1., color = 'k', ls = '--')
  
        # Indicate the fraction above and below
        if len(recovered):
          au = len(np.where(recovered > ranges[i][1])[0]) / len(recovered)
          al = len(np.where(recovered < ranges[i][0])[0]) / len(recovered)
          ax[i,j].annotate('%.2f' % al, xy = (0.01, 0.93), xycoords = 'axes fraction', 
                        xytext = (0.1, 0.93), ha = 'left', va = 'center', color = 'b',
                        arrowprops = dict(arrowstyle="->",color='b'))
          ax[i,j].annotate('%.2f' % au, xy = (0.99, 0.93), xycoords = 'axes fraction', 
                        xytext = (0.9, 0.93), ha = 'right', va = 'center', color = 'b',
                        arrowprops = dict(arrowstyle="->",color='b'))
        if len(control):  
          cu = len(np.where(control > ranges[i][1])[0]) / len(control)
          cl = len(np.where(control < ranges[i][0])[0]) / len(control)
          ax[i,j].annotate('%.2f' % cl, xy = (0.01, 0.86), xycoords = 'axes fraction', 
                        xytext = (0.1, 0.86), ha = 'left', va = 'center', color = 'r',
                        arrowprops = dict(arrowstyle="->",color='r'))
          ax[i,j].annotate('%.2f' % cu, xy = (0.99, 0.86), xycoords = 'axes fraction', 
                        xytext = (0.9, 0.86), ha = 'right', va = 'center', color = 'r',
                        arrowprops = dict(arrowstyle="->",color='r'))
                
        # Indicate the median
        if len(recovered):
          ax[i,j].annotate('M = %.2f' % np.median(recovered), xy = (0.35, 0.5), ha = 'right',
                        xycoords = 'axes fraction', color = 'b', fontsize = 16)
        if len(control):
          ax[i,j].annotate('M = %.2f' % np.median(control), xy = (0.65, 0.5), ha = 'left',
                        xycoords = 'axes fraction', color = 'r', fontsize = 16)
  
        # Tweaks
        ax[i,j].set_xticks(xticks[i])
        ax[i,j].set_xlim(xticks[i][0], xticks[i][-1])
        ax[i,j].set_ylim(-0.005, ymax[i])
        ax[i,j].set_xlabel(r'$D/D_0$', fontsize = 16)
        
        ax[i,j].get_yaxis().set_major_locator(MaxNLocator(5))
        for tick in ax[i,j].get_xticklabels() + ax[i,j].get_yticklabels():
          tick.set_fontsize(14)
        
    if show:
      pl.show()
    else:
      return fig, ax

def HDUCards(headers, hdu = 0):
  '''
  
  '''
  
  if hdu == 0:
    # Get info from the TPF Primary HDU Header
    tpf_header = headers[0]
    entries = ['TELESCOP', 'INSTRUME', 'OBJECT', 'KEPLERID', 'CHANNEL', 'MODULE', 
               'OUTPUT', 'CAMPAIGN', 'DATA_REL', 'OBSMODE', 'TTABLEID',
               'RADESYS', 'RA_OBJ', 'DEC_OBJ',  'EQUINOX', 'KEPMAG']
  elif (hdu == 1) or (hdu == 6):
    # Get info from the TPF BinTable HDU Header
    tpf_header = headers[1]
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
  elif (hdu == 3) or (hdu == 4) or (hdu == 5):
    # Get info from the TPF BinTable HDU Header
    tpf_header = headers[2]
    entries = ['TELESCOP', 'INSTRUME', 'OBJECT', 'KEPLERID', 'RADESYS', 'RA_OBJ', 
               'DEC_OBJ', 'EQUINOX', 'WCSAXES', 'CTYPE1', 'CTYPE2', 'CRPIX1', 
               'CRPIX2', 'CRVAL1', 'CRVAL2', 'CUNIT1', 'CUNIT2', 'CDELT1', 
               'CDELT2', 'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2', 'WCSNAMEP', 'WCSAXESP', 
               'CTYPE1P', 'CUNIT1P', 'CRPIX1P', 'CRVAL1P', 'CDELT1P', 'CTYPE2P', 
               'CUNIT2P', 'CRPIX2P', 'CRVAL2P', 'CDELT2P', 'NPIXSAP', 'NPIXMISS']
  else:
    return []
        
  cards = []
  cards.append(('COMMENT', '************************'))
  cards.append(('COMMENT', '*     MISSION INFO     *'))
  cards.append(('COMMENT', '************************'))
  for entry in entries:
    try:
      cards.append(tuple(tpf_header[entry]))
    except KeyError:
      pass
  return cards  

def TargetDirectory(ID, season, relative = False, **kwargs):
  '''
  
  '''
  
  if relative:
    path = ''
  else:
    path = EVEREST_DAT
  return os.path.join(path, 'k2', 'c%02d' % season, 
                     ('%09d' % ID)[:4] + '00000', 
                     ('%09d' % ID)[4:])

def FITSFile(ID, season, cadence = 'lc'):
  '''
  
  '''
  
  if cadence == 'lc':
    return 'hlsp_everest_k2_llc_%d-c%02d_kepler_v%s_lc.fits' % (ID, season, EVEREST_VERSION)
  else:
    return 'hlsp_everest_k2_llc_%d-c%02d_kepler_v%s_sclc.fits' % (ID, season, EVEREST_VERSION)
    
def FITSUrl(ID, season):
  '''
  
  '''
  
  mast_version = MAST_VERSION()
  assert mast_version == EVEREST_VERSION, "Version %s light curves not available on MAST." % EVEREST_VERSION
  url = MAST_ROOT + 'c%02d/' % season + ('%09d' % ID)[:4] + '00000/' + ('%09d/' % ID)[4:]
  return url
