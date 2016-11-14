#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`k2.py` - Main mission routines
---------------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .aux import *
from ...config import EVEREST_SRC, EVEREST_DAT, EVEREST_DEV
from ...utils import DataContainer, sort_like, AP_COLLAPSED_PIXEL, AP_SATURATED_PIXEL
from ...math import SavGol, Interpolate
try:
  import pyfits
except ImportError:
  try:
    import astropy.io.fits as pyfits
  except ImportError:
    raise Exception('Please install the `pyfits` package.')
import matplotlib.pyplot as pl
from matplotlib.ticker import ScalarFormatter
import k2plr as kplr; kplr_client = kplr.API()
from k2plr.config import KPLR_ROOT
import numpy as np
from tempfile import NamedTemporaryFile
import random
import os, sys, shutil
import logging
log = logging.getLogger(__name__)

__all__ = ['Setup', 'Season', 'Breakpoint', 'GetData', 'GetNeighbors', 
           'Statistics', 'TargetDirectory', 'HasShortCadence', 'InjectionStatistics']

def Setup():
  '''
  
  '''
  
  if not os.path.exists(os.path.join(EVEREST_DAT, 'k2', 'stats')):
    os.makedirs(os.path.join(EVEREST_DAT, 'k2', 'stats'))
  GetK2Stars(clobber = False)

def Season(EPIC, **kwargs):
  '''
  Returns the campaign number for a given EPIC target.
  
  '''
  
  return Campaign(EPIC, **kwargs)

def Breakpoint(EPIC, **kwargs):
  '''
  .. todo:: Populate these.
  
  '''
  
  campaign = Season(EPIC)
  breakpoints = {0: None,
                 1: None,
                 2: 2042,
                 3: None,
                 4: None,
                 5: None,
                 6: 2143,
                 7: None,
                 8: None,
                 9: None,
                 10: None,
                 11: None,
                 12: None,
                 13: None,
                 14: None,
                 15: None,
                 16: None,
                 17: None}
  if campaign in breakpoints:
    return breakpoints[campaign]
  else:
    return None

def GetData(EPIC, season = None, clobber = False, delete_raw = False, 
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
                  if tpf_big_aperture[n[0]][n[1]] == 1:
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
      sc_naninds = None
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
                        sc_fpix_err = sc_fpix_err, sc_qual = sc_qual, sc_naninds = sc_naninds,
                        sc_pc1 = sc_pc1, sc_pc2 = sc_pc2, sc_fitsheader = sc_fitsheader)
    f.flush()
    os.fsync(f.fileno())
    f.close()
    shutil.move(f.name, filename)
    
    if download_only:
      return
  
  # Load
  else:
    data = np.load(filename)
    cadn = data['cadn']
    time = data['time']
    fpix = data['fpix']
    fpix_err = data['fpix_err']
    qual = data['qual']
    apertures = data['apertures'][()]
    pc1 = data['pc1']
    pc2 = data['pc2']
    fitsheader = data['fitsheader']
    pixel_images = data['pixel_images']
    nearby = data['nearby']
    hires = data['hires']    
    sc_cadn = data['sc_cadn']
    sc_time = data['sc_time']
    sc_fpix = data['sc_fpix']
    sc_fpix_err = data['sc_fpix_err']
    sc_qual = data['sc_qual']
    sc_naninds = data['sc_naninds']
    sc_pc1 = data['sc_pc1']
    sc_pc2 = data['sc_pc2']
    sc_fitsheader = data['sc_fitsheader']
  
  # Select the "saturated aperture" to check if the star is saturated
  # If it is, we will use this aperture instead
  if saturated_aperture_name == 'custom':
    saturated_aperture = GetCustomAperture(data)
  else:
    if saturated_aperture_name is None:
      saturated_aperture_name = 'k2sff_19'
    saturated_aperture = apertures[saturated_aperture_name]
    assert saturated_aperture is not None, "Invalid aperture selected."
    
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
    assert aperture is not None, "Invalid aperture selected."

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
    if short_cadence:
      sc_fpixnew = []
      sc_ferrnew = []

    for j in range(aperture.shape[1]):
      if np.any(f97[:,j] > satflx):
        marked = False
        collapsed = np.zeros(len(fpix[:,0,0]))
        collapsed_err2 = np.zeros(len(fpix[:,0,0]))
        if short_cadence:
          sc_collapsed = np.zeros(len(sc_fpix[:,0,0]))
          sc_collapsed_err2 = np.zeros(len(sc_fpix[:,0,0]))
        for i in range(aperture.shape[0]):
          if aperture[i,j]:
            if not marked:
              aperture[i,j] = AP_COLLAPSED_PIXEL
              marked = True
            else:
              aperture[i,j] = AP_SATURATED_PIXEL
            collapsed += fpix[:,i,j]
            collapsed_err2 += fpix_err[:,i,j] ** 2
            if short_cadence:
              sc_collapsed += sc_fpix[:,i,j]
              sc_collapsed_err2 += sc_fpix_err[:,i,j] ** 2
        if np.any(collapsed):
          fpixnew.append(collapsed)
          ferrnew.append(np.sqrt(collapsed_err2))
          if short_cadence:
            sc_fpixnew.append(sc_collapsed)
            sc_ferrnew.append(np.sqrt(sc_collapsed_err2))
          ncol += 1
      else:
        for i in range(aperture.shape[0]):
          if aperture[i,j]:
            fpixnew.append(fpix[:,i,j])
            ferrnew.append(fpix_err[:,i,j])
            if short_cadence:
              sc_fpixnew.append(sc_fpix[:,i,j])
              sc_ferrnew.append(sc_fpix_err[:,i,j])
    fpix2D = np.array(fpixnew).T
    fpix_err2D = np.array(ferrnew).T
    if short_cadence:
      sc_fpix2D = np.array(sc_fpixnew).T
      sc_fpix_err2D = np.array(sc_ferrnew).T
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
    if short_cadence:
      sc_fpix2D = np.array([f[ap] for f in sc_fpix], dtype='float64')
      sc_fpix_err2D = np.array([p[ap] for p in sc_fpix_err], dtype='float64')
    
  # Compute the background
  binds = np.where(aperture ^ 1)
  if RemoveBackground(EPIC) and (len(binds[0]) > 0):
    bkg = np.nanmedian(np.array([f[binds] for f in fpix], dtype='float64'), axis = 1)
    # Uncertainty of the median: http://davidmlane.com/hyperstat/A106993.html
    bkg_err = 1.253 * np.nanmedian(np.array([e[binds] for e in fpix_err], 
                      dtype='float64'), axis = 1) / np.sqrt(len(binds[0]))
    bkg = bkg.reshape(-1, 1)
    bkg_err = bkg_err.reshape(-1, 1)
    if short_cadence:
      sc_bkg = np.nanmedian(np.array([f[binds] for f in sc_fpix], dtype='float64'), axis = 1)
      sc_bkg_err = 1.253 * np.nanmedian(np.array([e[binds] for e in sc_fpix_err], 
                           dtype='float64'), axis = 1) / np.sqrt(len(binds[0]))
      sc_bkg = sc_bkg.reshape(-1, 1)
      sc_bkg_err = sc_bkg_err.reshape(-1, 1)
  else:
    bkg = 0.
    bkg_err = 0.
    if short_cadence:
      sc_bkg = 0.
      sc_bkg_err = 0.
  
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
  t = np.delete(time, badmask)
  f = np.delete(flux, badmask)
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

  # Repeat for short cadence
  if short_cadence:
    sc_fpix = sc_fpix2D - sc_bkg
    sc_fpix_err = np.sqrt(sc_fpix_err2D ** 2 + sc_bkg_err ** 2)
    sc_flux = np.sum(sc_fpix, axis = 1)
    sc_ferr = np.sqrt(np.sum(sc_fpix_err ** 2, axis = 1))
    
    # NaN and bad masks
    sc_nanmask = np.where(np.isnan(sc_flux) | (sc_flux == 0))[0]                         
    sc_badmask = []
    for b in bad_bits:
      sc_badmask += list(np.where(sc_qual & 2 ** (b - 1))[0])
  
    # Flag >10 sigma outliers -- same thing.
    t = np.delete(sc_time, sc_badmask)
    f = np.delete(sc_flux, sc_badmask)
    f = SavGol(f)
    med = np.nanmedian(f)
    MAD = 1.4826 * np.nanmedian(np.abs(f - med))
    bad = np.where((f > med + 10. * MAD) | (f < med - 10. * MAD))[0]
    sc_badmask.extend([np.argmax(sc_time == t[i]) for i in bad])
  
    # Campaign 2 hack: the first day or two are screwed up
    if campaign == 2:
      sc_badmask.extend(np.where(sc_time < 2061.5)[0])
  
    # TODO: Fix time offsets in first half of 
    # Campaign 0. See note in everest 1.0 code
  
    # Finalize the mask
    sc_badmask = np.array(sorted(list(set(sc_badmask))))
  
    # Interpolate the nans
    sc_fpix = Interpolate(sc_time, sc_nanmask, sc_fpix)
    sc_fpix_err = Interpolate(sc_time, sc_nanmask, sc_fpix_err)

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
  
  if short_cadence:
    data.sc_cadn = sc_cadn
    data.sc_time = sc_time
    data.sc_fpix = sc_fpix
    data.sc_fpix_err = sc_fpix_err
    data.sc_nanmask = sc_nanmask
    data.sc_badmask = sc_badmask
    data.sc_quality = sc_qual
    data.sc_Xpos = sc_pc1
    data.sc_Ypos = sc_pc2
    data.sc_meta = sc_fitsheader
    data.sc_bkg = sc_bkg
  else:
    data.sc_cadn = None
    data.sc_time = None
    data.sc_fpix = None
    data.sc_fpix_err = None
    data.sc_nanmask = None
    data.sc_badmask = None
    data.sc_quality = None
    data.sc_Xpos = None
    data.sc_Ypos = None
    data.sc_meta = None
    data.sc_bkg = None
    
  return data

def GetNeighbors(EPIC, model = None, neighbors = 10, mag_range = (11., 13.), 
                 cdpp_range = None, aperture_name = 'k2sff_15', **kwargs):
  '''
  Return `neighbors` random bright stars on the same module as `EPIC`.
  
  '''
  
  # Get the IDs
  campaign = Season(EPIC)
  has_short_cadence = HasShortCadence(EPIC, season = campaign)
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
    if has_short_cadence:
      mag_lo = 8.
  if cdpp_range is None:
    cdpp_lo = -np.inf
    cdpp_hi = np.inf
  else:
    cdpp_lo = cdpp_range[0]
    cdpp_hi = cdpp_range[1]
  targets = []
  
  # First look for nearby targets, then relax the constraint
  for nearby in [True, False]:
  
    # Loop over all stars
    for star, kp, channel, sc in zip(epics, kepmags, channels, short_cadence):
    
      # Preliminary vetting
      if not (((channel in c) if nearby else True) and (kp < mag_hi) and (kp > mag_lo) and (sc if has_short_cadence else True)):
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
          cdpp = np.load(os.path.join(TargetDirectory(star, campaign), model + '.npz'))['cdpp6']
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
  outfile = os.path.join(EVEREST_DAT, 'k2', 'stats', '%s_c%02d.tsv' % (model, int(campaign)))
  if clobber or not os.path.exists(outfile):
    with open(outfile, 'w') as f:
      print("EPIC               Kp           Raw CDPP     Everest CDPP      Validation        Outliers       Saturated", file = f)
      print("---------          ------       ---------    ------------      ----------        --------       ---------", file = f)
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
          print("{:>09d} {:>15.3f} {:>15.3f} {:>15.3f} {:>15.3f} {:>15d} {:>15d}".format(stars[i], kpmgs[i], data['cdppr'][()], data['cdpp6'][()], data['cdppv'][()], len(data['outmask']), int(data['saturated'])), file = f)
        except:
          print("{:>09d} {:>15.3f} {:>15.3f} {:>15.3f} {:>15.3f} {:>15d} {:>15d}".format(stars[i], kpmgs[i], np.nan, np.nan, np.nan, 0, 0), file = f)
      print("")
  
  if plot:
  
    # Load all stars
    epic, kp, cdpp6r, cdpp6, cdpp6v, outliers, saturated = np.loadtxt(outfile, unpack = True, skiprows = 2)
    epic = np.array(epic, dtype = int)
    outliers = np.array(outliers, dtype = int)
    saturated = np.array(saturated, dtype = int)

    # Get only stars in this subcampaign
    inds = np.array([e in sub for e in epic])
    epic = epic[inds]
    kp = kp[inds]
    cdpp6r = cdpp6r[inds]
    cdpp6 = cdpp6[inds]
    cdpp6v = cdpp6v[inds]
    outliers = outliers[inds]
    saturated = saturated[inds]
    sat = np.where(saturated == 1)
    unsat = np.where(saturated == 0)
    
    # Control transparency
    alpha_kepler = 0.03
    alpha_unsat = min(0.1, 2000. / (1 + len(unsat[0])))
    alpha_sat = min(1., 20. / (1 + len(sat[0])))

    # Get the comparison model stats
    if compare_to.lower() == 'everest1':
      epic_1, cdpp6_1 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 
                                   'tables', 'c%02d_everest1.tsv' % int(campaign)), unpack = True)
      cdpp6_1 = sort_like(cdpp6_1, epic, epic_1)   
    elif compare_to.lower() == 'k2sc':
      epic_1, cdpp6_1 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 
                                   'tables', 'c%02d_k2sc.tsv' % int(campaign)), unpack = True)
      cdpp6_1 = sort_like(cdpp6_1, epic, epic_1)
    elif compare_to.lower() == 'k2sff':
      epic_1, cdpp6_1 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 
                                   'tables', 'c%02d_k2sff.tsv' % int(campaign)), unpack = True)
      cdpp6_1 = sort_like(cdpp6_1, epic, epic_1)
    elif compare_to.lower() == 'kepler':
      kic, kepler_kp, kepler_cdpp6 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 
                                         'tables', 'kepler.tsv'), unpack = True)
    else:
      compfile = os.path.join(EVEREST_DAT, 'k2', 'stats', '%s_c%02d.tsv' % (compare_to, int(campaign)))
      epic_1, _, _, cdpp6_1, _, _, _ = np.loadtxt(compfile, unpack = True, skiprows = 2)
      epic_1 = np.array(epic_1, dtype = int)
      inds = np.array([e in sub for e in epic_1])
      epic_1 = epic_1[inds]
      cdpp6_1 = cdpp6_1[inds]
      cdpp6_1 = sort_like(cdpp6_1, epic, epic_1) 
 
    # ------ 1. Plot cdpp vs. mag
    fig = pl.figure(figsize = (16, 5))
    fig.canvas.set_window_title('K2 Campaign %s: %s versus %s' % (campaign, model, compare_to))
    fig.subplots_adjust(left = 0.05, right = 0.95, bottom = 0.125, top = 0.9)
    if compare_to.lower() != 'kepler':
      ax = [pl.subplot2grid((120, 120), (0,  0), colspan=35, rowspan=120),
            pl.subplot2grid((120, 120), (0,  40), colspan=35, rowspan=120),
            pl.subplot2grid((120, 120), (0,  80), colspan=29, rowspan=120),
            pl.subplot2grid((120, 120), (0,  110), colspan=10, rowspan=120)]
    else:
      ax = [pl.subplot2grid((120, 120), (0,  0), colspan=70, rowspan=120),
            None,
            pl.subplot2grid((120, 120), (0,  75), colspan=34, rowspan=120),
            pl.subplot2grid((120, 120), (0,  110), colspan=10, rowspan=120)]
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
    ax[2].scatter(kp[unsat], outliers[unsat], color = 'k', marker = '.', alpha = alpha_unsat, picker = True)
    ax[2].scatter(kp[sat], outliers[sat], color = 'r', marker = '.', alpha = alpha_sat, picker = True)
    ax[3].hist(np.log10(outliers[outliers > 0]), 50, orientation = "horizontal", color = 'k', alpha = 0.5)
    ax[2].set_title('Number of Outliers', fontsize = 18)
    ax[2].set_xlabel('Kepler Magnitude', fontsize = 18)
    ax[2].set_xlim(8, 18)
    ax[2].set_yscale('log')
    ax[2].set_ylim(0.9, 1e3)
    ax[2].yaxis.set_major_formatter(ScalarFormatter())
    ax[3].set_ylim(*np.log10(ax[2].get_ylim()))
    ax[3].set_xticks([])
    ax[3].set_yticks([])
    ax3t = ax[3].twinx()
    ax3t.set_yscale('log')
    ax3t.set_ylim(0.9, 1e3)
    ax3t.yaxis.set_ticklabels([])
    
    # Pickable points
    Picker = StatsPicker([ax[0], ax[1], ax[2]], [kp, kp, kp], [cdpp6, y, outliers], epic, campaign, model = model, compare_to = compare_to)
    fig.canvas.mpl_connect('pick_event', Picker)
    
    # Show
    pl.show()

def TargetDirectory(ID, season, **kwargs):
  '''
  
  '''

  return os.path.join(EVEREST_DAT, 'k2', 'c%02d' % season, 
                     ('%09d' % ID)[:4] + '00000', 
                     ('%09d' % ID)[4:])

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

def InjectionStatistics(campaign = 0, clobber = False, model = 'nPLD', plot = True, **kwargs):
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
      print("EPIC               Depth              UControl           URecovered         MControl           MRecovered", file = f)
      print("----               -----              --------           ----------         --------           ----------", file = f)
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
            data = np.load(os.path.join(path, '%s_Inject_U%g' % (model, depth))
            assert depth == data['inject'][()]['depth'], ""
            ucontrol = data['inject'][()]['rec_depth_control']
            urecovered = data['inject'][()]['rec_depth']
        
            # Masked
            data = np.load(os.path.join(path, '%s_Inject_M%g' % (model, depth))
            assert depth == data['inject'][()]['depth'], ""
            mcontrol = data['inject'][()]['rec_depth_control']
            mrecovered = data['inject'][()]['rec_depth']
        
            # Log it
            print("{:>09d} {:>13.8f} {:>13.8f} {:>13.8f}".format(stars[i], depth, ucontrol, urecovered, mcontrol, mrecovered), file = f)
        
          except:
            pass
          
      print("")
  
  if plot:
  
    # Load the statistics
    epic, depth, ucontrol, urecovered, mcontrol, mrecovered = np.loadtxt(outfile, unpack = True, skiprows = 2)
    
    # Normalize to the injected depth
    ucontrol /= depth
    urecovered /= depth
    mcontrol /= depth
    mrecovered /= depth
    
    # Set up the plot
    fig, ax = pl.subplots(3,2, figsize = (10,12))
    fig.subplots_adjust(hspace = 0.25)
    ax[0,0].set_title(r'Default', fontsize = 18)
    ax[0,1].set_title(r'Masked', fontsize = 18)
    ax[0,0].set_ylabel(r'D$_0$ = 10$^{-2}$', rotation = 90, fontsize = 18, labelpad = 10)
    ax[1,0].set_ylabel(r'D$_0$ = 10$^{-3}$', rotation = 90, fontsize = 18, labelpad = 10)
    ax[2,0].set_ylabel(r'D$_0$ = 10$^{-4}$', rotation = 90, fontsize = 18, labelpad = 10)
    
    # Define some useful stuff for plotting
    depths = [1e-2, 1e-3, 1e-4]
    ranges = [(0.5, 1.5), (0.5, 1.5), (0.5, 1.5), (0.5, 1.5), (0., 2.), (0., 2.)]
    nbins = [30, 30, 30]
    ymax = [0.6, 0.4, 0.15]
    xticks = [[0.5, 0.75, 1., 1.25, 1.5], [0.5, 0.75, 1., 1.25, 1.5], [0., 0.5, 1., 1.5, 2.0]]
    
    # Plot
    for i in range(3): 
      
      # Indices for this plot
      idx = np.where(depth == depths[i])
      
      for j, control, recovered in zip([0, 1], [ucontrol[idx], mcontrol[idx]], [urecovered[idx], mrecovered[idx]]):
      
        # Control
        axis.hist(control, bins = nbins[i], range = ranges[i], color = 'r', 
                  histtype = 'step', weights = np.ones_like(control) / len(control))
  
        # Recovered
        axis.hist(recovered, bins = nbins[i], range = ranges[i], color = 'b', 
                  histtype = 'step', weights = np.ones_like(recovered) / len(recovered))
      
        # Indicate center
        axis.axvline(1., color = 'k', ls = '--')
  
        # Indicate the fraction above and below
        if len(recovered):
          au = len(np.where(recovered > ranges[i][1])[0]) / len(recovered)
          al = len(np.where(recovered < ranges[i][0])[0]) / len(recovered)
          axis.annotate('%.2f' % al, xy = (0.01, 0.95), xycoords = 'axes fraction', 
                        xytext = (0.1, 0.95), ha = 'left', va = 'center', color = 'b',
                        arrowprops = dict(arrowstyle="->",color='b'))
          axis.annotate('%.2f' % au, xy = (0.99, 0.95), xycoords = 'axes fraction', 
                        xytext = (0.9, 0.95), ha = 'right', va = 'center', color = 'b',
                        arrowprops = dict(arrowstyle="->",color='b'))
        if len(control):  
          cu = len(np.where(control > ranges[i][1])[0]) / len(control)
          cl = len(np.where(control < ranges[i][0])[0]) / len(control)
          axis.annotate('%.2f' % cl, xy = (0.01, 0.88), xycoords = 'axes fraction', 
                        xytext = (0.1, 0.88), ha = 'left', va = 'center', color = 'r',
                        arrowprops = dict(arrowstyle="->",color='r'))
          axis.annotate('%.2f' % cu, xy = (0.99, 0.88), xycoords = 'axes fraction', 
                        xytext = (0.9, 0.88), ha = 'right', va = 'center', color = 'r',
                        arrowprops = dict(arrowstyle="->",color='r'))
                
        # Indicate the median
        if len(recovered):
          axis.annotate('M = %.2f' % np.median(recovered), xy = (0.3, 0.5), ha = 'right',
                        xycoords = 'axes fraction', color = 'b', fontsize = 14)
        if len(control):
          axis.annotate('M = %.2f' % np.median(control), xy = (0.7, 0.5), ha = 'left',
                        xycoords = 'axes fraction', color = 'r', fontsize = 14)
  
        # Tweaks
        axis.set_xticks(xticks[i])
        axis.set_xlim(xticks[i][0], xticks[i][-1])
        axis.set_ylim(0, ymax[i])
        axis.set_xlabel(r'D/D$_0$', fontsize = 14)

    pl.show()
    