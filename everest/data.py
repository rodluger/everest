#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`data.py` - Download routines
-------------------------------------

These are routines for downloading and storing the raw `K2` data, as well
as information about planet candidates and eclipsing binaries.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os
EVEREST_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
from .sources import GetSources, Source
from .utils import MedianFilter
import kplr
from kplr.config import KPLR_ROOT
import numpy as np
import re
import urllib
from tempfile import NamedTemporaryFile
import shutil
import logging
log = logging.getLogger(__name__)

float_ = lambda x: float(x) if x != '' else np.nan

adapter = {
            'rowid': str,
            'epic_name': str,
            'tm_name': str,
            'epic_candname': str,
            'pl_name': str,
            'k2c_refdisp': str,
            'k2c_reflink': str,
            'k2c_disp': str,
            'k2c_note': str,
            'k2_campaign': str,
            'k2c_recentflag': str,
            'ra_str': str,
            'ra': str,
            'dec_str': str,
            'dec': str,
            'pl_orbper': float_,
            'pl_orbpererr1': float_,
            'pl_orbpererr2': float_,
            'pl_orbperlim': float_,
            'pl_tranmid': float_,
            'pl_tranmiderr1': float_,
            'pl_tranmiderr2': float_,
            'pl_tranmidlim': float_,
            'pl_trandep': float_,
            'pl_trandeperr1': float_,
            'pl_trandeperr2': float_,
            'pl_trandeplim': float_,
            'pl_trandur': float_,
            'pl_trandurerr1': float_,
            'pl_trandurerr2': float_,
            'pl_trandurlim': float_,
            'pl_imppar': float_,
            'pl_impparerr1': float_,
            'pl_impparerr2': float_,
            'pl_impparlim': float_,
            'pl_orbincl': float_,
            'pl_orbinclerr1': float_,
            'pl_orbinclerr2': float_,
            'pl_orbincllim': float_,
            'pl_ratdor': float_,
            'pl_ratdorerr1': float_,
            'pl_ratdorerr2': float_,
            'pl_ratdorlim': float_,
            'pl_ratror': float_,
            'pl_ratrorerr1': float_,
            'pl_ratrorerr2': float_,
            'pl_ratrorlim': float_,
            'pl_rade': float_,
            'pl_radeerr1': float_,
            'pl_radeerr2': float_,
            'pl_radelim': float_,
            'pl_radj': float_,
            'pl_radjerr1': float_,
            'pl_radjerr2': float_,
            'pl_radjlim': float_,
            'pl_eqt': float_,
            'pl_eqterr1': float_,
            'pl_eqterr2': float_,
            'pl_eqtlim': float_,
            'pl_fppprob': float_,
            'pl_fppproblim': float_,
            'st_plx': float_,
            'st_plxerr1': float_,
            'st_plxerr2': float_,
            'st_plxlim': float_,
            'st_dist': float_,
            'st_disterr1': float_,
            'st_disterr2': float_,
            'st_distlim': float_,
            'st_teff': float_,
            'st_tefferr1': float_,
            'st_tefferr2': float_,
            'st_tefflim': float_,
            'st_logg': float_,
            'st_loggerr1': float_,
            'st_loggerr2': float_,
            'st_logglim': float_,
            'st_metfe': float_,
            'st_metfeerr1': float_,
            'st_metfeerr2': float_,
            'st_metfelim': float_,
            'st_metratio': str,
            'st_rad': float_,
            'st_raderr1': float_,
            'st_raderr2': float_,
            'st_radlim': float_,
            'st_vsini': float_,
            'st_vsinierr1': float_,
            'st_vsinierr2': float_,
            'st_vsinilim': float_,
            'st_kep': float_,
            'st_keperr': float_,
            'st_keplim': float_,
            'st_bj': float_,
            'st_bjerr': float_,
            'st_bjlim': float_,
            'st_vj': float_,
            'st_vjerr': float_,
            'st_vjlim': float_,
            'st_us': float_,
            'st_userr': float_,
            'st_uslim': float_,
            'st_gs': float_,
            'st_gserr': float_,
            'st_gslim': float_,
            'st_rs': float_,
            'st_rserr': float_,
            'st_rslim': float_,
            'st_is': float_,
            'st_iserr': float_,
            'st_islim': float_,
            'st_zs': float_,
            'st_zserr': float_,
            'st_zslim': float_,
            'st_j2': float_,
            'st_j2err': float_,
            'st_j2lim': float_,
            'st_h2': float_,
            'st_h2err': float_,
            'st_h2lim': float_,
            'st_k2': float_,
            'st_k2err': float_,
            'st_k2lim': float_,
            'st_wise1': float_,
            'st_wise1err': float_,
            'st_wise1lim': float_,
            'st_wise2': float_,
            'st_wise2err': float_,
            'st_wise2lim': float_,
            'st_wise3': float_,
            'st_wise3err': float_,
            'st_wise3lim': float_,
            'st_wise4': float_,
            'st_wise4err': float_,
            'st_wise4lim': float_,
            'st_bmvj': float_,
            'st_bmvjerr': float_,
            'st_bmvjlim': float_,
            'st_jmh2': float_,
            'st_jmh2err': float_,
            'st_jmh2lim': float_,
            'st_hmk2': float_,
            'st_hmk2err': float_,
            'st_hmk2lim': float_,
            'st_jmk2': float_,
            'st_jmk2err': float_,
            'st_jmk2lim': float_
}

class k2data(object):
  '''
  A generic `K2` data container. Nothing fancy here.
  
  '''
  
  pass

def GetK2Data(EPIC, apnum = 15, delete_kplr_data = True):
  '''
  Download and save a single quarter of `K2` data.
  
  :param int EPIC: The 9-digit `EPIC` number of the target
  
  :param apnum: The number of the aperture in the `K2SFF <https://archive.stsci.edu/prepds/k2sff/>`_ \
                fits file to use for the photometry. Default `15`
  :type apnum: int
  
  :param delete_kplr_data: Delete the fits file downloaded with :py:mod:`kplr` \
                           after processing it? Default `True`
  :type delete_kplr_data: bool
  
  :returns: 
    A :class:`k2data` object containing the following attributes:
  
    - **campaign** - The `K2` campaign the target was observed in
    - **time** - The array of timestamps, in `BJD - 245833`
    - **cadn** - The long cadence number corresponding to each observation
    - **fpix** - A 3-dimensional array of shape `(nt, nx, ny)` containing the \
                 raw flux in the pixel at position `(x, y)` at each timestamp `t`
    - **perr** - The standard error on each of the data points in `fpix`
    - **apertures** - An array containing the 20 aperture images obtained from the \
                      `K2SFF <https://archive.stsci.edu/prepds/k2sff/>`_ fits files
    - **aperture** - *Deprecated*
    - **bkg** - An estimate of the background flux at each cadence
    - **bkgerr** - The standard error on `bkg`
    - **kepmag** - The `Kepler` magnitude of the target
    - **planets** - A list of :class:`K2Planets` objects containing known planets or \
                    planet candidates for this target
    - **EB** - `False` if target is not an eclipsing binary; otherwise, a :class:`K2EB` \
               object containing EB info taken from the `Villanova <http://keplerebs.villanova.edu/>`_ \
               eclipsing binary catalog
    - **nearby** - A list of :class:`everest.sources.Source` instances containing \
                   other `EPIC` targets within or close to this target's aperture
    
  '''
  
  filename = os.path.join(KPLR_ROOT, 'data', 'everest', str(EPIC), str(EPIC) + '.npz')
  
  try:
    # Grab the K2SFF info, mainly to get the apertures
    k2sff = kplr.K2SFF(EPIC)
  except:
    # If we can't get the K2SFF files, we can't run Everest (for now)
    return None

  try:
    data = np.load(filename)
    time = data['time']
    fpix = data['fpix']
    perr = data['perr']
    campaign = data['campaign']
    aperture = data['aperture']
    cadn = data['cadn']
    _nearby = data['nearby']
    nearby = [Source(**s) for s in _nearby]
    clobber = False
  except:
    clobber = True
      
  if clobber:
    if not os.path.exists(os.path.join(KPLR_ROOT, 'data', 'everest', str(EPIC))):
      os.makedirs(os.path.join(KPLR_ROOT, 'data', 'everest', str(EPIC)))
  
    client = kplr.API()
    star = client.k2_star(EPIC)
    tpf = star.get_target_pixel_files()[0]
    campaign = tpf.sci_campaign
    with tpf.open() as f:
      aperture = f[2].data
      qdata = f[1].data
        
    # Get the arrays
    time = np.array(qdata.field('TIME'), dtype='float64')
    cadn = np.array(qdata.field('CADENCENO'), dtype='int32')
    fpix = np.array(qdata.field('FLUX'), dtype='float64')
    fpix_opt = np.array([f[np.where(aperture & 1)] for f in fpix], dtype='float64')
    rawc = np.array(qdata.field('RAW_CNTS'), dtype='int32')
    perr = np.array(qdata.field('FLUX_ERR'), dtype='float64')
    qual = np.array(qdata.field('QUALITY'), dtype=int)
    colmotion = np.array(qdata.field('POS_CORR1'), dtype='float64')
    rowmotion = np.array(qdata.field('POS_CORR2'), dtype='float64')

    # Get bad timestamps
    t_nan_inds = list(np.where(np.isnan(time))[0])
    
    # Get bad flux values
    apidx = np.where(k2sff.apertures[apnum] & 1 & ~np.isnan(fpix[0]))
    flux = np.sum(np.array([f[apidx] for f in fpix], dtype='float64'), axis = 1)
    f_nan_inds = list(np.where(np.isnan(flux))[0])
    
    # Get flagged data points. Note that like K2SFF, we do not throw out all data
    # points flagged with bit #15, but treat them separately. See "LDE Flags" in 
    # http://keplerscience.arc.nasa.gov/k2-data-release-notes.html#k2-campaign-2
    bad_bits = [1,2,3,4,5,6,7,8,9,11,12,13,14,16,17]                          
    qual_inds = []
    for b in bad_bits:
      qual_inds += list(np.where(qual & 2 ** (b - 1))[0])
    
    # Treat bit #15 separately. Only remove a data point flagged with bit 15
    # if it's more than 20% away from the local median.
    med = MedianFilter(flux, 10)
    dev = np.abs((flux - med) / med)
    bit15_inds = list(set(list(np.where(dev > 0.2)[0])) & \
                      set(list(np.where(qual & 2 ** (15 - 1))[0])))

    # All the bad inds
    bad_inds = np.array(sorted(list(set(qual_inds + t_nan_inds + f_nan_inds + bit15_inds))))

    # Campaign 0 hack. The first half of campaign zero was not in fine
    # pointing. Everest actually does pretty well at detrending it, but I'm 
    # actually finding a **time** offset for data before t ~ 1940. For EPIC
    # 202072596 (an EB), the primary eclipses prior to this time are offset
    # by a few cadences -- they happen later than they should. No idea why
    # this is happening, but I get *much* better folded eclipses when I remove
    # the first part of C0. Here I simply remove whatever SFF removes.
    if time[0] < 1940.:
      bad_inds = np.append(bad_inds, np.where(time < k2sff.time[0]))

    # Campaign 2 hack. The first 1-2 days in C2 have very different noise
    # properties than the rest of the campaign, so we'll again trust the SFF
    # cuts.
    if time[0] < 2061 and time[0] > 2059:
      bad_inds = np.append(bad_inds, np.where(time < k2sff.time[0]))

    # Remove them
    time = np.delete(time, bad_inds)
    cadn = np.delete(cadn, bad_inds)
    fpix = np.delete(fpix, bad_inds, 0)
    rawc = np.delete(rawc, bad_inds, 0)
    perr = np.delete(perr, bad_inds, 0)
    colmotion = np.delete(colmotion, bad_inds)
    rowmotion = np.delete(rowmotion, bad_inds)
    
    # Get nearby targets
    _nearby = [s.__dict__ for s in GetSources(EPIC)]
  
    # Save
    np.savez_compressed(filename, time = time, fpix = fpix, perr = perr, cadn = cadn,
                        aperture = aperture, nearby = _nearby, campaign = campaign)
  
    # Make it into an object
    nearby = [Source(**s) for s in _nearby]
  
    # Delete the kplr tpf
    if delete_kplr_data:
      os.remove(os.path.join(KPLR_ROOT, 'data', 'k2', 'target_pixel_files', '%d' % EPIC, tpf._filename))
  
  # Get any K2 planets associated with this EPIC
  planets = []
  for planet in GetK2Planets():
    if planet.epic_name == 'EPIC %d' % EPIC:
      if planet.epic_candname not in [p.epic_candname for p in planets]:
        planets.append(planet)
  if len(planets):
    planets = sorted(planets, key = lambda x: x.epic_candname)
      
  # Create the object to return
  res = k2data()
  res.campaign = campaign
  res.time = time
  res.cadn = cadn
  res.fpix = fpix
  res.apertures = k2sff.apertures
  
  # Compute the background from the median outside the aperture
  binds = np.where(res.apertures[apnum] ^ 1)
  if len(binds[0]) > 0:
    res.bkg = np.nanmedian(np.array([f[binds] for f in fpix], dtype='float64'), axis = 1)
    res.bkgerr = np.nanmedian(np.array([p[binds] for p in perr], dtype='float64'), axis = 1)
  else:
    # Unable to determine the background!
    res.bkg = np.zeros_like(time)
    res.bkgerr = np.zeros_like(time)
    log.warn('Unable to compute the background flux for EPIC %d.' % EPIC)
    log.warn('Consider re-running this target with a smaller aperture.')
    
  res.perr = perr
  res.aperture = aperture
  res.planets = planets
  
  # Is this an EB?
  res.EB = False
  for eb in GetK2EBs():
    if eb.epic == EPIC:
      res.EB = eb
      continue
  
  # Get the kepler magnitude from the sources in ``nearby``
  # Sometimes the magnitude is not listed in the TPF, but
  # it's in the MAST database...
  foo = [s.kepmag for s in nearby if s.epic == EPIC]
  if len(foo):
    res.kepmag = foo[0]
  else:
    res.kepmag = np.nan
  
  # Get nearby sources
  res._nearby = _nearby
  res.nearby = nearby

  return res

class K2Planet(object):
  '''
  A generic `K2` planet candidate object.
  
  '''
  
  def __init__(self):
    self.epic_candname = ""
  
  def __repr__(self):
    return "<K2Planet: %s>" % self.epic_candname

class K2EB(object):
  '''
  A generic `K2` EB candidate object.
  
  '''
  
  def __init__(self):
    self.epic = ""
  
  def __repr__(self):
    return "<K2EB: %s>" % self.epic

def GetK2Stars(clobber = False):
  '''
  Download and return a `dict` of all `K2` stars organized by campaign. Saves each
  campaign to a `csv` file in the `/tables` directory.
  
  :param bool clobber: If `True`, download and overwrite existing files. Default `False`
  
  '''
  
  # Download
  if clobber:
    client = kplr.API()
    stars = client.k2_star_list()
    for campaign in stars.keys():
      with open(os.path.join(EVEREST_ROOT, 'tables', 'C%02d.csv' % campaign), 'w') as f:
        for star in stars[campaign]:
          print(star, file = f)
  
  # Return
  res = {}
  for campaign in range(100):
    f = os.path.join(EVEREST_ROOT, 'tables', 'C%02d.csv' % campaign)
    if os.path.exists(f):
      stars = np.loadtxt(f, dtype = int)
      res.update({campaign: stars})
  
  return res

def GetK2InjectionTestStars(clobber = False):
  '''
  Download and return a dict of 2000 `K2` stars, with 100 stars per magnitude 
  bin in the range 8-18. These are used for injection tests. The stars are
  saved in `tables/Injections.csv`.
  
  '''
  
  # Download
  if clobber:
    client = kplr.API()
    allstars = client.k2_star_mags(stars_per_mag = 200, mags = range(8,18))
    with open(os.path.join(EVEREST_ROOT, 'tables', 'Injections.csv'), 'w') as f:
      for stars in allstars: print(", ".join([str(s) for s in stars]), file = f)
  
  # Return the flattened list
  stars = np.loadtxt(os.path.join(EVEREST_ROOT, 'tables', 'Injections.csv'), 
                     dtype = int, delimiter = ',')
  return [item for sublist in stars for item in sublist]
  
def GetK2Planets():
  '''
  Returns a list of :class:`K2Planet` instances generated from the file
  `/tables/k2candidates.csv`. This file was downloaded from the
  `Exoplanet Archive <http://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=k2candidates>`_
  on February 26, 2016.
  
  '''
  
  # Read the CSV file
  with open(os.path.join(EVEREST_ROOT, 'tables', 'k2candidates.csv'), 'r') as f:
    lines = f.readlines()

  # Get columns
  columns = [('rowid', 'Row ID')]
  for line in lines:
    regex = re.match("# COLUMN (.*):[ X]+(.*)", line)
    if regex is not None:
      columns.append((regex.groups()[0], regex.groups()[1]))

  # Get the data
  planets = []
  for line in lines[len(columns) + 4:]:
    planet = K2Planet()
    line = line.replace('\n', '')
    entries = line.split(',')
    for c, e in zip(columns, entries):
      setattr(planet, c[0], adapter[c[0]](e))    
    planets.append(planet)

  return planets

def VillanovaBJDOffset(campaign):
  '''
  There's a strange time offset in the `Villanova` EB catalog for 
  some campaigns. This function returns the time offset for a given
  campaign, which was determined empirically. These numbers have not
  been thoroughly verified.
  
  '''
  
  if campaign == 0:
    return -54833.            # For 202060523, the offset is -54833.279. (!?!?)
  elif campaign == 1:
    return 0.                 # For 201158453, the offset is -1.94. (!?!?)
  elif campaign == 2:
    return -54833.
  elif campaign >= 3:
    return -2454833.

class EclipseTimes(object):
  '''
  A simple class that determines the times of all eclipses for a given EB.
  
  :param float t0: The time of first eclipse
  :param float period: The period in days
  :param float duration: The eclipse duration in days
  
  :returns: A :py:class:`numpy` array of the times of all transits between `start` \
            and `stop`
  
  '''
  
  def __init__(self, t0, period, duration):
    '''
    
    '''
    
    self.t0 = t0
    self.period = period
    self.duration = duration
    
  def __call__(self, start, end):
    '''

    '''
    if self.duration > 0:
      return np.arange(self.t0 + np.ceil((start - self.duration - self.t0) / self.period) 
                       * self.period, end + self.duration, self.period)
    else:
      return np.array([], dtype = float)

class EclipseMask(object):
  '''
  An eclipse masking object for EBs.
  
  :param `EclipseTimes` primary: An instance containing the times of primary eclipse
  :param `EclipseTimes` secondary: An instance containing the times of secondary eclipse
  
  :returns: The indices in `time` that contain the primary and secondary eclipses
  
  '''
  
  def __init__(self, primary, secondary):
    self.primary = primary
    self.secondary = secondary
    
  def __call__(self, time):
    '''
  
    '''
    
    p = []; s = []
    for t in self.primary(time[0], time[-1]): 
      p.extend(np.where(np.abs(time - t) < self.primary.duration / 2.)[0])
    for t in self.secondary(time[0], time[-1]): 
      s.extend(np.where(np.abs(time - t) < self.secondary.duration / 2.)[0])
    
    return sorted(set(p + s))
                   
def GetK2EBs(clobber = False):
  '''
  Grab all `K2` EBs from the pre-downloaded `Villanova` catalog, which is stored in
  `/tables/k2ebs.tsv`.
  
  :param bool clobber: If `True`, download and overwrite existing files. Default `False`
  
  '''
  
  # Download a new CSV file?
  if clobber or not os.path.exists(os.path.join(EVEREST_ROOT, 'tables', 'k2ebs.tsv')):
    url = 'http://keplerebs.villanova.edu/results/?q={"sort":"kic",' + \
          '"campaign":["0","1","2","3","4","5","6","7","8","9"],' + \
          '"kics":[],"etvlong":true,' + \
          '"cols":["camp","bjd0","kic","p","sep","pwidth","swidth"],' + \
          '"etvshort":true,"incat1":true,"kois":[]}&format=csv'
    r = urllib.request.Request(url)
    handler = urllib.request.urlopen(r)
    code = handler.getcode()
    if int(code) != 200:
      raise Exception("Error downloading Villanova EB data.")
    data = handler.read()
    f = NamedTemporaryFile("wb", delete=False)
    f.write(data)
    f.flush()
    os.fsync(f.fileno())
    f.close()
    shutil.move(f.name, os.path.join(EVEREST_ROOT, 'tables', 'k2ebs.tsv'))
  
  # Read the CSV file
  with open(os.path.join(EVEREST_ROOT, 'tables', 'k2ebs.tsv'), 'r') as f:
    lines = f.readlines()

  # Create a list of EB objects
  EBs = []
  for line in lines:
    if line.startswith('#') or len(line) < 5:
      continue
    line = line.replace(',\n', '')
    epic, campaign, period, bjd0, pwidth, swidth, sep = line.split(',')
    
    # Create our EB instance
    EB = K2EB()
    EB.epic = int(epic)
    EB.campaign = int(campaign.replace('K2C', ''))
    EB.period = float(period)
    EB.p0 = float(bjd0) + VillanovaBJDOffset(EB.campaign)
    EB.s0 = EB.p0 + float(sep) * EB.period
    EB.pdur = float(pwidth) * 2 * EB.period
    EB.sdur = float(swidth) * 2 * EB.period
    
    # The Villanova catalog lists some separations as 1
    # and some secondary durations as -1, presumably
    # when a secondary isn't detected. We're going to assume
    # the secondary is at primary + per / 2. for plotting
    # purposes, but we won't mask it.
    if float(sep) == 1. or float(swidth) < 0:
      EB.s0 = EB.p0 + EB.period / 2.
      EB.sdur = 0.
    
    # Some of the primary/secondary durations are just
    # pure nonsense -- for EPIC 201160662, for instance, the
    # primary duration + secondary duration equals 1.2 times
    # the period, i.e., primary and secondary overlap a bit
    # and the star is in permanent eclipse (?!?!)
    if EB.pdur + EB.sdur > 0.5 * EB.period:
      # We're just not going to bother with these EBs for now.
      continue

    # Get the times of primary eclipse center and secondary eclipse center
    EB.primary = EclipseTimes(EB.p0, EB.period, EB.pdur)
    EB.secondary = EclipseTimes(EB.s0, EB.period, EB.sdur)
    
    # Get the indices during all eclipses (for masking)
    EB.mask = EclipseMask(EB.primary, EB.secondary)
      
    # Append to the list
    EBs.append(EB)
  
  # Now read the user-defined list of updated EBs
  with open(os.path.join(EVEREST_ROOT, 'tables', 'k2ebs_updated.tsv'), 'r') as f:
    lines = f.readlines()
  
  for line in lines:
    if line.startswith('#') or len(line) < 5:
      continue
    line = line.replace(',\n', '')
    epic, campaign, period, bjd0, pwidth, swidth, sep = line.split(',')
    
    # Create our EB instance
    EB = K2EB()
    EB.epic = int(epic)
    EB.campaign = int(campaign.replace('K2C', ''))
    EB.period = float(period)
    EB.p0 = float(bjd0) + VillanovaBJDOffset(EB.campaign)
    EB.s0 = EB.p0 + float(sep) * EB.period
    EB.pdur = float(pwidth) * 2 * EB.period
    EB.sdur = float(swidth) * 2 * EB.period
    
    # The Villanova catalog lists some separations as 1
    # and some secondary durations as -1, presumably
    # when a secondary isn't detected. We're going to assume
    # the secondary is at primary + per / 2. for plotting
    # purposes, but we won't mask it.
    if float(sep) == 1. or float(swidth) < 0:
      EB.s0 = EB.p0 + EB.period / 2.
      EB.sdur = 0.
    
    # Get the times of primary eclipse center and secondary eclipse center
    EB.primary = EclipseTimes(EB.p0, EB.period, EB.pdur)
    EB.secondary = EclipseTimes(EB.s0, EB.period, EB.sdur)
    
    # Get the indices during all eclipses (for masking)
    EB.mask = EclipseMask(EB.primary, EB.secondary)
      
    # Append to the list
    if int(epic) not in [eb.epic for eb in EBs]:
      EBs.append(EB)
    else:
      i = np.argmax(np.array([eb.epic for eb in EBs], dtype = int) == int(epic))
      EBs[i] = EB
      
  return EBs