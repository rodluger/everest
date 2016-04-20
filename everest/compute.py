#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
compute.py
----------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os
EVEREST_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
from .data import GetK2Data
from .detrend import PLDBasis, PLDModel, PLDCoeffs, ComputeScatter
from .utils import InitLog, Mask, GetMasks, Breakpoint, SatSev, \
                   AcorSev, PadWithZeros, RMS, Outliers
from .kernels import KernelModels
from .gp import GetGP
from .transit import Transit
import george
import numpy as np
from scipy.signal import savgol_filter
import subprocess
import re
import logging
log = logging.getLogger(__name__)

def Compute(EPIC, run_name = 'default', clobber = False, apnum = 15, 
            outlier_sigma = 5, mask_times = [],
            ps_iter = 100, npc_arr = np.arange(25, 250, 10),
            inject = {}, log_level = logging.DEBUG, scatter_alpha = 0.,
            screen_level = logging.DEBUG, gp_iter = 1, **kwargs):
  '''
  
  '''
  
  # Grab the data
  k2star = GetK2Data(EPIC, apnum = apnum)
  if k2star is None:
    return None

  # Setup the output directory
  outdir = os.path.join(EVEREST_ROOT, 'output', 'C%02d' % k2star.campaign, 
                        '%d' % EPIC, run_name)
  if not os.path.exists(outdir):
    os.makedirs(outdir)
  
  # Setup logging
  InitLog(os.path.join(outdir, 'compute.log'), log_level = log_level, 
          screen_level = screen_level)
  
  # Check if we've done this already
  if os.path.exists(os.path.join(outdir, 'data.npz')) and not clobber:
    # Load from disk
    data = dict(np.load(os.path.join(outdir, 'data.npz')))
    # Replace the ``gp`` entry with an actual GP instance
    gpinfo = data['gpinfo'][()]
    kpars = gpinfo['kpars']
    knum = gpinfo['knum']
    kernel = KernelModels[knum]
    kernel[:] = kpars
    data['gp'] = george.GP(kernel.george_kernel())
    return data
  
  # Begin pre-processing
  log.info('Pre-processing the data...')
  time = k2star.time
  bkg = k2star.bkg
  bkgerr = k2star.bkgerr

  # Calculate the total background-corrected flux in the aperture
  # Note that some K2SFF apertures include pixels that are always NaN, so we
  # remove those here.
  apidx = np.where(k2star.apertures[apnum] & 1 & ~np.isnan(k2star.fpix[0]))
  fpix = np.array([f[apidx] for f in k2star.fpix], dtype='float64') - \
                   bkg.reshape(bkg.shape[0], 1)
  perr = np.array([f[apidx] for f in k2star.perr], dtype='float64')
  perr = np.sqrt(perr ** 2 + bkgerr.reshape(bkgerr.shape[0], 1) ** 2)
  flux = np.sum(fpix, axis = 1)
  ferr = np.sqrt(np.sum(perr ** 2, axis = 1))
  npix_total = fpix.shape[1]
  
  # 3rd order PLD on 75 pixels for a typical K2 campaign is over 2 GB
  # of arrays. When we parallelize this, we can run into memory errors
  # on the cluster. Also, stars with 60+ pixels are almost certainly
  # saturated.
  if fpix.shape[1] > 75:
    log.info('Too many pixels in aperture! Aborting to avoid memory errors.')
    with open(os.path.join(outdir, '%d.err' % EPIC), 'w') as f:
      print('ERROR: Too many pixels in aperture (%d).' % fpix.shape[1])
    return None
  
  # Inject transits?
  if len(inject):
  
    # Set missing params to default values
    inject.update({'t0': inject.get('t0', 0.)})
    inject.update({'per': inject.get('per', 3.56789)})
    inject.update({'dur': inject.get('dur', 0.1)})
    inject.update({'depth': inject.get('depth', 0.001)})
    inject.update({'mask': inject.get('mask', False)})
    
    # Remove all outliers. They screw up our linear transit fitting routine.
    # This won't actually skew our recovery results -- it's just simpler than
    # doing some sort of iterative outlier removal step in the recovery.
    log.info('Removing outliers...')
    outliers, _ = Outliers(time, flux, fpix = fpix, ferr = ferr, sigma = 3)
    time = np.delete(time, outliers)
    flux = np.delete(flux, outliers)
    ferr = np.delete(ferr, outliers)
    bkg = np.delete(bkg, outliers)
    bkgerr = np.delete(bkgerr, outliers)
    fpix = np.delete(fpix, outliers, axis = 0)
    perr = np.delete(perr, outliers, axis = 0)
        
    # Compute the transit model and inject
    log.info('Injecting transits...')
    transit_model = Transit(time, **inject)
    if inject['mask']:
      mask_times = sorted(set(mask_times) | set(time[np.where(transit_model < 1.)]))
    for i in range(fpix.shape[1]):
      fpix[:,i] *= transit_model    
    flux = np.sum(fpix, axis = 1)
  
  # Are there any new planet candidates in this lightcurve?
  # NOTE: This section was coded specifically for Ethan Kruse's 
  # transit search pipeline.
  new_candidates = []
  if kwargs.get('new_candidate', False):
    cnum = 1
    while os.path.exists(os.path.join(EVEREST_ROOT, 'newcandidates', 
                         'KIC%d.QATSinds%d' % (EPIC, cnum))):
      log.info('Masking new planet candidate #%d...' % (cnum))
      filename = os.path.join(EVEREST_ROOT, 'newcandidates', 
                              'KIC%d.QATSinds%d' % (EPIC, cnum))
      with open(filename, 'r') as f:
        # Ethan's duration is in hours, so convert to days. Also add a bit of a buffer
        new_tdur = float(f.readlines()[1][2:7]) * 1.2 / 24.
      _, _, _, new_ttimes, _ = np.loadtxt(filename).T
      # Ethan uses a different epoch; convert
      new_ttimes += (2455000 - 2454833)
      # Update our mask
      new_tmask = np.concatenate([time[np.where(np.abs(time - tn) < new_tdur)] 
                                  for tn in new_ttimes])
      mask_times = sorted(set(mask_times) | set(new_tmask))
      new_candidates.append([{'tdur': new_tdur, 'ttimes': new_ttimes, 
                              'tmask': new_tmask}])
      cnum += 1
  
  # Obtain transit and outlier masks
  log.info('Computing transit and outlier masks...')
  try:
    maskdata = np.load(os.path.join(outdir, 'mask.npz'))
    mask = maskdata['mask']
    trn_mask = maskdata['trn_mask']
    rem_mask = maskdata['rem_mask']
    keep_mask = maskdata['keep_mask']
  except:
    mask, trn_mask, rem_mask, keep_mask = \
    GetMasks(time, flux, fpix, ferr, outlier_sigma, planets = k2star.planets, 
             EB = k2star.EB, mask_times = mask_times)
    np.savez(os.path.join(outdir, 'mask.npz'), mask = mask,
             trn_mask = trn_mask, rem_mask = rem_mask, keep_mask = keep_mask)

  # Compute GP hyperparameters
  log.info('Computing the GP hyperparameters...')
  try:
    gpdata = np.load(os.path.join(outdir, 'gpdata.npz'))
    knum = gpdata['knum']
    kpars = gpdata['kpars']
    white = gpdata['white']
    amp = gpdata['amp']
    chisq = gpdata['chisq']
    powerspec = gpdata['powerspec']
    acor = gpdata['acor']
    kernfunc = gpdata['kernfunc']
    fpldgp = gpdata['fpldgp']
  except:
    res = GetGP(EPIC, time, fpix, ferr, mask = mask, niter = gp_iter)
    knum, kpars, white, amp, chisq, fpldgp = res['data']
    powerspec = res['powerspec']
    acor = res['acor']
    kernfunc = res['kernfunc']
    np.savez(os.path.join(outdir, 'gpdata.npz'), knum = knum, kpars = kpars, 
             white = white, amp = amp, chisq = chisq, powerspec = powerspec, 
             acor = acor, kernfunc = kernfunc, fpldgp = fpldgp)

  # Set up the Gaussian process
  kernel = KernelModels[knum]
  kernel[:] = kpars
  gp = george.GP(kernel.george_kernel())

  # Optimize the PLD order and number of PCA components
  # ``npc_pred`` is the full array of # of components
  npc_pred = np.arange(npc_arr[0], npc_arr[-1])
  try:
    best = np.load(os.path.join(outdir, 'best.npz'))
    pld_order = best['pld_order']
    npc = best['npc']
    bestij = best['bestij']
    X = best['X']
    masked_scatter = best['masked_scatter']
    unmasked_scatter = best['unmasked_scatter']
    msf = best['msf']
    usf = best['usf']
  except:
    # Compute PLD basis vectors and save them to disk
    log.info('Computing the basis vectors...')
    breakpoint = Breakpoint(k2star.campaign, time, mask)
    if breakpoint is None:
      breakpoints = []
    else:
      breakpoints = [breakpoint]
    if not os.path.exists(os.path.join(outdir, 'best.npz')):
      for pld_order in [1,2,3]:
        file = os.path.join(outdir, 'X_%02d.npz' % pld_order)
        if not os.path.exists(file):
          X = PLDBasis(fpix, time = time, pld_order = pld_order, 
                       max_components = npc_arr[-1], 
                       breakpoints = breakpoints)
          np.savez(file, X = X)
          del X
  
    # Reload them from disk
    tmp = [np.load(os.path.join(outdir, 'X_%02d.npz' % pld_order))['X'] 
                   for pld_order in [1,2,3]]
    X_arr = np.empty(len(tmp), dtype = object)
    X_arr[:] = tmp
  
    # Compute the scatter for each combination of (pld order, number of components)
    log.info('Minimizing the predictive scatter...')
    masked_scatter = np.zeros((3, len(npc_arr)))
    unmasked_scatter = np.zeros((3, len(npc_arr)))
    for i, pld_order in enumerate([1,2,3]):
      for j, npc in enumerate(npc_arr):
        log.debug('PLD order = %d, Number of components = %d...' % (pld_order, npc))
        masked_scatter[i,j], unmasked_scatter[i,j] = \
        ComputeScatter(X_arr[i][:,:npc], flux, time, ferr, 
                       gp, mask = mask, niter = ps_iter)
    
    # Find the params that minimize the scatter  
    i, j, msf, usf = MinimizeScatter(npc_arr, npc_pred, 
                     masked_scatter, unmasked_scatter, alpha = scatter_alpha)
    bestij = [i, j]
    pld_order = i + 1
    npc = npc_pred[j]
    X = X_arr[i][:,:npc]
    
    # Save and delete the large X files
    np.savez(os.path.join(outdir, 'best.npz'), X = X, pld_order = pld_order, 
             npc = npc, bestij = bestij, masked_scatter = masked_scatter,
             unmasked_scatter = unmasked_scatter, msf = msf, usf = usf)
    for pld_order, nchunk in bestij:
      os.remove(os.path.join(outdir, 'X_%02d.npz' % (pld_order)))

  # Now detrend with the best values
  log.info('Detrending using best solution...') 

  # Compute the PLD model and the PLD-de-trended flux
  M = Mask(mask)
  gp.compute(M(time), M(ferr))
  A = np.dot(M(X).T, gp.solver.apply_inverse(M(X)))
  B = np.dot(M(X).T, gp.solver.apply_inverse(M(flux)))
  C = np.linalg.solve(A, B)
  model = np.dot(C, X.T)
  fpld = flux - model
  fpld += np.median(flux)
  fpld_norm = fpld / np.median(M(fpld))  
  
  # Now compute the fully whitened flux
  fwhite = flux - model
  fwhite += np.median(flux)
  med = np.median(M(fwhite))
  outliers, _ = Outliers(M(time), M(fwhite), sigma = 5)
  O = Mask(outliers)
  gp.compute(O(M(time)), O(M(ferr)))
  mu, _ = gp.predict(O(M(fwhite)) - med, time)
  fwhite = (fwhite - mu) 
  fwhite_norm = fwhite / med

  # If we injected a transit, let's try to recover the depth
  if len(inject):
    inject.update({'everest': GetTransitDepth(time, fpld_norm, inject, 
                              buf = 5, order = 2)})
    # Save recovery info to disk
    with open(os.path.join(outdir, '%d.inj' % EPIC), 'w') as injfile:
      if inject['mask']:
        print('Masked:            True', file = injfile)
      else:
        print('Masked:            False', file = injfile)
      print('Injected period:   %.2f d' % inject['per'], file = injfile)
      print('Injected duration: %.2f d' % inject['dur'], file = injfile)
      print('Injected depth:    %.5f' % inject['depth'], file = injfile)
      print('Everest depth:     %.8f +/- %.8f' % (inject['everest'][0], 
                                inject['everest'][1]), file = injfile)

  # Get the CDPP of the raw and de-trended data
  rms_raw_simple, rms_raw_savgol, rms_evr_simple, \
  rms_evr_savgol, rms_pht = GetCDPP(flux, fpld)
  
  # Gauge the saturation, crowding, and acor severity
  crwdinfo, crwdsev = GetContaminants(EPIC, k2star.fpix, k2star.apertures, 
                                      apnum, k2star.kepmag, k2star.nearby)
  satsev = SatSev(fpix)
  acorsev = AcorSev(chisq)

  # Write the .pld file to disk
  log.info('Saving detrended lightcurve to disk...') 
  WritePLDFile(EPIC, k2star.kepmag, satsev, crwdsev, crwdinfo, chisq, 
               rms_raw_simple, rms_raw_savgol, rms_evr_simple, rms_evr_savgol, 
               rms_pht, time, fpld, np.sqrt(ferr ** 2 + white ** 2), outdir)

  # Get the git info (if we're in a git repo)
  try:
    branches = subprocess.check_output(['git', 
               'branch']).decode('utf-8').replace('\n', '')
    git_branch = re.findall('\*\s([a-zA-Z0-9_]*)', branches)[0]
    git_hash = subprocess.check_output(['git', 'rev-parse', 
               '--verify', 'HEAD']).decode('utf-8').replace('\n', '')
  except:
    git_branch = None
    git_hash = None
  
  # Save everything to one massive .npz file
  data = dict(time = time, flux = flux, ferr = ferr, fpix = fpix, perr = perr, 
             fpix_full = k2star.fpix, apertures = k2star.apertures, mask = mask,
             trn_mask = trn_mask, rem_mask = rem_mask, keep_mask = keep_mask, 
             apnum = apnum, bkg = bkg, inject = inject, new_candidates = new_candidates, 
             kepmag = k2star.kepmag, EB = k2star.EB,
             planets = [planet.__dict__ for planet in k2star.planets], 
             nearby = k2star._nearby, fpld = fpld_norm, fwhite = fwhite_norm, 
             C = C, X = X, gp = None, gpinfo = dict(knum = knum, kpars = kpars, 
             white = white), fpldgp = fpldgp, pld_order = pld_order,
             rms = [rms_raw_simple, rms_raw_savgol, rms_evr_simple, 
             rms_evr_savgol, rms_pht], satsev = satsev, crwdsev = crwdsev, 
             acorsev = acorsev, chisq = chisq, npc_arr = npc_arr, npc_pred = npc_pred,
             masked_scatter = masked_scatter, unmasked_scatter = unmasked_scatter, 
             msf = msf, usf = usf, bestij = bestij,
             acor = acor, powerspec = powerspec, white = white, amp = amp, 
             kernfunc = kernfunc, EPIC = EPIC, run_name = run_name, 
             git_hash = git_hash, git_branch = git_branch, outdir = outdir,
             campaign = k2star.campaign)
  np.savez_compressed(os.path.join(outdir, 'data.npz'), **data)
  
  # Finally, delete the old .npz files
  for file in ['best.npz', 'gpdata.npz', 'mask.npz']:
    os.remove(os.path.join(outdir, file))
  
  # Replace the ``gp`` entry with the actual GP instance and return
  data['gp'] = gp
  return data

def WritePLDFile(EPIC, kepmag, satsev, crwdsev, crwdinfo, kchisq, r1, r2, r3, r4, r5,
                 time, fpld, ferr, outdir):
  '''
  
  '''
  
  # Info for the output file headers
  header = [
   "#",
   "# EPIC %d" % EPIC,
   "# Kp:  %.3f" % kepmag,
   "#",
   "# SATURATION FLAG:    %d" % satsev,
   "# CROWDING SEVERITY:  %d" % crwdsev,
   "# ACOR CHISQ:         %.2f" % kchisq,
   "#",
   "# RAW PRECISION:    {:8.2f} ppm / {:8.2f} ppm".format(r1, r2),
   "# ANCHOR PRECISION: {:8.2f} ppm / {:8.2f} ppm".format(r3, r4),
   "# PHOTON LIMIT:     {:8.2f} ppm".format(r5),
   "#"
  ]

  with open(os.path.join(outdir, '%d.pld' % EPIC), 'w') as pldfile:
    for line in header:
      print(line, file = pldfile)
    print(crwdinfo, file = pldfile)
    for t, f, s in zip(time, fpld, ferr):
      print('%.10e\t%.10e\t%.10e' % (t, f, s), file = pldfile)  

def GetCDPP(flux, fpld):
  '''
  Returns the 6-hr CDPP (RMS) of the raw and de-trended data, with and without
  a Savitsky-Golay filter applied, as well as the approximate photon limit
  
  '''
  
  # 1. Raw flux
  rms_raw_simple = RMS(flux / np.nanmedian(flux), remove_outliers = True)
  flux_savgol = flux - savgol_filter(flux, 49, 2) + np.nanmedian(flux)
  rms_raw_savgol = RMS(flux_savgol / np.nanmedian(flux_savgol), remove_outliers = True)
  
  # 2. Everest flux
  rms_evr_simple = RMS(fpld / np.nanmedian(fpld), remove_outliers = True)
  fpld_savgol = fpld - savgol_filter(fpld, 49, 2) + np.nanmedian(fpld)
  rms_evr_savgol = RMS(fpld_savgol / np.nanmedian(fpld_savgol), remove_outliers = True)
  
  # 3. Photon limit
  rms_pht = 1.e6 / np.sqrt(np.mean(flux) * 21600)
  
  return rms_raw_simple, rms_raw_savgol, rms_evr_simple, rms_evr_savgol, rms_pht
  
def MinimizeScatter(npc_arr, npc_pred, masked_scatter, unmasked_scatter, alpha = 0.):
  '''
  
  '''

  # Some global minima/maxima, etc.
  ytop = np.inf
  mms = [np.inf, (0, 0)]
  
  # Masked and unmasked "scatter functions": smooth GP-generated curves
  # that approximate the scatter as a function of the number of components
  msf = np.zeros((3, len(npc_pred))) * np.nan
  usf = np.zeros((3, len(npc_pred))) * np.nan
  
  for i in range(3):
      
    # Fit the scatter with a GP. GP params are hard-coded for now.
    sig = 1.4826 * np.nanmedian(np.abs(masked_scatter[i] - np.nanmedian(masked_scatter[i])))
    amp = 10 * sig
    tau = 100
    gp_scatter = george.GP(amp ** 2 * george.kernels.Matern32Kernel(tau ** 2))
    try:
      gp_scatter.compute(npc_arr, np.ones_like(npc_arr) * sig)
    except np.linalg.linalg.LinAlgError:
      continue
    
    # NOTE: I'm wrapping this in an exception catcher since sometimes there could
    # be NaNs in ``masked_scatter`` and ``unmasked_scatter``. Not sure how george 
    # handles those, so let's play it safe... Eventually I should treat this more 
    # rigorously.
    try:
      # Predicted scatter in masked regions
      msf[i], _ = gp_scatter.predict(masked_scatter[i] - np.nanmedian(masked_scatter[i]), npc_pred)
      msf[i] += np.median(masked_scatter[i])
    
      # Computed scatter in unmasked regions
      usf[i], _ = gp_scatter.predict(unmasked_scatter[i] - np.nanmedian(unmasked_scatter[i]), npc_pred)
      usf[i] += np.median(unmasked_scatter[i])
    except:
      continue
       
    # The minimum masked scatter (mms). Sometimes it might be desirable to
    # also minimize the *difference* between the masked and unmasked scatter
    # to avoid choosing outliers on the scatter plot. In that case, set alpha
    # to something > 0 (but <= 1).
    met = msf[i] + alpha * np.abs(msf[i] - usf[i])
    j = np.nanargmin(met)
    if met[j] < mms[0]:
      mms = [met[j], (i, j)]

  # Finally, get the best run
  i, j = mms[1]

  return i, j, msf, usf

def GetTransitDepth(time, flux, inject, buf = 5, order = 3):
  '''
  Recover the injected transit depth.
  
  '''
  
  t0 = inject['t0']
  per = inject['per']
  dur = inject['dur']
  transit_model = (Transit(time, **inject) - 1) / inject['depth']
  
  # Count the transits
  t0 += np.ceil((time[0] - dur - t0) / per) * per
  ttimes0 = np.arange(t0, time[-1] + dur, per)
  tinds = []
  for tt in ttimes0:
    # Get indices for this chunk
    inds = np.where(np.abs(time - tt) < buf * dur / 2.)[0]
    # Ensure there's a transit in this chunk, and that
    # there are enough points for the polynomial fit
    if np.any(transit_model[inds] < 0.) and len(inds) > order:
      tinds.append(inds)

  # Our design matrix
  sz = (order + 1) * len(tinds)
  X = np.empty([0, 1 + sz], dtype = float)
  Y = np.array([], dtype = float)
  T = np.array([], dtype = float)
  
  # Loop over all transits
  for i, inds in enumerate(tinds):
    # Get the transit model
    trnvec = transit_model[inds].reshape(-1, 1)
    # Normalize the time array
    t = time[inds]
    t = (t - t[0]) / (t[-1] - t[0])
    # Cumulative arrays
    T = np.append(T, time[inds])
    Y = np.append(Y, flux[inds] / np.median(flux[inds]))
    # Polynomial vector
    polyvec = np.array([t ** o for o in range(0, order + 1)]).T
    # Update the design matrix with this chunk
    lzeros = np.zeros((len(t), i * (order + 1)))
    rzeros = np.zeros((len(t), sz - (i + 1) * (order + 1)))
    chunk = np.hstack((trnvec, lzeros, polyvec, rzeros))
    X = np.vstack((X, chunk))

  # Get the relative depth  
  A = np.dot(X.T, X)
  B = np.dot(X.T, Y)
  C = np.linalg.solve(A, B)
  depth = C[0]
  
  # Get the uncertainties
  sig = 1.4826 * np.median(np.abs(flux - np.median(flux))) / np.median(flux)
  cov = sig ** 2 * np.linalg.solve(A, np.eye(A.shape[0]))
  err = np.sqrt(np.diag(cov))
  depth_err = err[0]
  
  # Get the detrended, folded data
  D = (Y - np.dot(C[1:], X[:,1:].T) + np.median(Y)) / np.median(Y)
  fold = lambda t: (t - inject['t0'] - inject['per'] / 2.) % inject['per'] - inject['per'] / 2.  
  
  return depth, depth_err, fold(T), D

def GetContaminants(EPIC, fpix, apertures, apnum, kepmag, nearby):
  '''
  Returns nearby sources that are at risk of contaminating a given target.
  
  '''
  
  _, ny, nx = fpix.shape
  apidx = np.where(apertures[apnum] & 1)
  contour = np.zeros((ny,nx))
  contour[apidx] = 1
  contour = np.lib.pad(contour, 1, PadWithZeros)
  neighbors = []  
  for source in nearby:
  
    # Get the distance from a source to the edge of the aperture
    # Not terribly robust, but sufficient for an estimate
    if source.epic != EPIC:
      dist = np.inf
      m, n = int(np.round(source.y - source.y0)), int(np.round(source.x - source.x0))
      for dm in range(-3,4):
        for dn in range(-3,4):
          if m + dm >= 0 and m + dm < ny and n + dn >= 0 and n + dn < nx:
            if contour[m + dm, n + dn]:
              d = np.sqrt(dm ** 2 + dn ** 2)
              if d < dist:
                dist = d
      if not np.isinf(dist):
        neighbors.append((dist, source.epic, source.kepmag))

  # Generate some info about possible contaminant sources
  neighbors = sorted(neighbors)
  crwdsev = 0
  crwdinfo = ["#   EPIC       D (PIX)   ΔKp (mag)   SEVERITY (0-5)",
              "# ---------    -------   ---------   --------------"]
  for d, e, k in neighbors:
    dk = k - kepmag
    if not np.isnan(dk):
      # Assess the severity of the contamination
      if d == 0 and dk < 0:
        note = "5 (PLD FAILURE)"
        if crwdsev < 5: crwdsev = 5
      elif d == 0 and dk < 3:
        note = "4 (LIKELY OVERFITTING)"
        if crwdsev < 4: crwdsev = 4
      elif d == 0 and dk < 5:
        note = "3 (MINOR OVERFITTING)"
        if crwdsev < 3: crwdsev = 3
      elif d == 0 and dk > 5:
        note = "1"
        if crwdsev < 1: crwdsev = 1
      elif d <= 1. and dk < 0:
        note = "4 (LIKELY OVERFITTING)"
        if crwdsev < 4: crwdsev = 4
      elif d <= 1. and dk < 3:
        note = "2"
        if crwdsev < 2: crwdsev = 2
      elif d <= 1. and dk < 5:
        note = "1"
        if crwdsev < 1: crwdsev = 1
      elif d <= 2. and dk < 0:
        note = "3 (MINOR OVERFITTING)"
        if crwdsev < 3: crwdsev = 3
      elif d <= 2. and dk < 5:
        note = "2"
        if crwdsev < 2: crwdsev = 2
      elif d <= 3. and dk < 0:
        note = "2"
        if crwdsev < 2: crwdsev = 2
      else:
        note = "0"
      # Make it into a string
      if dk < 0:
        dk = "%.3f" % dk
      else:
        dk = "+%.3f" % dk
    else:
      dk = "????"
      note = ""
    crwdinfo.append("# %d    %.2f      %s      %s" % (e, d, dk, note))
  crwdinfo.append("#")
  crwdinfo = "\n".join(crwdinfo)
  
  return crwdinfo, crwdsev