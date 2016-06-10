#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`plot.py` - Plotting routines
-------------------------------------

These are routines for generating the output plots during an :py:mod:`everest`
run.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .detrend import PLDCoeffs, PLDModel, PLDBasis
from .utils import RMS, Mask, PadWithZeros, LatexExp
from .data import RemoveBackground, Campaign
from .transit import Transit
from .sources import Source
from .config import EVEREST_DAT, EVEREST_SRC
from .crowding import Contamination
from k2plr.config import KPLR_ROOT
import os
import matplotlib
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator, ScalarFormatter
from scipy.ndimage import zoom
from scipy.interpolate import interp1d
import numpy as np
import logging
log = logging.getLogger(__name__)

def Plot(data):
  '''
  The main plotting function. This calls all other plotting functions below.
  
  :param dict data: The `data` dictionary returned by a call to :py:func:`everest.compute.Compute`
  
  '''
  
  EPIC = data['EPIC']
  outdir = data['outdir'][()]
  outdir = os.path.join(EVEREST_DAT, outdir[outdir.find(os.path.join('everest', 'output')) + 8:])
  jpeg_quality = int(data['jpeg_quality'][()])
  if matplotlib.get_backend() == 'MacOSX':
    ext = 'png'
  else:
    try:
      ext = data['fig_ext'][()]
    except KeyError:
      ext = 'jpg'
  
  # Plot the contamination
  log.info('Plotting the contamination...')
  if not os.path.exists(os.path.join(outdir, 'aper.%s' % ext)):
    PlotContamination(EPIC, run_name = data['run_name'], apnum = data['apnum'], 
                      fig_ext = ext, jpeg_quality = jpeg_quality)
  
  # Plot the apertures
  log.info('Plotting the apertures...')
  if not os.path.exists(os.path.join(outdir, 'aper.%s' % ext)):
    fig, ax = PlotApertures(EPIC, data)
    fig.savefig(os.path.join(outdir, 'aper.%s' % ext), quality = jpeg_quality)
    pl.close()

  # Plot the outliers
  log.info('Plotting the outliers...')
  if not os.path.exists(os.path.join(outdir, 'outliers.%s' % ext)):
    fig, ax = PlotOutliers(EPIC, data)
    fig.savefig(os.path.join(outdir, 'outliers.%s' % ext), quality = jpeg_quality)
    pl.close()  

  # Plot the apertures
  log.info('Plotting the GP...')
  if not os.path.exists(os.path.join(outdir, 'acor.%s' % ext)):
    fig, ax = PlotGP(EPIC, data)
    fig.savefig(os.path.join(outdir, 'acor.%s' % ext), quality = jpeg_quality)
    pl.close()

  # Plot the scatter curves
  if len(data['masked_scatter']):
    log.info('Plotting the scatter curve...')
    if not os.path.exists(os.path.join(outdir, 'scatter.%s' % ext)):
      fig, ax = PlotScatter(EPIC, data)
      fig.savefig(os.path.join(outdir, 'scatter.%s' % ext), quality = jpeg_quality)
      pl.close()

  # Plot the detrended data
  log.info('Plotting the detrended data...')
  if not os.path.exists(os.path.join(outdir, 'detrended.%s' % ext)):
    fig, ax = PlotDetrended(EPIC, data)
    fig.savefig(os.path.join(outdir, 'detrended.%s' % ext), quality = jpeg_quality)
    pl.close()
  
  # ----
  # DEBUG: REPLOT OLD FIGURES
  elif os.path.getmtime(os.path.join(outdir, 'detrended.%s' % ext)) < 1464740277:
    fig, ax = PlotDetrended(EPIC, data)
    fig.savefig(os.path.join(outdir, 'detrended.%s' % ext), quality = jpeg_quality)
    pl.close()
  # /DEBUG
  # ----
  
  # Plot the folded data
  if data['mask_candidates'] or len(data['inject'][()]):
    log.info('Plotting the folded data...')
    try:
      PlotFolded(EPIC, data)
    except:
      log.error('An error occurred while plotting the folded data.')
  
  log.info('Done!')

def PlotContamination(EPIC, run_name = 'default', apnum = 15, 
                      fig_ext = 'jpg', jpeg_quality = 30):
  '''
  
  '''
  
  campaign = Campaign(EPIC)
  filename = os.path.join(KPLR_ROOT, 'data', 'everest', str(EPIC), str(EPIC) + '.npz')
  data = np.load(filename)
  fpix = data['fpix']
  perr = data['perr']
  _nearby = data['nearby']
  nearby = [Source(**s) for s in _nearby]
  apertures = data['apertures']
  apidx = np.where(apertures[apnum] & 1 & ~np.isnan(fpix[0])) 
  bkidx = np.where(apertures[apnum] ^ 1) 
  contamination, fig, ax = Contamination(EPIC, fpix, perr, apidx, bkidx, nearby, plot = True)
  fig.savefig(os.path.join(EVEREST_DAT, 'output', 'C%02d' % campaign, '%d' % EPIC, 
              run_name, 'contamination.%s' % fig_ext), quality = jpeg_quality)
  pl.close()
  
def PlotScatter(EPIC, data):
  '''
  Plots the cross-validation scatter plot for a given run.
  
  :param int EPIC: The 9-digit `EPIC` number of the target
  
  :param dict data: The `data` dictionary returned by a call to :py:func:`everest.compute.Compute`
  
  '''
  
  flux = data['flux']
  npc_arr = data['npc_arr']
  npc_pred = data['npc_pred']
  masked_scatter = data['masked_scatter']
  unmasked_scatter = data['unmasked_scatter']
  msf = data['msf']
  usf = data['usf']
  besti = data['besti']
  r1, r2, r3, r4, r5 = data['rms']
  
  # Plot the scatter as a function of the number of components
  fig, ax = pl.subplots(1, figsize = (12, 8))
  fig.subplots_adjust(hspace = 0.3, wspace = 0.2, left = 0.1, right = 0.95)
  
  # Plot the scatter curves
  ax.plot(npc_arr, masked_scatter, 'r.', alpha = 0.25)
  ax.plot(npc_arr, unmasked_scatter, 'b.', alpha = 0.25)
  ax.plot(npc_pred, msf, 'r-', label = 'Masked')
  ax.plot(npc_pred, usf, 'b-', label = 'Unmasked') 
  
  # Add labels   
  ax.set_ylabel('Precision (ppm)', fontsize = 18)
  ax.set_xlabel('Number of principal components', fontsize = 18)
  pl.suptitle('EPIC %d (%.1f ppm)' % (EPIC, r1), fontsize = 25)
  
  # Appearance
  ax.margins(0.01, 0.01)
  ax.legend(loc = 'best', fontsize = 12)
  for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(12)
  
  # Finally, mark the best run
  ax.axvline(npc_pred[besti], color = 'k', ls = '--')
  ax.axhline(msf[besti], color = 'k', ls = '--')
  ax.set_axis_bgcolor('#e6e6ff')

  return fig, ax

def PlotGP(EPIC, data):
  '''
  Plots the GP optimization diagnostic plot for a given run.
  
  :param int EPIC: The 9-digit `EPIC` number of the target
  
  :param dict data: The `data` dictionary returned by a call to :py:func:`everest.compute.Compute`
  
  '''
  
  time = data['time']
  fpld = data['fpldgp']
  ferr = data['ferr']
  mask = data['mask']
  acor, lags, sigma, count = data['acor']
  period, PS, sig1, pers = data['powerspec']
  chisq  = data['chisq']
  white = data['white']
  amp = data['amp']
  kernfunc, kernstr = data['kernfunc']
  gp_iter = data['gp_iter']
  
  time = np.delete(time, mask, axis = 0)
  dt = np.median(time[1:] - time[:-1])
  tfull = np.arange(time[0], time[-1], dt)
  pfull = interp1d(time, fpld)(tfull)
  
  # Plot
  fig, ax = pl.subplots(3, figsize = (16, 12))
  fig.subplots_adjust(left = 0.1, right = 0.95, hspace = 0.3)
  
  # The detrended data
  ax[0].plot(time, fpld, 'k.', alpha = 0.3)
  ax[0].plot(tfull, pfull, 'r-', alpha = 0.3)
  ax[0].margins(0.01, None)
  
  # The periodogram
  ax[1].plot(period, PS, '-', c='black', lw=1, zorder=1)
  ax[1].plot([period[0], period[-1]], [sig1, sig1], ':', c='black')
  if len(pers):
    ax[1].axvline(pers[0], color = 'b', alpha = 0.5, label = 'Peak periods')
    [ax[1].axvline(p, color = 'b', alpha = 0.5) for p in pers[1:]]
  ax[1].legend(loc = 'upper right')
  ax[1].margins(0, None)
  
  # The autocorrelation function and fit
  ax[2].plot(lags, acor, 'k-', lw = 2, alpha = 1.)
  ax[2].fill_between(lags, acor - sigma, acor + sigma, alpha = 0.1, color = 'k')
  ax[2].plot(lags, kernfunc, '-', alpha = 1., label = kernstr)
  ax[2].legend(loc = 'upper right', fontsize = 14)
  ax[2].annotate(r'$\chi^2 = %.2f\ (%d\times)$' % (chisq, gp_iter), xy = (0.02, 0.96), xycoords = 'axes fraction',
                 fontsize = 14, verticalalignment = 'top', horizontalalignment = 'left')
  ax[2].annotate(r'$W = %.1f\sigma_w$' % (white / np.median(ferr)) + '\n' + 
                 r'$A = %.1f\sigma_a$' % (amp /  np.std(fpld)), 
                 xy = (0.02, 0.04), xycoords = 'axes fraction',
                 fontsize = 14, verticalalignment = 'bottom', 
                 horizontalalignment = 'left')                  
  ax[2].margins(0, None)

  # Labels
  pl.suptitle('EPIC %d' % EPIC, fontsize = 25)
  ax[0].set_ylabel('PLD-detrended data', fontsize = 18)
  ax[0].set_xlabel('Time (days)', fontsize = 18)
  ax[1].set_ylabel('Power', fontsize = 18)
  ax[1].set_xlabel('Period (days)', fontsize = 18)
  ax[2].set_ylabel('Autocorrelation', fontsize = 18)
  ax[2].set_xlabel('Time Lag (days)', fontsize = 18)

  return fig, ax

def PlotOutliers(EPIC, data):
  '''
  Plots the results of the iterative outlier detection step.
  
  :param int EPIC: The 9-digit `EPIC` number of the target
  
  :param dict data: The `data` dictionary returned by a call to :py:func:`everest.compute.Compute`
  
  '''

  time = data['time']
  flux = data['flux']
  mask = data['mask']
  campaign = data['campaign']
  if RemoveBackground(EPIC):
    bkg = data['bkg']
  else:
    bkg = np.zeros_like(time)
  trn_mask = np.array(data['trn_mask'], dtype = int)
  remove = np.array(data['out_mask'], dtype = int)
  
  fig, ax = pl.subplots(2, figsize = (16, 8))
  fig.subplots_adjust(left = 0.1, right = 0.95, hspace = 0.02)
  
  ax[0].plot(time, bkg, 'r.', alpha = 0.3, label = 'Background')
  ax[0].set_ylabel('Bkg Flux', fontsize = 18)
  ax[0].set_xticklabels([])
  ax[0].margins(0.01, None)
  
  ax[1].plot(np.delete(time, trn_mask), np.delete(flux, trn_mask), 'k.', alpha = 0.3)
  ax[1].plot(time[trn_mask], flux[trn_mask], 'b.', label = 'Masked (user)', alpha = 0.5)
  ax[1].plot(time[remove], flux[remove], 'r.', label = 'Removed from lightcurve', )
  ax[1].legend(loc = 'upper left', fontsize = 12, numpoints = 1)
  ax[1].margins(0.01, 0.05)
  ylim = ax[1].get_ylim()
  ax[1].set_ylim(ylim[0], ylim[1] + 0.05 * (ylim[1] - ylim[0]))
  ax[1].set_xlabel('Time (days)', fontsize = 18)
  ax[1].set_ylabel('Corr Flux', fontsize = 18)
  pl.suptitle('EPIC %d' % EPIC, fontsize = 25)
  
  return fig, ax

def PlotApertures(EPIC, data):
  '''
  Plots the target postage stamp with the chosen aperture overlaid. The colors
  (dark blue: low, dark red: high) represent the time-integrated flux over the
  entire campaign. This gives us an idea of how suitable the chosen static aperture
  is for capturing all the target's flux.
  
  :param int EPIC: The 9-digit `EPIC` number of the target
  
  :param dict data: The `data` dictionary returned by a call to :py:func:`everest.compute.Compute`
  
  '''
  
  # Grab the data
  nearby = [Source(**s) for s in data['nearby']]
  apertures = data['apertures']
  apnum = data['apnum']
  fpix = data['fpix_full']
  kepmag = data['kepmag']
  _, ny, nx = fpix.shape
  mxf = np.log10(np.nansum(fpix, axis = 0))

  # Plot each aperture
  fig = pl.figure(figsize = (16,8))
  fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.9, bottom = 0.05)
  axes = [pl.subplot2grid((20,50), (5 * i, 5 * j), colspan=4, rowspan=4) 
          for i in range(4) for j in range(5)]
  axes += [pl.subplot2grid((20,50), (0, 25), colspan = 25, rowspan = 19)]

  for i, ax in enumerate(axes):

    # Plot the data
    ax.imshow(mxf, interpolation = 'nearest', alpha = 0.75)

    # Aperture contours (a bit of a hack!)
    if i < 20:
      apidx = np.where(apertures[i] & 1)
    else:
      apidx = np.where(apertures[apnum] & 1)
    contour = np.zeros((ny,nx))
    contour[apidx] = 1

    # Add padding around the contour mask so that the edges get drawn
    contour = np.lib.pad(contour, 1, PadWithZeros)

    # Zoom in to make the contours look vertical/horizontal
    highres = zoom(contour, 100, order=0, mode='nearest') 
    extent = np.array([0, nx, 0, ny]) + \
             np.array([0, -1, 0, -1]) + \
             np.array([-1, 1, -1, 1])

    # Apply the contours
    ax.contour(highres, levels=[0.5], extent=extent, 
               origin='lower', colors='k', linewidths=3)

    # Another hack, to remove the padding
    ax.set_xlim(-0.7, nx - 0.3)
    ax.set_ylim(-0.7, ny - 0.3)
  
    # Appearance tweaks
    if i < 20:
      ax.set_title('%02d' % i, fontsize = 12)
    else:
      # Overplot nearby sources
      neighbors = []  
      def size(k):
        s = 1000 * 2 ** (kepmag - k)
        return s
      for source in nearby:
        ax.scatter(source.x - source.x0, source.y - source.y0, 
                   s = size(source.kepmag), 
                   c = ['g' if source.epic == EPIC else 'r'],
                   alpha = 0.5,
                   edgecolor = 'k')
        ax.scatter(source.x - source.x0, source.y - source.y0, 
                   marker = r'$%.1f$' % source.kepmag, color = 'w',
                   s = 500)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
  
  pl.suptitle('EPIC %d' % EPIC, fontsize = 25, y = 0.98)
  
  return fig, axes

def PlotDetrended(EPIC, data):
  '''
  Plots the de-trended data.
  
  :param int EPIC: The 9-digit `EPIC` number of the target
  
  :param dict data: The `data` dictionary returned by a call to :py:func:`everest.compute.Compute`
  
  '''
  
  # Detrended data
  time = data['time']
  flux = data['flux']
  ferr = data['ferr']
  fpix = data['fpix']
  perr = data['perr']
  bkg = data['bkg'] * fpix.shape[1]
  mask = data['mask']
  trn_mask = np.array(data['trn_mask'], dtype = int)
  rem_mask = np.array(data['out_mask'], dtype = int)
  apnum = data['apnum']
  inject = data['inject'][()]
  new_candidates = data['new_candidates'][()]
  mask_candidates = data['mask_candidates']
  C = data['C']
  pld_order = data['pld_order']
  fpld = data['fpld']
  fwhite = data['fwhite']
  npars = '%d free parameters, %d data points' % (np.product(C.shape), fwhite.shape[0])
  r1, r2, r3, r4, r5 = data['rms']
  kepmag = data['kepmag']
  campaign = data['campaign']
  planets = data['planets']
  breakpoints = data['breakpoints']
                            
  # Plot
  if mask_candidates or len(inject):
    fig, ax = pl.subplots(3, figsize = (16, 12))
  else:
    fig, ax = pl.subplots(2, figsize = (16, 12))
  fig.subplots_adjust(left = 0.08, right = 0.98, hspace = 0.15, wspace = 0.125)
  
  # 1. Raw flux, background removed
  l01, = ax[0].plot(time, flux, 'r.', alpha = 0.3, label = 'Raw Flux')
  l02, = ax[0].plot(time[trn_mask], flux[trn_mask], 'b.', alpha = 0.5, 
               label = 'User mask')
  l03, = ax[0].plot(time[rem_mask], flux[rem_mask], 'k.', alpha = 1, 
               label = 'Outlier mask')
  l04, = ax[0].plot([np.nan], [np.nan], 'w.', 
                    label = 'Raw:   %.1f ppm / %.1f ppm' % (r1, r2)) 
  l05, = ax[0].plot([np.nan], [np.nan], 'w.', 
                    label = 'Photon limit:   %.1f ppm' % (r5)) 
               
  # Adjust y limits to make room for the legend at the top left
  ylim = list(ax[0].get_ylim())
  top = np.max(flux[:len(flux) // 3])
  if (ylim[1] - top) < (ylim[1] - ylim[0]) / 3:
    ylim[1] += (ylim[1] - ylim[0]) / 3
    ax[0].set_ylim(ylim)
      
  # 2. PLD detrended, removing instrumental noise only
  ax[1].plot(time, fpld, 'b.', alpha = 0.3)
  l11, = ax[1].plot([np.nan], [np.nan], 'w.', label = 'Everest: %.1f ppm / %.1f ppm' % (r3, r4))
  
  # Adjust y limits, but don't let a few outliers throw us off
  m = np.median(fpld)
  lo = fpld[fpld < m]
  lo = m - 2. * (m - lo[np.argsort(lo)][10])
  hi = fpld[fpld > m]
  hi = m + 2. * (hi[np.argsort(hi)][-10] - m)
  top = np.max(fpld[:len(fpld) // 3])
  if (hi - top) < (hi - lo) / 3:
    hi += (hi - lo) / 3
  ax[1].set_ylim(lo, hi)
  
  if mask_candidates or len(inject):
    # 3. The fully whitened PLD-detrended data, only if we masked candidate
    # transits and want to see what they look like in the whitened data      
    ax[2].plot(time, fwhite, 'b.', alpha = 0.3)
  
  if mask_candidates:
    # Adjust y limits based on the Everest-detrended data
    y = np.delete(fwhite, rem_mask)
    lo = y[y < 1]
    lo = 1. - 2. * (1. - lo[np.argsort(lo)][10])
    hi = y[y > 1]
    hi = 1. + 2. * (hi[np.argsort(hi)][-10] - 1.)
    rng = hi - lo
    lo -= 0.1 * rng
    hi += 0.4 * rng
    ax[2].set_ylim(lo, hi)

  # Adjust limits
  ax[0].margins(0.01, None)
  ax[1].set_xlim(*ax[0].get_xlim())
  if mask_candidates or len(inject):
    ax[2].set_xlim(*ax[0].get_xlim())
  
  # Legends and labels
  if len(trn_mask):
    leg1 = ax[0].legend(handles=[l01, l02, l03], loc = 'upper right', 
                        fontsize = 12, numpoints = 3)
  else:
    leg1 = ax[0].legend(handles=[l01, l03], loc = 'upper right', 
                        fontsize = 12, numpoints = 3)
  ax[0].add_artist(leg1)
  ax[0].legend(handles=[l04, l05], loc='upper left', fontsize = 12, numpoints = 1, handlelength = 0, handletextpad = 0.5)
  ax[1].legend(handles=[l11], loc = 'upper left', fontsize = 12, numpoints = 1, handlelength = 0, handletextpad = 0.5)
  ax[0].set_ylabel('Raw', fontsize = 18)
  ax[1].set_ylabel('Everest', fontsize = 18)
  ax[1].annotate('PLD order = %d\n%s' % (pld_order, npars), 
                   xy = (0.975, 0.95),
                   horizontalalignment = 'right', 
                   verticalalignment = 'top',
                   xycoords = 'axes fraction', fontsize = 12, color = 'k')
  if mask_candidates or len(inject):
    ax[2].set_ylabel('Whitened', fontsize = 18)
  ax[-1].set_xlabel('Time (days)', fontsize = 18)
  pl.suptitle(r'EPIC %d' % EPIC, fontsize = 25)
  for axis in ax:
    axis.ticklabel_format(useOffset=False)
  for axis in ax[:-1]:
    axis.xaxis.set_ticklabels([])
  for label in ax[-1].xaxis.get_ticklabels():
    label.set_fontsize(14)
    
  # Mark the known planet candidates
  if mask_candidates or len(inject):
    for axis in ax:
      tlim = axis.get_xlim()
      xticks = []
      xticklabels = []
      for planet, color in zip(planets, ['r', 'b', 'g', 'k', 'y'] * 10):
        per = planet['pl_orbper']
      
        # This will look horrible if the period is too short!
        if per < 0.75:
          continue
      
        t0 = planet['pl_tranmid'] - 2454833
        t0 += np.ceil((tlim[0] - t0) / per) * per
        for ttime in np.arange(t0, tlim[-1], per):
          xticks.append(ttime)
          xticklabels.append(int(planet['epic_candname'][-2:]))
          axis.axvline(ttime, color = color, ls = '--', alpha = 0.15)
    
      # Mark planet injections (if any)
      if len(inject):
        per = inject['per']
        if per < 0.75:
          continue
        t0 = inject['t0']
        t0 += np.ceil((tlim[0] - t0) / per) * per
        for ttime in np.arange(t0, tlim[-1], per):
          xticks.append(ttime)
          xticklabels.append('J')
          axis.axvline(ttime, color = 'k', ls = '--', alpha = 0.15) 
    
      # Mark new planet candidates (if any)
      if new_candidates:
        for new_candidate in new_candidates:
          for ttime in new_candidate['ttimes']:
            xticks.append(ttime)
            xticklabels.append('N')
            axis.axvline(ttime, color = 'b', ls = '--', alpha = 0.15) 
    
      if len(xticks):
        axtn = axis.twiny()
        axtn.set_xlim(tlim)
        axtn.set_xticks(xticks)
        axtn.set_xticklabels(xticklabels, fontsize = 8)
  
  # Indicate the campaign and Kp
  fig.text(0.03, 0.975,"C%02d" % campaign, ha="left", va="top", fontsize=24, color='b')
  fig.text(0.08, 0.975,"Kp%.1f" % kepmag, ha="left", va="top", fontsize=24, color='b')
  
  # Mark the breakpoints
  # Indicate where the lightcurve was split
  s_old = 0
  for axis in ax:
    fill = False
    for s in np.append(breakpoints, [time[-1]]):
      if fill:
        axis.axvspan(s_old, s, alpha=0.05, color='k')
      fill = not fill
      s_old = s
  return fig, ax  

def PlotFolded(EPIC, data):
  '''
  Plots the data folded on the period(s) of any planet/EB candidates. This routine
  typically won't get called, unless you're running a transit injection test or
  you set `mask_candidates` to `True` when calling :py:func:`everest.compute.Compute`.
  
  :param int EPIC: The 9-digit `EPIC` number of the target
  
  :param dict data: The `data` dictionary returned by a call to :py:func:`everest.compute.Compute`
  
  '''
  
  # Detrended data
  time = data['time']
  flux = data['flux']
  ferr = data['ferr']
  fpix = data['fpix']
  perr = data['perr']
  bkg = data['bkg'] * fpix.shape[1]
  mask = data['mask']
  trn_mask = np.array(data['trn_mask'], dtype = int)
  rem_mask = np.array(data['out_mask'], dtype = int)
  apnum = data['apnum']
  inject = data['inject'][()]
  new_candidates = data['new_candidates'][()]
  C = data['C']
  pld_order = data['pld_order']
  fpld = data['fpld']
  fwhite = data['fwhite']
  npars = '%d free parameters, %d data points' % (np.product(C.shape), fwhite.shape[0])
  r1, r2, r3, r4, r5 = data['rms']
  kepmag = data['kepmag']
  campaign = data['campaign']
  planets = data['planets']
  EB = data['EB'][()]
  outdir = data['outdir'][()]
  outdir = os.path.join(EVEREST_DAT, outdir[outdir.find(os.path.join('everest', 'output')) + 8:])
  jpeg_quality = int(data['jpeg_quality'][()])
  if matplotlib.get_backend() == 'MacOSX':
    ext = 'png'
  else:
    try:
      ext = data['fig_ext'][()]
    except KeyError:
      ext = 'jpg'
      
  # Is the star is an eclipsing binary?
  if EB:
    for n, eclipse, t0, dur in zip([1,2], ['Primary', 'Secondary'], ['p0', 's0'], ['pdur', 'sdur']):
      fig, ax = pl.subplots(1, figsize = (16, 8))
      fig.subplots_adjust(wspace = 0.1, hspace = 0.25, left = 0.075, right = 0.95)
      ax.set_ylabel('Flux', fontsize = 18)
      ax.set_title('EPIC %d: %s Eclipse' % (EPIC, eclipse), fontsize = 18)
      ax.set_xlabel('Time (days)', fontsize = 18)
      foldp = lambda t: (t - EB[t0] - EB['period'] / 2.) % EB['period'] - EB['period'] / 2.
      ax.plot(foldp(time), fwhite, 'b.', alpha = 0.4)
        
      # Fix y lims
      try:
        x = np.delete(fwhite, rem_mask)[np.abs(foldp(np.delete(time, rem_mask))) < 1]
        lo = x[x < 1]
        lo = 1. - 2. * (1. - lo[np.argsort(lo)][10])
        hi = x[x > 1]
        hi = 1. + 2. * (hi[np.argsort(hi)][-10] - 1.)
        ax.set_ylim(lo, hi)
      except:
        log.warn('Unable to set y lims in folded plot.')

      # Mark the eclipse duration 
      ax.axvline(-EB[dur] / 2., color = 'k', ls = '--', alpha = 0.2)
      ax.axvline(EB[dur] / 2., color = 'k', ls = '--', alpha = 0.2)

      # Tweak appearance
      ax.set_xlim(-4 * EB[dur], 4 * EB[dur])
      ax.ticklabel_format(useOffset=False)

      # Save this figure
      fig.savefig(os.path.join(outdir, 'folded_EB%02d.%s' % (n, ext)), quality = jpeg_quality)
      pl.close()

  # Is this a planet host?
  if len(planets) or len(inject) or len(new_candidates):
  
    # Plot
    npl = len(planets) + (len(inject) > 0) + len(new_candidates)
    pcount = 1
    
    for i, planet in enumerate(planets):
      
      # Set up the figure
      fig, ax = pl.subplots(1, figsize = (16, 8))
      fig.subplots_adjust(wspace = 0.1, hspace = 0.25, left = 0.075, right = 0.95)
      ax.set_ylabel('Flux', fontsize = 18)
      
      # Labels
      ax.set_title("%s (%s)" % (planet['epic_candname'], planet['k2c_disp']), fontsize = 18)
      ax.set_xlabel('Time (days)', fontsize = 18)
  
      # This is true for some of the false positives in the catalog
      if np.isnan(planet['pl_orbper']):
        continue
      fold = lambda t: (t - (planet['pl_tranmid'] - 2454833) - planet['pl_orbper'] / 2.) % \
                        planet['pl_orbper'] - planet['pl_orbper'] / 2.
      ax.plot(fold(time), fwhite, 'b.', alpha = 0.4)

      # Fix y lims
      try:
        x = np.delete(fwhite, rem_mask)[np.abs(fold(np.delete(time, rem_mask))) < 1]
        lo = x[x < 1]
        lo = 1. - 2. * (1. - lo[np.argsort(lo)][10])
        hi = x[x > 1]
        hi = 1. + 2. * (hi[np.argsort(hi)][-10] - 1.)
        ax.set_ylim(lo, hi)
      except:
        log.warn('Unable to set y lims in folded plot.')
      
      # Mark our assumed transit duration
      if not np.isnan(planet['pl_trandur']):
        tdur = planet['pl_trandur'] * 1.2
      else:
        # We assume 4 hours
        tdur = 4. / 24.
      ax.axvline(-tdur / 2., color = 'k', ls = '--', alpha = 0.2)
      ax.axvline(tdur / 2., color = 'k', ls = '--', alpha = 0.2)
  
      # Tweak appearance
      ax.set_xlim(-4 * tdur, 4 * tdur)
      ax.ticklabel_format(useOffset=False)
    
      # Save this figure
      fig.savefig(os.path.join(outdir, 'folded_%02d.%s' % (pcount, ext)), quality = jpeg_quality)
      pcount += 1
      pl.close()
    
    # Is this a transit injection?
    if len(inject):
    
      # Get the recovered depth and folded data for each model
      true_depth = inject['depth']
      adepth, adepth_err, at, ad = inject['everest']
    
      # Set up the figure
      fig, ax = pl.subplots(1, figsize = (16, 8))
      fig.subplots_adjust(wspace = 0.1, hspace = 0.25, left = 0.075, right = 0.95)
      ax.set_ylabel('Flux', fontsize = 18) 
      ax.set_title("Injected Transit (depth = %s)" % LatexExp(true_depth), fontsize = 18)
      ax.set_xlabel('Time (days)', fontsize = 18)
  
      # Plot
      ax.plot(at, ad, 'k.', alpha = 0.4)
      
      # Plot the transit model
      tprime = np.linspace(-1, 1, 1000)
      inject_0 = dict(inject); inject_0['t0'] = 0.
      ax.plot(tprime, Transit(tprime, **inject_0), 'r-')
      
      # Fix y lims
      try:
        lo = ad[ad < 1]
        lo = 1. - 2. * (1. - lo[np.argsort(lo)][10])
        hi = ad[ad > 1]
        hi = 1. + 2. * (hi[np.argsort(hi)][-10] - 1.)
        ax.set_ylim(lo, hi)
      except:
        log.warn('Unable to set y lims in folded plot.')
      
      # Mark our assumed transit duration
      tdur = inject['dur'] * 1.2
      ax.axvline(-tdur / 2., color = 'k', ls = '--', alpha = 0.2)
      ax.axvline(tdur / 2., color = 'k', ls = '--', alpha = 0.2)
      
      # Indicate the recovered depths
      ax.annotate(r"Recovered depth: %s $\pm$ %s (%.3f$\times$)" % (LatexExp(adepth), LatexExp(adepth_err), adepth / true_depth), 
                     xy = (0.98, 0.05), xycoords = 'axes fraction', ha = 'right', va = 'bottom')
  
      # Tweak appearance
      ax.set_xlim(max(min(at), -1), min(max(at), 1))
      ax.ticklabel_format(useOffset=False)
    
      # Save this figure
      fig.savefig(os.path.join(outdir, 'folded_%02d.%s' % (pcount, ext)), quality = jpeg_quality)
      pcount += 1
      pl.close()
    
    # Is it a new candidate?
    for new_candidate in new_candidates:
    
      # Set up the figure
      fig, ax = pl.subplots(1, figsize = (16, 8))
      fig.subplots_adjust(wspace = 0.1, hspace = 0.25, left = 0.075, right = 0.95)
      ax.set_ylabel('Flux', fontsize = 18)
      ax.set_title("New Candidate", fontsize = 18)
      ax.set_xlabel('Time (days)', fontsize = 18)
      ttimes = new_candidate['ttimes']
      fold = lambda t: np.array([ti - ttimes[np.argmin(np.abs(ttimes - ti))] 
                                 for ti in np.atleast_1d(t)])
      ax.plot(fold(time), fwhite, 'b.', alpha = 0.4)

      # Fix y lims
      try:
        x = np.delete(fwhite, rem_mask)[np.abs(fold(np.delete(time, rem_mask))) < 1]
        lo = x[x < 1]
        lo = 1. - 2. * (1. - lo[np.argsort(lo)][10])
        hi = x[x > 1]
        hi = 1. + 2. * (hi[np.argsort(hi)][-10] - 1.)
        axtop.set_ylim(lo, hi)
        axmid.set_ylim(lo, hi)
      except:
        log.warn('Unable to set y lims in folded plot.')
      
      # Mark the transit duration
      tdur = new_candidate['tdur'] * 1.2
      ax.axvline(-tdur / 2., color = 'k', ls = '--', alpha = 0.2)
      ax.axvline(tdur / 2., color = 'k', ls = '--', alpha = 0.2)
      
      # Tweak appearance
      ax.set_xlim(-4 * tdur, 4 * tdur)
      ax.ticklabel_format(useOffset=False)
    
      # Save this figure
      fig.savefig(os.path.join(outdir, 'folded_%02d.%s' % (pcount, ext)), quality = jpeg_quality)
      pcount += 1
      pl.close() 