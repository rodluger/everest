#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
utils.py
--------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from scipy.signal import medfilt
import numpy as np
import sys, traceback, pdb
import logging
log = logging.getLogger(__name__)

class Mask(object):
  '''
  A simple masking object to remove transits, outliers, etc.
  
  '''
  
  def __init__(self, mask = [], axis = 0):
    self.mask = mask
    self.axis = axis
    
  def __call__(self, x):
    if x is not None:
      return np.delete(x, self.mask, axis = self.axis)
    else:
      return None

class NoPILFilter(logging.Filter):
  '''
  The PIL image module has a nasty habit of sending all sorts of 
  unintelligible information to the logger. Let's filter that out
  here.
  
  '''
  
  def filter(self, record):
    return not record.name == 'PIL.PngImagePlugin'

def InitLog(file_name = None, log_level = logging.DEBUG, screen_level = logging.CRITICAL):
  '''
  A little routine to initialize the logging functionality.
  
  '''
  
  # Initialize the logging
  root = logging.getLogger()
  root.handlers = []
  root.setLevel(logging.DEBUG)

  # File handler
  if file_name is not None:
    fh = logging.FileHandler(file_name)
    fh.setLevel(log_level)
    fh_formatter = logging.Formatter('%(asctime)s %(levelname)s [%(name)s]: %(message)s', datefmt="%m/%d/%y %H:%M:%S")
    fh.setFormatter(fh_formatter)
    fh.addFilter(NoPILFilter())    
    root.addHandler(fh)

  # Screen handler
  sh = logging.StreamHandler(sys.stdout)
  sh.setLevel(screen_level)
  sh_formatter = logging.Formatter('%(levelname)s [%(name)s]: %(message)s')
  sh.setFormatter(sh_formatter)
  sh.addFilter(NoPILFilter()) 
  root.addHandler(sh)

def PadWithZeros(vector, pad_width, iaxis, kwargs):
    '''
    Pad an array with zeros.
    
    '''
    
    vector[:pad_width[0]] = 0
    vector[-pad_width[1]:] = 0
    return vector

def Breakpoint(campaign, time, mask = []):
  '''
  Return the timestamp of the breakpoint for a given campaign, if any.
  
  '''
  
  # We don't want to break during a transit!
  M = Mask(mask)
  mtime = M(time)

  # K2 Campaign 1: force lightcurve split at t ~ 2017.5,
  # which is a mid-campaign data gap
  if campaign == 1:
    return mtime[np.argmin(np.abs(mtime - 2017.5))]
  else:
    return None

def RMS(y, win = 13, remove_outliers = False):
  '''
  Return the CDPP (rms) in ppm based on the median running standard deviation for
  a window size of 13 cadences (~6 hours) as in VJ14.
  
  '''

  if remove_outliers:
    # Remove 5-sigma outliers from data 
    # smoothed on a 1 day timescale
    ys = y - Smooth(y, 50)
    M = np.nanmedian(ys)
    MAD = 1.4826 * np.nanmedian(np.abs(ys - M))
    out = []                                                                 
    for i, _ in enumerate(y):
      if (ys[i] > M + 5 * MAD) or (ys[i] < M - 5 * MAD):
        out.append(i)    
    out = np.array(out, dtype = int)
    y = np.delete(y, out)
    
  return 1.e6 * np.nanmedian([np.std(yi)/np.sqrt(win) for yi in Chunks(y, win, all = True)])

def Outliers(time, flux, fpix = None, ferr = None, mask = [], kernel_size = 5, 
             sigma = 5):
  '''
  Return ``remove`` and ``keep``, indices of outliers we should remove from and
  keep in the light curve, respectively.
  
  '''
  
  # Apply the mask
  mask = Mask(mask)
  t = mask(time)
  f = mask(flux)
  p = mask(fpix)
  e = mask(ferr)
      
  # Outliers indices
  otimes = []
  rtimes = []
  ktimes = []
  res = []
  
  # First we apply a median filter to the data
  fs = f - MedianFilter(f, kernel_size)

  # Then we perform a median absolute deviation (MAD) cut
  M = np.median(fs)
  MAD = 1.4826 * np.median(np.abs(fs - M))
  for ti, fi in zip(t, fs):
    if (fi > M + sigma * MAD) or (fi < M - sigma * MAD):
      # Outlier times
      otimes.append(ti)
      # Residuals
      res.append(np.abs(fi - M))
  # Sort the outliers by how bad they are
  _, otimes = np.array(sorted(zip(res, otimes))[::-1]).T
  
  # If the user provided fpix and ferr, we're going to
  # do PLD iteratively to identify outliers that contribute
  # to increasing the scatter.
  if len(otimes) and fpix is not None and ferr is not None:
  
    # Define an RMS function for the detrended data
    # We're applying first order PLD with a 30th order polynomial
    def GetRMS(ti, pi, fi, ei, poly_order = 30):
      frac = pi / fi.reshape(-1, 1)
      tprime = (ti - ti[0]) / (ti[-1] - ti[0])
      poly = np.vstack([tprime ** n for n in range(1, poly_order + 1)]).T
      x = np.hstack([frac, poly])
      Kinv = np.diag(1. / ei ** 2)
      A = np.dot(np.dot(x.T, Kinv), x)
      B = np.dot(np.dot(x.T, Kinv), fi)
      c = np.linalg.solve(A, B)
      m = np.dot(c, x.T)
      return RMS((fi - m) / np.median(fi))
    
    # Loop over all the outliers
    for otime in otimes:
    
      # Get the outlier index in the **masked** array
      outlier = np.argmax(t == otime)
      
      # Get the indices in the vicinity (300 cadences) of the outlier
      a = max(0, outlier - 150)
      b = min(outlier + 150, len(time))
      
      # Get the baseline RMS of the **masked** chunk
      RMS0 = GetRMS(t[a:b], p[a:b], f[a:b], e[a:b])
    
      # The **masked** arrays, with the outlier removed
      tr = np.delete(t, outlier)[a:b - 1]
      pr = np.delete(p, outlier, axis = 0)[a:b - 1]
      fr = np.delete(f, outlier)[a:b - 1]
      er = np.delete(e, outlier)[a:b - 1]
    
      # Remove the outlier if the RMS improved
      if GetRMS(tr, pr, fr, er) < RMS0:
        rtimes.append(otime)  
      else:
        ktimes.append(otime)
        
    rtimes = np.array(rtimes, dtype = float)
    ktimes = np.array(ktimes, dtype = float)
  
  else:
    rtimes = np.array(otimes)
    ktimes = np.array([], dtype = float)
  
  # Grab the indices from the times and return
  remove = np.array(sorted([np.argmax(time == t) for t in rtimes]))
  keep = np.array(sorted([np.argmax(time == t) for t in ktimes]))

  return remove, keep
      
def Chunks(l, n, all = False):
  '''
  Returns a generator of consecutive ``n``-sized chunks of list ``l``.
  If ``all`` is ``True``, returns **all** ``n``-sized chunks in ``l``
  by iterating over the starting point.
  
  '''
  
  if all:
    jarr = range(0, n - 1)
  else:
    jarr = [0]
  
  for j in jarr:
    for i in range(j, len(l), n):
      if i + 2 * n <= len(l):
        yield l[i:i+n]
      else:
        if not all:
          yield l[i:]
        break

def Smooth(x, window_len = 100, window = 'hanning'):
  '''
  Smooth data by convolving on a given timescale.
    
  '''
  
  if window_len == 0:
    return np.zeros_like(x)
  s = np.r_[2 * x[0] - x[window_len - 1::-1], x, 2 * x[-1] - x[-1:-window_len:-1]]
  if window == 'flat':
    w = np.ones(window_len, 'd')
  else:  
    w = eval('np.' + window + '(window_len)')
  y = np.convolve(w / w.sum(), s, mode = 'same')
  return y[window_len:-window_len + 1]

def MedianFilter(x, kernel_size = 5):
  '''
  A wrapper around ``scipy.signal.medfilt``.
  
  '''
  
  if kernel_size % 2 == 0:
    kernel_size += 1
  return medfilt(x, kernel_size = kernel_size)

def LatexExp(f, minexp = 3):
  '''
  Returns the TeX version of a number in scientific notation.
  
  '''
  
  if np.abs(np.log10(f)) >= minexp:
    float_str = "{0:.1e}".format(f)
    base, exponent = float_str.split("e")
    return r"${0} \times 10^{{{1}}}$".format(base, int(exponent))
  else:
    if f >= 1:
      return r"$%.1f$" % f
    elif f >= 0.1:
      return r"$%.2f$" % f
    else:
      return r"$%.3f$" % f
      
def LatexExpSq(f, minexp = 3):
  '''
  Same as ``LatexExp``, but specifically for printing a kernel
  amplitude or timescale, where the squaring is made explicit.
  
  '''
  
  if np.abs(np.log10(f)) >= minexp:
    float_str = "{0:.1e}".format(f)
    base, exponent = float_str.split("e")
    return r"$\left({0} \times 10^{{{1}}}\right)^2$".format(base, int(exponent))
  else:
    return LatexExp(f, minexp)[:-1] + '^2$'

def GetMasks(time, flux, fpix, ferr, outlier_sigma, planets = [], 
             EB = None, mask_times = []):
  '''
  
  '''
  
  # Get the transit masks
  if len(planets):
    mask = []
    for planet in planets:
    
      # Get the transit duration
      if not np.isnan(planet.pl_trandur):
        tdur = planet.pl_trandur * 1.2
      else:
        # Assume 4 hours for safety...
        tdur = 4. / 24.

      # Get the transit times
      per = planet.pl_orbper
      t0 = planet.pl_tranmid - 2454833
      t0 += np.ceil((time[0] - tdur - t0) / per) * per
      ttimes = np.arange(t0, time[-1] + tdur, per)

      for t in ttimes:
        mask.extend(np.where(np.abs(time - t) < tdur / 2.)[0])
    
    mask = sorted(mask)

  else:
    mask = []

  # Get eclipsing binary masks
  if EB:
    mask = sorted(set(mask + EB.mask(time)))
    
  # Enforce user-defined masks
  if len(mask_times):
    m = [np.argmax(np.abs(time - t) < 0.001) for t in mask_times]
    mask_pld = sorted(set(mask_pld + m))

  # Mask additional astrophysical outliers (including transits!) in the SAP flux 
  # for PLD to work properly. If we don't do this, PLD will actually attempt to
  # correct for these outliers at the expense of increasing the white noise in the
  # PLD fit.
  trn_mask = list(mask)
  rem_mask, keep_mask = Outliers(time, flux, fpix = fpix, ferr = ferr, mask = trn_mask, 
                                 sigma = outlier_sigma)
  mask = sorted(set(trn_mask + list(rem_mask)))
    
  return mask, trn_mask, rem_mask, keep_mask

def ExceptionHook(exctype, value, tb):
  '''
  A custom exception handler.
  
  '''
  
  for line in traceback.format_exception_only(exctype, value):
    log.error(line.replace('\n', ''))
  for line in traceback.format_tb(tb):
    log.error(line.replace('\n', ''))
  sys.__excepthook__(exctype, value, tb)

def ExceptionHookPDB(exctype, value, tb):
  '''
  A custom exception handler, with PDB post-mortem for debugging.
  
  '''
  
  for line in traceback.format_exception_only(exctype, value):
    log.error(line.replace('\n', ''))
  for line in traceback.format_tb(tb):
    log.error(line.replace('\n', ''))
  sys.__excepthook__(exctype, value, tb)
  pdb.pm()