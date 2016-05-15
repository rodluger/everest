#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`utils.py` - General utils
----------------------------------

General utility functions called from various parts of the code.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from scipy.signal import medfilt
import numpy as np
import sys, traceback, pdb
import george
import logging
log = logging.getLogger(__name__)

class Mask(object):
  '''
  A simple masking object to remove transits, outliers, etc.
  
  :param list mask: The indices to mask
  :param int axis: The axis to mask. Default `0`
  
  :returns: The array `x` with the indices `mask` masked
  
  '''
  
  def __init__(self, mask = [], axis = 0):
    self.mask = mask
    self.axis = axis
    
  def __call__(self, x):
    if x is not None:
      return np.delete(x, self.mask, axis = self.axis)
    else:
      return None

class FunctionWrapper(object):
  '''
  A simple function wrapper class. Stores `args` and `kwargs` and
  allows an arbitrary function to be called with a single parameter `x`
  
  '''
  
  def __init__(self, f, *args, **kwargs):
  
    self.f = f
    self.args = args
    self.kwargs = kwargs
  
  def __call__(self, x):
  
    return self.f(x, *self.args, **self.kwargs)

class NoPILFilter(logging.Filter):
  '''
  The PIL image module has a nasty habit of sending all sorts of 
  unintelligible information to the logger. We filter that out here.
  
  '''
  
  def filter(self, record):
    return not record.name == 'PIL.PngImagePlugin'

def InitLog(file_name = None, log_level = logging.DEBUG, screen_level = logging.CRITICAL):
  '''
  A little routine to initialize the logging functionality.
  
  :param str file_name: The name of the file to log to. Default `None` (set internally by :py:mod:`everest`)
  :param int log_level: The file logging level (`0-50`). Default `10` (debug)
  :param int log_level: The screen logging level (`0-50`). Default `50` (critical)
  
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
  Pads an array with zeros. This is used primarily for plotting aperture contours.
  
  '''
  
  vector[:pad_width[0]] = 0
  vector[-pad_width[1]:] = 0
  return vector

def Breakpoints(campaign, time, mask = []):
  '''
  Return the timestamp of the breakpoint for a given campaign, if any.
  
  :param int campaign: The `K2` campaign number
  :param ndarray time: The array of timestamps
  :param list mask: The light curve indices we're going to mask. You can provide this \
                    to ensure that the breakpoint does not occur in a masked region \
                    of the data. Useful if these regions contain transits!
  
  '''
  
  # We don't want to break during a transit!
  M = Mask(mask)
  mtime = M(time)

  # K2 Campaign 1: force lightcurve split at t ~ 2017.5,
  # which is a mid-campaign data gap
  if campaign == 1:
    return [mtime[np.argmin(np.abs(mtime - 2017.5))]]
  else:
    return []

def RMS(y, win = 13, remove_outliers = False):
  '''
  Return the CDPP (the 6-hr rms) in ppm based on the median running standard deviation for
  a window size of 13 cadences (~6 hours) as in VJ14.
  
  :param ndarray y: The array whose CDPP is to be computed
  :param int win: The window size in cadences. Default `13` (6.5 hours)
  :param bool remove_outliers: Clip outliers at 5 sigma before computing the CDPP? Default `False`
  
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

def MADOutliers(t, f, sigma = 5, kernel_size = 5):
  '''
  Performs a `median absolute deviation <https://en.wikipedia.org/wiki/Median_absolute_deviation>`_
  cut to identify outliers. A median filter is first applied to the data to remove
  trends.
  
  :param ndarray t: The independent variable array (time)
  :param ndarray t: The dependent variable array (flux)
  :param int sigma: The tolerance in standard deviations. Default `5`
  :param int kernel_size: The size of the kernel used in the median filter. Default `5`
  
  :returns: The indices of the outliers in `f`
  
  '''
  
  tout = []
  fs = f - MedianFilter(f, kernel_size)
  M = np.median(fs)
  MAD = 1.4826 * np.median(np.abs(fs - M))
  for ti, fi in zip(t, fs):
    if (fi > M + sigma * MAD) or (fi < M - sigma * MAD):
      tout.append(ti)
  outliers = np.array(sorted([np.argmax(t == ti) for ti in tout]))
  
  return outliers
      
def Chunks(l, n, all = False):
  '''
  Returns a generator of consecutive `n`-sized chunks of list `l`.
  If `all` is `True`, returns **all** `n`-sized chunks in `l`
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
  
  :param ndarray x: The data array
  :param int window_len: The size of the smoothing window. Default `100`
  :param str window: The window type. Default `hanning`
  
  
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
  A silly wrapper around :py:func:`scipy.signal.medfilt`.
  
  '''
  
  if kernel_size % 2 == 0:
    kernel_size += 1
  return medfilt(x, kernel_size = kernel_size)

def LatexExp(f, minexp = 3):
  '''
  Returns the TeX version of a number `f` in scientific notation.
  
  :param float f: The number
  :param int minexp: The minimum absolute value of the log of `f` above which \
                     scientific notation is used. Default `3`
  
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
  Same as :py:func:`LatexExp`, but specifically for printing a kernel
  amplitude or timescale, where the squaring is made explicit.
  
  :param float f: The number
  :param int minexp: The minimum absolute value of the log of `f` above which \
                     scientific notation is used. Default `3`
  
  '''
  
  if np.abs(np.log10(f)) >= minexp:
    float_str = "{0:.1e}".format(f)
    base, exponent = float_str.split("e")
    return r"$\left({0} \times 10^{{{1}}}\right)^2$".format(base, int(exponent))
  else:
    return LatexExp(f, minexp)[:-1] + '^2$'

def ExceptionHook(exctype, value, tb):
  '''
  A custom exception handler that logs errors to file.
  
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