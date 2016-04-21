#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
gp.py
-----

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .detrend import PLDCoeffs, PLDModel, PLDBasis
from .utils import Mask, LatexExp, MedianFilter, Smooth, Chunks
from .kernels import KernelModels
import os
import warnings
warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered')
warnings.filterwarnings('ignore', r'Polyfit may be poorly conditioned')
warnings.filterwarnings('ignore', r'Covariance of the parameters could not be estimated')
from statsmodels.tsa.stattools import acf
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit, fmin_l_bfgs_b
from astroML.time_series import lomb_scargle, lomb_scargle_bootstrap
from itertools import combinations
import numpy as np
import george
import logging
log = logging.getLogger(__name__)

def IsAlias(periods, p, tol = 0.05):
  '''
  Is ``p`` a 1/2x or 2x alias of one of ``periods``?
  
  '''
  
  for period in periods:
    for fac in [0.5, 2.0]:
      if np.abs(fac * period - p) < tol * period:
        return True
  return False

def GetGP(EPIC, time, fpix, ferr, mask = [], niter = 1):
  '''
  Fit various kernels to the autocorrelation of the PLD-de-trended data
  and return a bunch of stuff.
  
  TODO: If there are too many gaps, the autocorrelation function won't be
        correct, since it assumes evenly sampled data. We should eventually fill
        in the gaps with synthetic data.
  
  '''
  
  # Setup the mask
  mask = Mask(mask)
  
  # Hard-coded stuff
  maxp = 3
  maxt = 20.
  minchisq = 30.
  default_periods = [50., 25., 10.]
  
  # A rather arbitrary standard deviation function for curve_fit below.
  # We basically want the tightest fit closest to a lag of zero.
  sigma = lambda n: np.sqrt(np.linspace(0.05 ** 2, 0.5 ** 2, n))
  
  # Mask the transits and outliers now
  time = mask(time)
  fpix = mask(fpix)
  ferr = mask(ferr)
  flux = np.sum(fpix, axis = 1)
  
  # The first time through, de-trend with a 2-day Matern 3/2 kernel to get approximate 
  # stellar signal. We're using 1st order PLD and splitting the light curve into 10 chunks
  nchunks = 10
  amp = np.median([np.std(yi) for yi in 
                   Chunks(flux, int(2. / np.median(time[1:] - time [:-1])))])
  gp = george.GP(amp ** 2 * george.kernels.Matern32Kernel(2. ** 2))
  
  # We're going to do this iteratively, refining the GP each time
  for iter in range(niter):
    
    # De-trend
    brkpts = [time[int(i)] for i in np.linspace(0, len(time), nchunks + 1)[1:-1]]
    X, _ = PLDBasis(fpix, time = time, pld_order = 1, max_components = 50, 
                    breakpoints = brkpts)
    C = PLDCoeffs(X, flux, time, ferr, gp)
    M = PLDModel(C, X)
    fpld = flux - M
  
    # Now remove long-period trends and smooth the data to get ``fpld_med``
    fpld_med = Smooth(fpld - np.poly1d(np.polyfit(time, fpld, 3))(time), 13)
  
    # Run a periodogram analysis to get the peak periods. We'll use these as initial
    # guesses for the kernel curve fitting routine below.
    period = np.linspace(0.1, maxt, 1000)
    omega = 2 * np.pi / period
    PS = lomb_scargle(time, fpld_med, ferr, omega, generalized=True)
    D = lomb_scargle_bootstrap(time, fpld_med, ferr, omega, generalized=True,
                               N_bootstraps=100, random_state=0)
    sig1 = np.percentile(D, 99)
  
    # Get the peak periods (where the power is high and the derivative changes sign)
    PS_ = np.array(PS) 
    PS_[PS_ < sig1] = sig1
    deriv = PS_[1:] - PS_[:-1]
    peaks = np.where((np.sign(deriv[1:]) == -1)  & (np.sign(deriv[:-1]) == 1))[0]
    strength = [np.max([PS_[peak - 1], PS_[peak], PS_[peak + 1]]) for peak in peaks]
    pers = [foo[0] for foo in sorted(zip(period[peaks], strength), key = lambda x: -x[1])]
    
    # Remove any easy aliases
    aliases = []
    for i in range(1, len(pers)):
      if IsAlias(pers[:i], pers[i]):
        aliases.append(i)
    if len(aliases):
      pers = list(np.delete(pers, aliases))
  
    # Add the default periods if we have less than 3
    # and limit to ``maxp`` periods
    if len(pers) < 3:
      pers += default_periods
      pers = pers[:3]
    pers = pers[:maxp]
  
    # Make a list of period combinations of length ``n``
    pcomb = lambda n: list(combinations(pers, n))
  
    # Get the autocorrelation function of the de-trended data up to ``maxt`` days
    dt = np.median(time[1:] - time[:-1])
    tfull = np.arange(time[0], time[-1], dt)
    pfull = interp1d(time, fpld)(tfull)
    acor = acf(pfull, nlags = maxt / dt)
    lags = np.arange(len(acor)) * dt
    sigma = sigma(len(lags))
  
    # Fit various kernels to the autocorrelation function. These are all in
    # ``KernelModels`` in the ``kernel`` module.
    # Find the chisq based on only the first half (``maxt`` / 2 ~ 10 days) of the data.
    # We don't care that much about longer timescales.
    kernel = None
    chisq = np.inf
    knum = None
    count = 0
    white_guess = np.median([np.std(c) for c in Chunks(fpld, 13)]) / np.std(fpld)
  
    # We will repeat this a few times until we get a chisq
    # that's acceptable.
    while count < 5:
      for kn, k in enumerate(KernelModels):
        # Get indices of each of the param types
        perpars = np.where(k.parnames == 'per')[0]
        amppars = np.where(k.parnames == 'amp')[0]
        taupars = np.where(k.parnames == 'tau')[0]
    
        # Initialize our guess array
        p0 = np.zeros_like(k.params)
        p0[amppars[0]] = white_guess
        p0[amppars[1:]] = 1. / np.sqrt(len(amppars) - 1.)
    
        # Find all the possible period combinations
        pcombs = pcomb(len(perpars))
        if len(perpars) == 0: pcombs = [0.]
  
        # Loop over them
        for p in pcombs:
          if len(perpars): 
            # Add a larger perturbation every time we have to repeat this
            p0[perpars] = np.array(p) + count * np.array(p) * 0.1 * \
                          np.random.randn(len(perpars))
      
          # Randomize a timescale; make the standard deviation larger
          # every time we have to repeat this.
          p0[taupars] = np.abs(15. + (2 * count + 1.) * np.random.randn(len(taupars)))
      
          # Fit the kernel function to the autocorrelation function
          def kfunc(t, *params):
            k[:] = np.abs(params)
            return k(t)
            
          try:  
            params, _ = curve_fit(kfunc, lags, acor, sigma = sigma, 
                                  p0 = p0, maxfev = 10000)
          except:
            # Some of the fits may not converge; just move on to the next one
            continue
        
          # Now compute the chi^2 of the model. If it's better than the previous
          # value, store this kernel
          c = np.sum(((acor[:len(lags) // 2] - kfunc(lags[:len(lags) // 2], *params)) 
                       / sigma[:len(lags) // 2]) ** 2)
          if c < chisq:
            k[:] = np.abs(params)
            kernel = k
            knum = kn
            chisq = c
    
      count += 1
      if chisq < minchisq:
        break
    
    # Were we able to get the params?
    if kernel is None:
      raise Exception('Unable to determine GP params!')
  
    # Now we try to figure out the best kernel amplitude and white
    # noise component by maximizing the likelihood of the
    # de-trended data.
    amppars = np.where(kernel.parnames == 'amp')[0]
    kpars = np.array(kernel.params)

    # We'll maximize the likelihood using george
    log.info('Optimizing the kernel amplitude...')
    def _NegLL(x):
      w, a = x
      kp = np.array(kpars)
      kp[amppars[0]] = w
      kp[amppars[1:]] *= a
      kernel[:] = kp
      gp = george.GP(kernel.george_kernel())
      gp.compute(time, ferr)
      nll = -gp.lnlikelihood(fpld)
      return nll
    # The median 6-hr standard deviation will be our white noise amplitude guess
    # The full standard deviation will be our red kernel guess
    init = [np.median([np.std(c) for c in Chunks(fpld, 13)]), np.std(fpld)]
    bounds = [[0.1 * init[0], 10 * init[0]], [1., 100 * init[1]]]
    x = fmin_l_bfgs_b(_NegLL, init, approx_grad = True, bounds = bounds) 
    white, amp = x[0]
    kpars[amppars[0]] = white
    kpars[amppars[1:]] *= amp
    kernel[:] = kpars
    
    # Recalculate the chisq here. Shouldn't have to do this, but I was getting
    # strange values otherwise...
    chisq = np.sum(((acor - kernel(lags) / kernel(0.)) ** 2 / 
                     sigma ** 2)[:len(lags) // 2])
    
    # Next time around, we'll use the new GP and 5 chunks
    gp = george.GP(kernel.george_kernel())
    nchunks = 5
    
  # We're going to output some additional info for plotting
  return dict(data = [knum, kpars, white, amp, chisq, fpld],
              powerspec = [period, PS, sig1, pers],
              acor = [acor, lags, sigma, count],
              kernfunc = [kernel(lags) / kernel(0.), str(kernel)])