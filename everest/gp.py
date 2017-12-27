#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`gp.py` - Gaussian Processes
------------------------------------

Routines for optimizing the GP hyperparameters for a given light curve.

'''

from __future__ import division, print_function, absolute_import, \
     unicode_literals
from .mathutils import Chunks
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import fmin_l_bfgs_b
from scipy.signal import savgol_filter
import numpy as np
np.random.seed(48151623)
from distutils.version import LooseVersion, StrictVersion
import george
from george.kernels import Matern32Kernel, ExpSine2Kernel
try:
    george_version = StrictVersion(george.__version__)
except ValueError:
    george_version = LooseVersion(george.__version__)
    comp_version = LooseVersion("0.3.0")
else:
    comp_version = StrictVersion("0.3.0")
if george_version < comp_version:
    from george.kernels import WhiteKernel
    OLDGEORGE = True
else:
    OLDGEORGE = False
import logging
log = logging.getLogger(__name__)


def GP(kernel, kernel_params, white=False):
    '''

    '''

    if kernel == 'Basic':
        w, a, t = kernel_params
        if white:
            if OLDGEORGE:
                return george.GP(WhiteKernel(w ** 2) +
                                 a ** 2 * Matern32Kernel(t ** 2))
            else:
                return george.GP(a ** 2 * Matern32Kernel(t ** 2),
                                 white_noise=np.log(w ** 2),
                                 fit_white_noise=True)
        else:
            return george.GP(a ** 2 * Matern32Kernel(t ** 2))
    elif kernel == 'QuasiPeriodic':
        w, a, g, p = kernel_params
        if white:
            if OLDGEORGE:
                return george.GP(WhiteKernel(w ** 2) +
                                 a ** 2 * ExpSine2Kernel(g, p))
            else:
                return george.GP(a ** 2 * ExpSine2Kernel(g, p),
                                 white_noise=np.log(w ** 2),
                                 fit_white_noise=True)
        else:
            return george.GP(a ** 2 * ExpSine2Kernel(g, p))
    else:
        raise ValueError('Invalid value for `kernel`.')


def GetCovariance(kernel, kernel_params, time, errors):
    '''
    Returns the covariance matrix for a given light curve
    segment.

    :param array_like kernel_params: A list of kernel parameters \
          (white noise amplitude, red noise amplitude, and red noise timescale)
    :param array_like time: The time array (*N*)
    :param array_like errors: The data error array (*N*)

    :returns: The covariance matrix :py:obj:`K` (*N*,*N*)

    '''

    # NOTE: We purposefully compute the covariance matrix
    # *without* the GP white noise term
    K = np.diag(errors ** 2)
    K += GP(kernel, kernel_params, white=False).get_matrix(time)
    return K


def GetKernelParams(time, flux, errors, kernel='Basic', mask=[],
                    giter=3, gmaxf=200, guess=None):
    '''
    Optimizes the GP by training it on the current de-trended light curve.
    Returns the white noise amplitude, red noise amplitude,
    and red noise timescale.

    :param array_like time: The time array
    :param array_like flux: The flux array
    :param array_like errors: The flux errors array
    :param array_like mask: The indices to be masked when training the GP. \
           Default `[]`
    :param int giter: The number of iterations. Default 3
    :param int gmaxf: The maximum number of function evaluations. Default 200
    :param tuple guess: The guess to initialize the minimization with. \
           Default :py:obj:`None`

    '''

    log.info("Optimizing the GP...")

    # Save a copy of time and errors for later
    time_copy = np.array(time)
    errors_copy = np.array(errors)

    # Apply the mask
    time = np.delete(time, mask)
    flux = np.delete(flux, mask)
    errors = np.delete(errors, mask)

    # Remove 5-sigma outliers to be safe
    f = flux - savgol_filter(flux, 49, 2) + np.nanmedian(flux)
    med = np.nanmedian(f)
    MAD = 1.4826 * np.nanmedian(np.abs(f - med))
    mask = np.where((f > med + 5 * MAD) | (f < med - 5 * MAD))[0]
    time = np.delete(time, mask)
    flux = np.delete(flux, mask)
    errors = np.delete(errors, mask)

    # Initial guesses and bounds
    white = np.nanmedian([np.nanstd(c) for c in Chunks(flux, 13)])
    amp = np.nanstd(flux)
    tau = 30.0
    if kernel == 'Basic':
        if guess is None:
            guess = [white, amp, tau]
        bounds = [[0.1 * white, 10. * white],
                  [1., 10000. * amp],
                  [0.5, 100.]]
    elif kernel == 'QuasiPeriodic':
        if guess is None:
            guess = [white, amp, tau, 1., 20.]
        bounds = [[0.1 * white, 10. * white],
                  [1., 10000. * amp],
                  [1e-5, 1e2],
                  [0.02, 100.]]
    else:
        raise ValueError('Invalid value for `kernel`.')

    # Loop
    llbest = -np.inf
    xbest = np.array(guess)
    for i in range(giter):

        # Randomize an initial guess
        iguess = [np.inf for g in guess]
        for j, b in enumerate(bounds):
            tries = 0
            while (iguess[j] < b[0]) or (iguess[j] > b[1]):
                iguess[j] = (1 + 0.5 * np.random.randn()) * guess[j]
                tries += 1
                if tries > 100:
                    iguess[j] = b[0] + np.random.random() * (b[1] - b[0])
                    break

        # Optimize
        x = fmin_l_bfgs_b(NegLnLike, iguess, approx_grad=False,
                          bounds=bounds, args=(time, flux, errors, kernel),
                          maxfun=gmaxf)
        log.info('Iteration #%d/%d:' % (i + 1, giter))
        log.info('   ' + x[2]['task'].decode('utf-8'))
        log.info('   ' + 'Function calls: %d' % x[2]['funcalls'])
        log.info('   ' + 'Log-likelihood: %.3e' % -x[1])
        if kernel == 'Basic':
            log.info('   ' + 'White noise   : %.3e (%.1f x error bars)' %
                     (x[0][0], x[0][0] / np.nanmedian(errors)))
            log.info('   ' + 'Red amplitude : %.3e (%.1f x stand dev)' %
                     (x[0][1], x[0][1] / np.nanstd(flux)))
            log.info('   ' + 'Red timescale : %.2f days' % x[0][2])
        elif kernel == 'QuasiPeriodic':
            log.info('   ' + 'White noise   : %.3e (%.1f x error bars)' %
                     (x[0][0], x[0][0] / np.nanmedian(errors)))
            log.info('   ' + 'Red amplitude : %.3e (%.1f x stand dev)' %
                     (x[0][1], x[0][1] / np.nanstd(flux)))
            log.info('   ' + 'Gamma         : %.3e' % x[0][2])
            log.info('   ' + 'Period        : %.2f days' % x[0][3])
        if -x[1] > llbest:
            llbest = -x[1]
            xbest = np.array(x[0])

    return xbest


def NegLnLike(x, time, flux, errors, kernel):
    '''
    Returns the negative log-likelihood function and its gradient.

    '''

    gp = GP(kernel, x, white=True)
    gp.compute(time, errors)
    if OLDGEORGE:
        nll = -gp.lnlikelihood(flux)
        # NOTE: There was a bug on this next line! Used to be
        #
        #    ngr = -gp.grad_lnlikelihood(flux) / gp.kernel.pars
        #
        # But I think we want
        #
        # dlogL/dx =     dlogL/dlogx^2       * dlogx^2/dx^2 * dx^2/dx
        #          = gp.grad_lnlikelihood()  *     1/x^2    *   2x
        #          = 2 * gp.grad_lnlikelihood() / x
        #          = 2 * gp.grad_lnlikelihood() / np.sqrt(x^2)
        #          = 2 * gp.grad_lnlikelihood() / np.sqrt(gp.kernel.pars)
        #
        # (with a negative sign out front for the negative gradient).
        # So we probably weren't optimizing the GP correctly! This affects
        # all campaigns through C13. It's not a *huge* deal, since the sign
        # of the gradient was correct and the model isn't that sensitive to
        # the value of the hyperparameters, but it may have contributed to
        # the poor performance on super variable stars. In most cases it means
        # the solver takes longer to converge and isn't as good at finding
        # the minimum.
        ngr = -2 * gp.grad_lnlikelihood(flux) / np.sqrt(gp.kernel.pars)
    else:
        nll = -gp.log_likelihood(flux)
        ngr = -2 * gp.grad_log_likelihood(flux) / \
            np.sqrt(np.exp(gp.get_parameter_vector()))

    return nll, ngr
