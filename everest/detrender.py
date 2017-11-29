#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`detrender.py` - De-trending models
-------------------------------------------

This module contains the generic models used to de-trend light curves for
the various supported missions. Most of the functionality is implemented in
:py:class:`Detrender`, and specific de-trending methods are implemented as
subclasses. The default :py:obj:`everest` model is :py:class:`nPLD`.

'''

from __future__ import division, print_function, absolute_import, \
    unicode_literals
from . import missions
from .basecamp import Basecamp
from .config import EVEREST_DAT
from .utils import InitLog, Formatter, AP_SATURATED_PIXEL, AP_COLLAPSED_PIXEL
from .mathutils import Chunks, Scatter, SavGol, Interpolate
from .fits import MakeFITS
from .gp import GetCovariance, GetKernelParams, GP
from .dvs import DVS, CBV
import os
import sys
import numpy as np
import george
from scipy.optimize import fmin_powell
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
from PyPDF2 import PdfFileReader, PdfFileWriter
import traceback
import logging
log = logging.getLogger(__name__)

__all__ = ['Detrender', 'rPLD', 'nPLD', 'iPLD', 'pPLD']


class Detrender(Basecamp):
    '''
    A generic *PLD* model with scalar matrix *L2* regularization. Includes
    functionality for loading pixel-level light curves, identifying outliers,
    generating the data covariance matrix, computing the regularized pixel
    model, and plotting the results.
    Specific models are implemented as subclasses.

    **General:**

    :param ID: The target star ID (*EPIC*, *KIC*, or *TIC* number, \
           for instance)
    :param str cadence: The cadence of the observations. Default :py:obj:`lc`
    :param bool clobber: Overwrite existing :py:obj:`everest` models? Default \
           :py:obj:`False`
    :param bool clobber_tpf: Download and overwrite the saved raw TPF data? \
           Default :py:obj:`False`
    :param bool debug: De-trend in debug mode? If :py:obj:`True`, prints all \
           output to screen and enters :py:obj:`pdb` post-mortem mode for \
           debugging when an error is raised. Default :py:obj:`False`
    :param str mission: The name of the mission. Default `k2`

    **Detrender:**

    :param str aperture_name: The name of the aperture to use. These are \
           defined in the datasets and are mission specific. Defaults to \
           the mission default
    :param int bpad: When light curve breakpoints are set, the light curve \
           chunks must be stitched together at the end. To prevent kinks \
           and/or discontinuities, the chunks are made to overlap by \
           :py:obj:`bpad` cadences on either end. The chunks are then \
           mended and the overlap is  discarded. Default 100
    :param breakpoints: Add light curve breakpoints when de-trending? If \
           :py:obj:`True`, splits the light curve into chunks and de-trends \
           each one separately, then stitches them back and the end. This is \
           useful for missions like *K2*, where the light curve noise \
           properties are very different at the beginning and end of each \
           campaign. The cadences at which breakpoints are inserted are \
           specified in the :py:func:`Breakpoints` function \
           of each mission. Alternatively, the user may specify a list of \
           cadences at which to break up the light curve. Default \
           :py:obj:`True`
    :param int cbv_num: The number of CBVs to regress on during \
           post-processing. Default 1
    :param int cbv_niter: The number of :py:obj:`SysRem` iterations to \
           perform when computing CBVs. Default 50
    :param int cbv_win: The filter window size (in cadences) for smoothing \
           the CBVs. Default 999
    :param int cbv_order: The filter order for smoothing CBVs. Default 3
    :param int cdivs: The number of light curve subdivisions when \
           cross-validating. During each iteration, one of these subdivisions \
           will be masked and used as the validation set. Default 3
    :param str cv_min: The quantity to be minimized during cross-validation. \
           Default `MAD` (median absolute deviation). Can also be set to \
           `TV` (total variation).
    :param int giter: The number of iterations when optimizing the GP. \
           During each iteration, the minimizer is initialized with a \
           perturbed guess; after :py:obj:`giter` iterations, the step with \
           the highest likelihood is kept. Default 3
    :param int gmaxf: The maximum number of function evaluations when \
           optimizing the GP. Default 200
    :param float gp_factor: When computing the initial kernel parameters, \
           the red noise amplitude is set to the standard deviation of the \
           data times this factor. Larger values generally help with \
           convergence, particularly for very variable stars. Default 100
    :param array_like kernel_params: The initial value of the \
           :py:obj:`Matern-3/2` kernel parameters \
           (white noise amplitude in flux units, red noise amplitude in \
           flux units, and timescale in days). Default :py:obj:`None` \
           (determined from the data)
    :param bool get_hires: Download a high resolution image of the target? \
           Default :py:obj:`True`
    :param bool get_nearby: Retrieve the location of nearby sources? \
           Default :py:obj:`True`
    :param array_like lambda_arr: The array of :math:`\Lambda` values to \
           iterate over during the cross-validation step. :math:`\Lambda` \
           is the regularization parameter, or the standard deviation of \
           the Gaussian prior on the weights for each order of PLD. \
           Default ``10 ** np.arange(0,18,0.5)``
    :param float leps: The fractional tolerance when optimizing \
           :math:`\Lambda`. The chosen value of :math:`\Lambda` will be \
           within this amount of the minimum of the CDPP curve. \
           Default 0.05
    :param int max_pixels: The maximum number of pixels. Very large apertures \
           are likely to cause memory errors, particularly for high order \
           PLD. If the chosen aperture exceeds this many \
           pixels, a different aperture is chosen from the dataset. If no \
           apertures with fewer than this many pixels are available, an error \
            is thrown. Default 75
    :param bool optimize_gp: Perform the GP optimization steps? \
           Default :py:obj:`True`
    :param float osigma: The outlier standard deviation threshold. Default 5
    :param int oiter: The maximum number of steps taken during iterative \
           sigma clipping. Default 10
    :param planets: Any transiting planets/EBs that should be explicitly \
           masked during cross-validation. It is not \
           usually necessary to specify these at the cross-validation stage, \
           since deep transits are masked as outliers and shallow transits \
           do not affect the lambda optimization. However, it *is* necessary \
           to mask deep transits in short cadence mode, since these can \
           heavily bias the cross-validation scheme to lower values of \
           lambda, leading to severe underfitting. \
           This parameter should be a tuple or a list of tuples in the \
           form (`t0`, `period`, `duration`) \
           for each of the planets to be masked (all values in days).
    :param int pld_order: The pixel level decorrelation order. Default `3`. \
           Higher orders may cause memory errors
    :param str saturated_aperture_name: If the target is found to be \
           saturated, de-trending is performed \
           on this aperture instead. Defaults to the mission default
    :param float saturation_tolerance: The tolerance when determining whether \
           or not to collapse a column in the aperture. The column collapsing \
           is implemented in the individual mission modules. Default -0.1, \
           i.e., if a target is 10% shy of the nominal saturation level, it
           is considered to be saturated.
    :param transit_model: An instance or list of instances of \
           :py:class:`everest.transit.TransitModel`. If specified, \
           :py:obj:`everest` will include these in the regression when \
           calculating the PLD coefficients. The final instrumental light \
           curve model will **not** include the transit fits -- they are used \
           solely to obtain unbiased PLD coefficients. The best fit transit \
           depths from the fit are stored \
           in the :py:obj:`transit_depth` attribute of the model. \
           Default :py:obj:`None`.
    '''

    def __init__(self, ID, **kwargs):
        '''

        '''

        # Initialize logging
        self.ID = ID
        if kwargs.get('season', None) is not None:
            self._season = kwargs.get('season')
            if hasattr(self._season, '__len__'):
                raise AttributeError(
                    "Please choose a campaign/season for this target: %s."
                    % self._season)
        self._data = kwargs.get('data', None)
        self.cadence = kwargs.get('cadence', 'lc').lower()
        if self.cadence not in ['lc', 'sc']:
            raise ValueError("Invalid cadence selected.")
        self.mission = kwargs.get('mission', 'k2')
        self.clobber = kwargs.get('clobber', False)
        self.debug = kwargs.get('debug', False)
        self.is_parent = kwargs.get('is_parent', False)
        if not self.is_parent:
            screen_level = kwargs.get('screen_level', logging.CRITICAL)
            log_level = kwargs.get('log_level', logging.DEBUG)
            InitLog(self.logfile, log_level, screen_level, self.debug)
            log.info("Initializing %s model for %d." % (self.name, self.ID))

        # If this is a short cadence light curve, get the
        # GP params from the long cadence model. It would
        # take way too long and too much memory to optimize
        # the GP based on the short cadence light curve
        if self.cadence == 'sc':
            kernel_params = kwargs.get('kernel_params', None)
            if kernel_params is None:
                log.info("Loading long cadence model...")
                kwcpy = dict(kwargs)
                kwcpy.pop('cadence', None)
                kwcpy.pop('clobber', None)
                lc = self.__class__(ID, is_parent=True, **kwcpy)
                kernel_params = np.array(lc.kernel_params)
                del lc
            kwargs.update(
                {'kernel_params': kernel_params, 'optimize_gp': False})

        # Read general model kwargs
        self.lambda_arr = kwargs.get('lambda_arr', 10 ** np.arange(0, 18, 0.5))
        if self.lambda_arr[0] != 0:
            self.lambda_arr = np.append(0, self.lambda_arr)
        self.leps = kwargs.get('leps', 0.05)
        self.osigma = kwargs.get('osigma', 5)
        self.oiter = kwargs.get('oiter', 10)
        self.cdivs = kwargs.get('cdivs', 3)
        self.giter = kwargs.get('giter', 3)
        self.gmaxf = kwargs.get('gmaxf', 200)
        self.optimize_gp = kwargs.get('optimize_gp', True)
        self.kernel_params = kwargs.get('kernel_params', None)
        self.kernel = kwargs.get('kernel', 'Basic')
        assert self.kernel in ['Basic', 'QuasiPeriodic'], \
            "Kwarg `kernel` must be one of `Basic` or `QuasiPeriodic`."
        self.clobber_tpf = kwargs.get('clobber_tpf', False)
        self.bpad = kwargs.get('bpad', 100)
        self.aperture_name = kwargs.get('aperture', None)
        self.saturated_aperture_name = kwargs.get('saturated_aperture', None)
        self.max_pixels = kwargs.get('max_pixels', 75)
        self.saturation_tolerance = kwargs.get('saturation_tolerance', -0.1)
        self.gp_factor = kwargs.get('gp_factor', 100.)
        self.get_hires = kwargs.get('get_hires', True)
        self.get_nearby = kwargs.get('get_nearby', True)
        self.planets = kwargs.get('planets', [])
        if type(self.planets) is tuple and len(self.planets) == 3 and \
                not hasattr(self.planets[0], '__len__'):
                    self.planets = [self.planets]
        for planet in self.planets:
            assert len(planet) == 3, \
                "Planets must be provided as (`t0`, `per`, `dur`) tuples."
        # Handle breakpointing. The breakpoint is the *last* index of each
        # light curve chunk.
        bkpts = kwargs.get('breakpoints', True)
        if bkpts is True:
            self.breakpoints = np.append(self._mission.Breakpoints(
                self.ID, season=self.season, cadence=self.cadence), [999999])
        elif hasattr(bkpts, '__len__'):
            self.breakpoints = np.append(bkpts, [999999])
        else:
            self.breakpoints = np.array([999999])
        nseg = len(self.breakpoints)
        self.cv_min = kwargs.get('cv_min', 'mad').lower()
        assert self.cv_min in ['mad', 'tv'], "Invalid value for `cv_min`."
        self.cbv_num = kwargs.get('cbv_num', 1)
        self.cbv_niter = kwargs.get('cbv_niter', 50)
        self.cbv_win = kwargs.get('cbv_win', 999)
        self.cbv_order = kwargs.get('cbv_order', 3)

        # Get the pld order
        pld_order = kwargs.get('pld_order', 3)
        assert (pld_order > 0), "Invalid value for the de-trending order."
        self.pld_order = pld_order

        # Get the transit model
        self._transit_model = kwargs.get('transit_model', None)

        # Initialize model params
        self.lam_idx = -1
        self.lam = [[1e5] + [None for i in range(self.pld_order - 1)]
                    for b in range(nseg)]
        self.reclam = None
        self.recmask = []
        self.X1N = None
        self.XCBV = None
        self.cdpp_arr = np.array([np.nan for b in range(nseg)])
        self.cdppr_arr = np.array([np.nan for b in range(nseg)])
        self.cdppv_arr = np.array([np.nan for b in range(nseg)])
        self.cdpp = np.nan
        self.cdppr = np.nan
        self.cdppv = np.nan
        self.cdppg = np.nan
        self.neighbors = []
        self.loaded = False
        self._weights = None

        # Initialize plotting
        self.dvs = DVS(len(self.breakpoints), pld_order=self.pld_order)

        # Check for saved model
        if self.load_model():
            return

        # Setup (subclass-specific)
        self.setup(**kwargs)

        # Run
        self.run()

    @property
    def name(self):
        '''
        Returns the name of the current :py:class:`Detrender` subclass.

        '''

        if self.cadence == 'lc':
            return self.__class__.__name__
        else:
            return '%s.sc' % self.__class__.__name__

    @name.setter
    def name(self, value):
        '''

        '''

        raise NotImplementedError("Can't set this property.")

    def setup(self, **kwargs):
        '''
        A subclass-specific routine.

        '''

        pass

    def cv_precompute(self, mask, b):
        '''
        Pre-compute the matrices :py:obj:`A` and :py:obj:`B`
        (cross-validation step only)
        for chunk :py:obj:`b`.

        '''

        # Get current chunk and mask outliers
        m1 = self.get_masked_chunk(b)
        flux = self.fraw[m1]
        K = GetCovariance(self.kernel, self.kernel_params,
                          self.time[m1], self.fraw_err[m1])
        med = np.nanmedian(flux)

        # Now mask the validation set
        M = lambda x, axis = 0: np.delete(x, mask, axis=axis)
        m2 = M(m1)
        mK = M(M(K, axis=0), axis=1)
        f = M(flux) - med

        # Pre-compute the matrices
        A = [None for i in range(self.pld_order)]
        B = [None for i in range(self.pld_order)]
        for n in range(self.pld_order):
            # Only compute up to the current PLD order
            if self.lam_idx >= n:
                X2 = self.X(n, m2)
                X1 = self.X(n, m1)
                A[n] = np.dot(X2, X2.T)
                B[n] = np.dot(X1, X2.T)
                del X1, X2

        if self.transit_model is None:
            C = 0
        else:
            C = np.zeros((len(m2), len(m2)))
            mean_transit_model = med * \
                np.sum([tm.depth * tm(self.time[m2])
                        for tm in self.transit_model], axis=0)
            f -= mean_transit_model
            for tm in self.transit_model:
                X2 = tm(self.time[m2]).reshape(-1, 1)
                C += tm.var_depth * np.dot(X2, X2.T)
                del X2

        return A, B, C, mK, f, m1, m2

    def cv_compute(self, b, A, B, C, mK, f, m1, m2):
        '''
        Compute the model (cross-validation step only) for chunk :py:obj:`b`.

        '''

        A = np.sum([l * a for l, a in zip(self.lam[b], A)
                    if l is not None], axis=0)
        B = np.sum([l * b for l, b in zip(self.lam[b], B)
                    if l is not None], axis=0)
        W = np.linalg.solve(mK + A + C, f)
        if self.transit_model is None:
            model = np.dot(B, W)
        else:
            w_pld = np.concatenate([l * np.dot(self.X(n, m2).T, W)
                                    for n, l in enumerate(self.lam[b])
                                    if l is not None])
            model = np.dot(np.hstack(
                [self.X(n, m1) for n, l in enumerate(self.lam[b])
                 if l is not None]), w_pld)
        model -= np.nanmedian(model)

        return model

    def get_outliers(self):
        '''
        Performs iterative sigma clipping to get outliers.

        '''

        log.info("Clipping outliers...")
        log.info('Iter %d/%d: %d outliers' %
                 (0, self.oiter, len(self.outmask)))

        def M(x): return np.delete(x, np.concatenate(
            [self.nanmask, self.badmask, self.transitmask]), axis=0)
        t = M(self.time)
        outmask = [np.array([-1]), np.array(self.outmask)]

        # Loop as long as the last two outlier arrays aren't equal
        while not np.array_equal(outmask[-2], outmask[-1]):

            # Check if we've done this too many times
            if len(outmask) - 1 > self.oiter:
                log.error('Maximum number of iterations in ' +
                          '``get_outliers()`` exceeded. Skipping...')
                break

            # Check if we're going in circles
            if np.any([np.array_equal(outmask[-1], i) for i in outmask[:-1]]):
                log.error('Function ``get_outliers()`` ' +
                          'is going in circles. Skipping...')
                break

            # Compute the model to get the flux
            self.compute()

            # Get the outliers
            f = SavGol(M(self.flux))
            med = np.nanmedian(f)
            MAD = 1.4826 * np.nanmedian(np.abs(f - med))
            inds = np.where((f > med + self.osigma * MAD) |
                            (f < med - self.osigma * MAD))[0]

            # Project onto unmasked time array
            inds = np.array([np.argmax(self.time == t[i]) for i in inds])
            self.outmask = np.array(inds, dtype=int)

            # Add them to the running list
            outmask.append(np.array(inds))

            # Log
            log.info('Iter %d/%d: %d outliers' %
                     (len(outmask) - 2, self.oiter, len(self.outmask)))

    def optimize_lambda(self, validation):
        '''
        Returns the index of :py:attr:`self.lambda_arr` that minimizes the
        validation scatter in the segment with minimum at the lowest value
        of :py:obj:`lambda`, with
        fractional tolerance :py:attr:`self.leps`.

        :param numpy.ndarray validation: The scatter in the validation set \
               as a function of :py:obj:`lambda`

        '''

        maxm = 0
        minr = len(validation)
        for n in range(validation.shape[1]):
            # The index that minimizes the scatter for this segment
            m = np.nanargmin(validation[:, n])
            if m > maxm:
                # The largest of the `m`s.
                maxm = m
            # The largest index with validation scatter within
            # `self.leps` of the minimum for this segment
            r = np.where((validation[:, n] - validation[m, n]) /
                         validation[m, n] <= self.leps)[0][-1]
            if r < minr:
                # The smallest of the `r`s
                minr = r
        return min(maxm, minr)

    def fobj(self, y, y0, t, gp, mask):
        '''

        '''

        if self.cv_min == 'mad':
            # Note that we're computing the MAD, not the
            # standard deviation, as this handles extremely variable
            # stars much better!
            gpm, _ = gp.predict(y - y0, t[mask])
            fdet = (y[mask] - gpm) / y0
            scatter = 1.e6 * \
                (1.4826 * np.nanmedian(np.abs(fdet - np.nanmedian(fdet))) /
                 np.sqrt(len(mask)))
            return scatter
        elif self.cv_min == 'tv':
            # We're going to minimize the total variation instead
            return 1.e6 * np.sum(np.abs(np.diff(y[mask]))) / len(mask) / y0

    def cross_validate(self, ax, info=''):
        '''
        Cross-validate to find the optimal value of :py:obj:`lambda`.

        :param ax: The current :py:obj:`matplotlib.pyplot` axis instance to \
               plot the cross-validation results.
        :param str info: The label to show in the bottom right-hand corner \
               of the plot. Default `''`

        '''

        # Loop over all chunks
        ax = np.atleast_1d(ax)
        for b, brkpt in enumerate(self.breakpoints):

            log.info("Cross-validating chunk %d/%d..." %
                     (b + 1, len(self.breakpoints)))
            med_training = np.zeros_like(self.lambda_arr)
            med_validation = np.zeros_like(self.lambda_arr)

            # Mask for current chunk
            m = self.get_masked_chunk(b)

            # Check that we have enough data
            if len(m) < 3 * self.cdivs:
                self.cdppv_arr[b] = np.nan
                self.lam[b][self.lam_idx] = 0.
                log.info(
                    "Insufficient data to run cross-validation on this chunk.")
                continue

            # Mask transits and outliers
            time = self.time[m]
            flux = self.fraw[m]
            ferr = self.fraw_err[m]
            med = np.nanmedian(flux)

            # The precision in the validation set
            validation = [[] for k, _ in enumerate(self.lambda_arr)]

            # The precision in the training set
            training = [[] for k, _ in enumerate(self.lambda_arr)]

            # Setup the GP
            gp = GP(self.kernel, self.kernel_params, white=False)
            gp.compute(time, ferr)

            # The masks
            masks = list(Chunks(np.arange(0, len(time)),
                                len(time) // self.cdivs))

            # Loop over the different masks
            for i, mask in enumerate(masks):

                log.info("Section %d/%d..." % (i + 1, len(masks)))

                # Pre-compute (training set)
                pre_t = self.cv_precompute([], b)

                # Pre-compute (validation set)
                pre_v = self.cv_precompute(mask, b)

                # Iterate over lambda
                for k, lam in enumerate(self.lambda_arr):

                    # Update the lambda matrix
                    self.lam[b][self.lam_idx] = lam

                    # Training set
                    model = self.cv_compute(b, *pre_t)
                    training[k].append(
                        self.fobj(flux - model, med, time, gp, mask))

                    # Validation set
                    model = self.cv_compute(b, *pre_v)
                    validation[k].append(
                        self.fobj(flux - model, med, time, gp, mask))

            # Finalize
            training = np.array(training)
            validation = np.array(validation)
            for k, _ in enumerate(self.lambda_arr):

                # Take the mean
                med_validation[k] = np.nanmean(validation[k])
                med_training[k] = np.nanmean(training[k])

            # Compute best model
            i = self.optimize_lambda(validation)
            v_best = med_validation[i]
            t_best = med_training[i]
            self.cdppv_arr[b] = v_best / t_best
            self.lam[b][self.lam_idx] = self.lambda_arr[i]
            log.info("Found optimum solution at log(lambda) = %.1f." %
                     np.log10(self.lam[b][self.lam_idx]))

            # Plotting: There's not enough space in the DVS to show the
            # cross-val results for more than three light curve segments.
            if len(self.breakpoints) <= 3:

                # Plotting hack: first x tick will be -infty
                lambda_arr = np.array(self.lambda_arr)
                lambda_arr[0] = 10 ** (np.log10(lambda_arr[1]) - 3)

                # Plot cross-val
                for n in range(len(masks)):
                    ax[b].plot(np.log10(lambda_arr),
                               validation[:, n], 'r-', alpha=0.3)

                ax[b].plot(np.log10(lambda_arr),
                           med_training, 'b-', lw=1., alpha=1)
                ax[b].plot(np.log10(lambda_arr),
                           med_validation, 'r-', lw=1., alpha=1)
                ax[b].axvline(np.log10(self.lam[b][self.lam_idx]),
                              color='k', ls='--', lw=0.75, alpha=0.75)
                ax[b].axhline(v_best, color='k', ls='--', lw=0.75, alpha=0.75)
                ax[b].set_ylabel(r'Scatter (ppm)', fontsize=5)
                hi = np.max(validation[0])
                lo = np.min(training)
                rng = (hi - lo)
                ax[b].set_ylim(lo - 0.15 * rng, hi + 0.15 * rng)
                if rng > 2:
                    ax[b].get_yaxis().set_major_formatter(Formatter.CDPP)
                    ax[b].get_yaxis().set_major_locator(
                        MaxNLocator(4, integer=True))
                elif rng > 0.2:
                    ax[b].get_yaxis().set_major_formatter(Formatter.CDPP1F)
                    ax[b].get_yaxis().set_major_locator(MaxNLocator(4))
                else:
                    ax[b].get_yaxis().set_major_formatter(Formatter.CDPP2F)
                    ax[b].get_yaxis().set_major_locator(MaxNLocator(4))

                # Fix the x ticks
                xticks = [np.log10(lambda_arr[0])] + list(np.linspace(
                    np.log10(lambda_arr[1]), np.log10(lambda_arr[-1]), 6))
                ax[b].set_xticks(xticks)
                ax[b].set_xticklabels(['' for x in xticks])
                pad = 0.01 * \
                    (np.log10(lambda_arr[-1]) - np.log10(lambda_arr[0]))
                ax[b].set_xlim(np.log10(lambda_arr[0]) - pad,
                               np.log10(lambda_arr[-1]) + pad)
                ax[b].annotate('%s.%d' % (info, b), xy=(0.02, 0.025),
                               xycoords='axes fraction',
                               ha='left', va='bottom', fontsize=7, alpha=0.25,
                               fontweight='bold')

        # Finally, compute the model
        self.compute()

        # Tidy up
        if len(ax) == 2:
            ax[0].xaxis.set_ticks_position('top')
        for axis in ax[1:]:
            axis.spines['top'].set_visible(False)
            axis.xaxis.set_ticks_position('bottom')

        if len(self.breakpoints) <= 3:

            # A hack to mark the first xtick as -infty
            labels = ['%.1f' % x for x in xticks]
            labels[0] = r'$-\infty$'
            ax[-1].set_xticklabels(labels)
            ax[-1].set_xlabel(r'Log $\Lambda$', fontsize=5)

        else:

            # We're just going to plot lambda as a function of chunk number
            bs = np.arange(len(self.breakpoints))
            ax[0].plot(bs + 1, [np.log10(self.lam[b][self.lam_idx])
                                for b in bs], 'r.')
            ax[0].plot(bs + 1, [np.log10(self.lam[b][self.lam_idx])
                                for b in bs], 'r-', alpha=0.25)
            ax[0].set_ylabel(r'$\log\Lambda$', fontsize=5)
            ax[0].margins(0.1, 0.1)
            ax[0].set_xticks(np.arange(1, len(self.breakpoints) + 1))
            ax[0].set_xticklabels([])

            # Now plot the CDPP and approximate validation CDPP
            cdpp_arr = self.get_cdpp_arr()
            cdppv_arr = self.cdppv_arr * cdpp_arr
            ax[1].plot(bs + 1, cdpp_arr, 'b.')
            ax[1].plot(bs + 1, cdpp_arr, 'b-', alpha=0.25)
            ax[1].plot(bs + 1, cdppv_arr, 'r.')
            ax[1].plot(bs + 1, cdppv_arr, 'r-', alpha=0.25)
            ax[1].margins(0.1, 0.1)
            ax[1].set_ylabel(r'Scatter (ppm)', fontsize=5)
            ax[1].set_xlabel(r'Chunk', fontsize=5)
            if len(self.breakpoints) < 15:
                ax[1].set_xticks(np.arange(1, len(self.breakpoints) + 1))
            else:
                ax[1].set_xticks(np.arange(1, len(self.breakpoints) + 1, 2))

    def finalize(self):
        '''
        This method is called at the end of the de-trending, prior to
        plotting the final results.
        Subclass it to add custom functionality to individual models.

        '''

        pass

    def get_ylim(self):
        '''
        Computes the ideal y-axis limits for the light curve plot. Attempts to
        set the limits equal to those of the raw light curve, but if more than
        1% of the flux lies either above or below these limits, auto-expands
        to include those points. At the end, adds 5% padding to both the
        top and the bottom.

        '''

        bn = np.array(
            list(set(np.concatenate([self.badmask, self.nanmask]))), dtype=int)
        fraw = np.delete(self.fraw, bn)
        lo, hi = fraw[np.argsort(fraw)][[3, -3]]
        flux = np.delete(self.flux, bn)
        fsort = flux[np.argsort(flux)]
        if fsort[int(0.01 * len(fsort))] < lo:
            lo = fsort[int(0.01 * len(fsort))]
        if fsort[int(0.99 * len(fsort))] > hi:
            hi = fsort[int(0.99 * len(fsort))]
        pad = (hi - lo) * 0.05
        ylim = (lo - pad, hi + pad)
        return ylim

    def plot_lc(self, ax, info_left='', info_right='', color='b'):
        '''
        Plots the current light curve. This is called at several stages to
        plot the de-trending progress as a function of the different
        *PLD* orders.

        :param ax: The current :py:obj:`matplotlib.pyplot` axis instance
        :param str info_left: Information to display at the left of the \
               plot. Default `''`
        :param str info_right: Information to display at the right of the \
               plot. Default `''`
        :param str color: The color of the data points. Default `'b'`

        '''

        # Plot
        if (self.cadence == 'lc') or (len(self.time) < 4000):
            ax.plot(self.apply_mask(self.time), self.apply_mask(self.flux),
                    ls='none', marker='.', color=color,
                    markersize=2, alpha=0.5)
            ax.plot(self.time[self.transitmask], self.flux[self.transitmask],
                    ls='none', marker='.', color=color,
                    markersize=2, alpha=0.5)
        else:
            ax.plot(self.apply_mask(self.time), self.apply_mask(
                    self.flux), ls='none', marker='.', color=color,
                    markersize=2, alpha=0.03, zorder=-1)
            ax.plot(self.time[self.transitmask], self.flux[self.transitmask],
                    ls='none', marker='.', color=color,
                    markersize=2, alpha=0.03, zorder=-1)
            ax.set_rasterization_zorder(0)
        ylim = self.get_ylim()

        # Plot the outliers, but not the NaNs
        badmask = [i for i in self.badmask if i not in self.nanmask]

        def O1(x): return x[self.outmask]

        def O2(x): return x[badmask]
        if self.cadence == 'lc':
            ax.plot(O1(self.time), O1(self.flux), ls='none',
                    color="#777777", marker='.', markersize=2, alpha=0.5)
            ax.plot(O2(self.time), O2(self.flux),
                    'r.', markersize=2, alpha=0.25)
        else:
            ax.plot(O1(self.time), O1(self.flux), ls='none', color="#777777",
                    marker='.', markersize=2, alpha=0.25, zorder=-1)
            ax.plot(O2(self.time), O2(self.flux), 'r.',
                    markersize=2, alpha=0.125, zorder=-1)
        for i in np.where(self.flux < ylim[0])[0]:
            if i in badmask:
                color = "#ffcccc"
            elif i in self.outmask:
                color = "#cccccc"
            elif i in self.nanmask:
                continue
            else:
                color = "#ccccff"
            ax.annotate('', xy=(self.time[i], ylim[0]), xycoords='data',
                        xytext=(0, 15), textcoords='offset points',
                        arrowprops=dict(arrowstyle="-|>", color=color))
        for i in np.where(self.flux > ylim[1])[0]:
            if i in badmask:
                color = "#ffcccc"
            elif i in self.outmask:
                color = "#cccccc"
            elif i in self.nanmask:
                continue
            else:
                color = "#ccccff"
            ax.annotate('', xy=(self.time[i], ylim[1]), xycoords='data',
                        xytext=(0, -15), textcoords='offset points',
                        arrowprops=dict(arrowstyle="-|>", color=color))

        # Plot the breakpoints
        for brkpt in self.breakpoints[:-1]:
            if len(self.breakpoints) <= 5:
                ax.axvline(self.time[brkpt], color='r', ls='--', alpha=0.5)
            else:
                ax.axvline(self.time[brkpt], color='r', ls='-', alpha=0.025)

        # Appearance
        if len(self.cdpp_arr) == 2:
            ax.annotate('%.2f ppm' % self.cdpp_arr[0], xy=(0.02, 0.975),
                        xycoords='axes fraction',
                        ha='left', va='top', fontsize=10)
            ax.annotate('%.2f ppm' % self.cdpp_arr[1], xy=(0.98, 0.975),
                        xycoords='axes fraction',
                        ha='right', va='top', fontsize=10)
        elif len(self.cdpp_arr) < 6:
            for n in range(len(self.cdpp_arr)):
                if n > 0:
                    x = (self.time[self.breakpoints[n - 1]] - self.time[0]
                         ) / (self.time[-1] - self.time[0]) + 0.02
                else:
                    x = 0.02
                ax.annotate('%.2f ppm' % self.cdpp_arr[n], xy=(x, 0.975),
                            xycoords='axes fraction',
                            ha='left', va='top', fontsize=8)
        else:
            ax.annotate('%.2f ppm' % self.cdpp, xy=(0.02, 0.975),
                        xycoords='axes fraction',
                        ha='left', va='top', fontsize=10)
        ax.annotate(info_right, xy=(0.98, 0.025), xycoords='axes fraction',
                    ha='right', va='bottom', fontsize=10, alpha=0.5,
                    fontweight='bold')
        ax.annotate(info_left, xy=(0.02, 0.025), xycoords='axes fraction',
                    ha='left', va='bottom', fontsize=8)
        ax.set_xlabel(r'Time (%s)' % self._mission.TIMEUNITS, fontsize=5)
        ax.margins(0.01, 0.1)
        ax.set_ylim(*ylim)
        ax.get_yaxis().set_major_formatter(Formatter.Flux)

    def plot_final(self, ax):
        '''
        Plots the final de-trended light curve.

        '''

        # Plot the light curve
        bnmask = np.array(
            list(set(np.concatenate([self.badmask, self.nanmask]))), dtype=int)

        def M(x): return np.delete(x, bnmask)
        if (self.cadence == 'lc') or (len(self.time) < 4000):
            ax.plot(M(self.time), M(self.flux), ls='none',
                    marker='.', color='k', markersize=2, alpha=0.3)
        else:
            ax.plot(M(self.time), M(self.flux), ls='none', marker='.',
                    color='k', markersize=2, alpha=0.03, zorder=-1)
            ax.set_rasterization_zorder(0)
        # Hack: Plot invisible first and last points to ensure
        # the x axis limits are the
        # same in the other plots, where we also plot outliers!
        ax.plot(self.time[0], np.nanmedian(M(self.flux)), marker='.', alpha=0)
        ax.plot(self.time[-1], np.nanmedian(M(self.flux)), marker='.', alpha=0)

        # Plot the GP (long cadence only)
        if self.cadence == 'lc':
            gp = GP(self.kernel, self.kernel_params, white=False)
            gp.compute(self.apply_mask(self.time),
                       self.apply_mask(self.fraw_err))
            med = np.nanmedian(self.apply_mask(self.flux))
            y, _ = gp.predict(self.apply_mask(self.flux) - med, self.time)
            y += med
            ax.plot(M(self.time), M(y), 'r-', lw=0.5, alpha=0.5)

            # Compute the CDPP of the GP-detrended flux
            self.cdppg = self._mission.CDPP(self.apply_mask(
                self.flux - y + med), cadence=self.cadence)

        else:

            # We're not going to calculate this
            self.cdppg = 0.

        # Appearance
        ax.annotate('Final', xy=(0.98, 0.025), xycoords='axes fraction',
                    ha='right', va='bottom', fontsize=10, alpha=0.5,
                    fontweight='bold')
        ax.margins(0.01, 0.1)

        # Get y lims that bound 99% of the flux
        flux = np.delete(self.flux, bnmask)
        N = int(0.995 * len(flux))
        hi, lo = flux[np.argsort(flux)][[N, -N]]
        fsort = flux[np.argsort(flux)]
        pad = (hi - lo) * 0.1
        ylim = (lo - pad, hi + pad)
        ax.set_ylim(ylim)
        ax.get_yaxis().set_major_formatter(Formatter.Flux)

    def plot_cbv(self, ax, flux, info, show_cbv=False):
        '''
        Plots the final CBV-corrected light curve.

        '''

        # Plot the light curve
        bnmask = np.array(
            list(set(np.concatenate([self.badmask, self.nanmask]))), dtype=int)

        def M(x): return np.delete(x, bnmask)
        if self.cadence == 'lc':
            ax.plot(M(self.time), M(flux), ls='none', marker='.',
                    color='k', markersize=2, alpha=0.45)
        else:
            ax.plot(M(self.time), M(flux), ls='none', marker='.',
                    color='k', markersize=2, alpha=0.03, zorder=-1)
            ax.set_rasterization_zorder(0)
        # Hack: Plot invisible first and last points to ensure
        # the x axis limits are the
        # same in the other plots, where we also plot outliers!
        ax.plot(self.time[0], np.nanmedian(M(flux)), marker='.', alpha=0)
        ax.plot(self.time[-1], np.nanmedian(M(flux)), marker='.', alpha=0)

        # Show CBV fit?
        if show_cbv:
            ax.plot(self.time, self._mission.FitCBVs(
                self) + np.nanmedian(flux), 'r-', alpha=0.2)

        # Appearance
        ax.annotate(info, xy=(0.98, 0.025), xycoords='axes fraction',
                    ha='right', va='bottom', fontsize=10, alpha=0.5,
                    fontweight='bold')
        ax.margins(0.01, 0.1)

        # Get y lims that bound 99% of the flux
        flux = np.delete(flux, bnmask)
        N = int(0.995 * len(flux))
        hi, lo = flux[np.argsort(flux)][[N, -N]]
        fsort = flux[np.argsort(flux)]
        pad = (hi - lo) * 0.2
        ylim = (lo - pad, hi + pad)
        ax.set_ylim(ylim)
        ax.get_yaxis().set_major_formatter(Formatter.Flux)
        ax.set_xlabel(r'Time (%s)' % self._mission.TIMEUNITS, fontsize=9)
        for tick in ax.get_xticklabels() + ax.get_yticklabels():
            tick.set_fontsize(7)

    def load_tpf(self):
        '''
        Loads the target pixel file.

        '''

        if not self.loaded:
            if self._data is not None:
                data = self._data
            else:
                data = self._mission.GetData(
                         self.ID, season=self.season,
                         cadence=self.cadence,
                         clobber=self.clobber_tpf,
                         aperture_name=self.aperture_name,
                         saturated_aperture_name=self.saturated_aperture_name,
                         max_pixels=self.max_pixels,
                         saturation_tolerance=self.saturation_tolerance,
                         get_hires=self.get_hires,
                         get_nearby=self.get_nearby)
                if data is None:
                    raise Exception("Unable to retrieve target data.")
            self.cadn = data.cadn
            self.time = data.time
            self.model = np.zeros_like(self.time)
            self.fpix = data.fpix
            self.fraw = np.sum(self.fpix, axis=1)
            self.fpix_err = data.fpix_err
            self.fraw_err = np.sqrt(np.sum(self.fpix_err ** 2, axis=1))
            self.nanmask = data.nanmask
            self.badmask = data.badmask
            self.transitmask = np.array([], dtype=int)
            self.outmask = np.array([], dtype=int)
            self.aperture = data.aperture
            self.aperture_name = data.aperture_name
            self.apertures = data.apertures
            self.quality = data.quality
            self.Xpos = data.Xpos
            self.Ypos = data.Ypos
            self.mag = data.mag
            self.pixel_images = data.pixel_images
            self.nearby = data.nearby
            self.hires = data.hires
            self.saturated = data.saturated
            self.meta = data.meta
            self.bkg = data.bkg

            # Update the last breakpoint to the correct value
            self.breakpoints[-1] = len(self.time) - 1

            # Get PLD normalization
            self.get_norm()

            self.loaded = True

    def load_model(self, name=None):
        '''
        Loads a saved version of the model.

        '''

        if self.clobber:
            return False

        if name is None:
            name = self.name
        file = os.path.join(self.dir, '%s.npz' % name)
        if os.path.exists(file):
            if not self.is_parent:
                log.info("Loading '%s.npz'..." % name)
            try:
                data = np.load(file)
                for key in data.keys():
                    try:
                        setattr(self, key, data[key][()])
                    except NotImplementedError:
                        pass

                # HACK: Backwards compatibility. Previous version stored
                # the CDPP in the `cdpp6`
                # and `cdpp6_arr` attributes. Let's move them over.
                if hasattr(self, 'cdpp6'):
                    self.cdpp = self.cdpp6
                    del self.cdpp6
                if hasattr(self, 'cdpp6_arr'):
                    self.cdpp_arr = np.array(self.cdpp6_arr)
                    del self.cdpp6_arr
                if hasattr(self, 'gppp'):
                    self.cdppg = self.gppp
                    del self.gppp

                # HACK: At one point we were saving the figure instances,
                # so loading the .npz
                # opened a plotting window. I don't think this is the case
                # any more, so this
                # next line should be removed in the future...
                pl.close()

                return True
            except:
                log.warn("Error loading '%s.npz'." % name)
                exctype, value, tb = sys.exc_info()
                for line in traceback.format_exception_only(exctype, value):
                    ln = line.replace('\n', '')
                    log.warn(ln)
                os.rename(file, file + '.bad')

        if self.is_parent:
            raise Exception(
                'Unable to load `%s` model for target %d.'
                % (self.name, self.ID))

        return False

    def save_model(self):
        '''
        Saves all of the de-trending information to disk in an `npz` file
        and saves the DVS as a `pdf`.

        '''

        # Save the data
        log.info("Saving data to '%s.npz'..." % self.name)
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
        d.pop('transit_model', None)
        d.pop('_transit_model', None)
        np.savez(os.path.join(self.dir, self.name + '.npz'), **d)

        # Save the DVS
        pdf = PdfPages(os.path.join(self.dir, self.name + '.pdf'))
        pdf.savefig(self.dvs.fig)
        pl.close(self.dvs.fig)
        d = pdf.infodict()
        d['Title'] = 'EVEREST: %s de-trending of %s %d' % (
            self.name, self._mission.IDSTRING, self.ID)
        d['Author'] = 'Rodrigo Luger'
        pdf.close()

    def exception_handler(self, pdb):
        '''
        A custom exception handler.

        :param pdb: If :py:obj:`True`, enters PDB post-mortem \
               mode for debugging.

        '''

        # Grab the exception
        exctype, value, tb = sys.exc_info()

        # Log the error and create a .err file
        errfile = os.path.join(self.dir, self.name + '.err')
        with open(errfile, 'w') as f:
            for line in traceback.format_exception_only(exctype, value):
                ln = line.replace('\n', '')
                log.error(ln)
                print(ln, file=f)
            for line in traceback.format_tb(tb):
                ln = line.replace('\n', '')
                log.error(ln)
                print(ln, file=f)

        # Re-raise?
        if pdb:
            raise

    def update_gp(self):
        '''
        Calls :py:func:`gp.GetKernelParams` to optimize the GP and obtain the
        covariance matrix for the regression.

        '''

        self.kernel_params = GetKernelParams(self.time, self.flux,
                                             self.fraw_err,
                                             mask=self.mask,
                                             guess=self.kernel_params,
                                             kernel=self.kernel,
                                             giter=self.giter,
                                             gmaxf=self.gmaxf)

    def init_kernel(self):
        '''
        Initializes the covariance matrix with a guess at
        the GP kernel parameters.

        '''

        if self.kernel_params is None:
            X = self.apply_mask(self.fpix / self.flux.reshape(-1, 1))
            y = self.apply_mask(self.flux) - np.dot(X, np.linalg.solve(
                np.dot(X.T, X), np.dot(X.T, self.apply_mask(self.flux))))
            white = np.nanmedian([np.nanstd(c) for c in Chunks(y, 13)])
            amp = self.gp_factor * np.nanstd(y)
            tau = 30.0
            if self.kernel == 'Basic':
                self.kernel_params = [white, amp, tau]
            elif self.kernel == 'QuasiPeriodic':
                self.kernel_params = [white, amp, 1., 20.]

    def mask_planets(self):
        '''

        '''

        for i, planet in enumerate(self.planets):
            log.info('Masking planet #%d...' % (i + 1))
            t0, period, dur = planet
            mask = []
            t0 += np.ceil((self.time[0] - dur - t0) / period) * period
            for t in np.arange(t0, self.time[-1] + dur, period):
                mask.extend(np.where(np.abs(self.time - t) < dur / 2.)[0])
            self.transitmask = np.array(
                list(set(np.concatenate([self.transitmask, mask]))))

    def run(self):
        '''
        Runs the de-trending step.

        '''

        try:

            # Load raw data
            log.info("Loading target data...")
            self.load_tpf()
            self.mask_planets()
            self.plot_aperture([self.dvs.top_right() for i in range(4)])
            self.init_kernel()
            M = self.apply_mask(np.arange(len(self.time)))
            self.cdppr_arr = self.get_cdpp_arr()
            self.cdpp_arr = np.array(self.cdppr_arr)
            self.cdppv_arr = np.array(self.cdppr_arr)
            self.cdppr = self.get_cdpp()
            self.cdpp = self.cdppr
            self.cdppv = self.cdppr

            log.info("%s (Raw): CDPP = %s" % (self.name, self.cdpps))
            self.plot_lc(self.dvs.left(), info_right='Raw', color='k')

            # Loop
            for n in range(self.pld_order):
                self.lam_idx += 1
                self.get_outliers()
                if n > 0 and self.optimize_gp:
                    self.update_gp()
                self.cross_validate(self.dvs.right(), info='CV%d' % n)
                self.cdpp_arr = self.get_cdpp_arr()
                self.cdppv_arr *= self.cdpp_arr
                self.cdpp = self.get_cdpp()
                self.cdppv = np.nanmean(self.cdppv_arr)
                log.info("%s (%d/%d): CDPP = %s" %
                         (self.name, n + 1, self.pld_order, self.cdpps))
                self.plot_lc(self.dvs.left(), info_right='LC%d' % (
                    n + 1), info_left='%d outliers' % len(self.outmask))

            # Save
            self.finalize()
            self.plot_final(self.dvs.top_left())
            self.plot_info(self.dvs)
            self.save_model()

        except:

            self.exception_handler(self.debug)

    def publish(self, **kwargs):
        '''
        Correct the light curve with the CBVs, generate a
        cover page for the DVS figure,
        and produce a FITS file for publication.

        '''

        try:

            # HACK: Force these params for publication
            self.cbv_win = 999
            self.cbv_order = 3
            self.cbv_num = 1

            # Get the CBVs
            self._mission.GetTargetCBVs(self)

            # Plot the final corrected light curve
            cbv = CBV()
            self.plot_info(cbv)
            self.plot_cbv(cbv.body(), self.fcor, 'Corrected')
            self.plot_cbv(cbv.body(), self.flux, 'De-trended', show_cbv=True)
            self.plot_cbv(cbv.body(), self.fraw, 'Raw')

            # Save the CBV pdf
            pdf = PdfPages(os.path.join(self.dir, 'cbv.pdf'))
            pdf.savefig(cbv.fig)
            pl.close(cbv.fig)
            d = pdf.infodict()
            d['Title'] = 'EVEREST: %s de-trending of %s %d' % (
                self.name, self._mission.IDSTRING, self.ID)
            d['Author'] = 'Rodrigo Luger'
            pdf.close()

            # Now merge the two PDFs
            assert os.path.exists(os.path.join(
                self.dir, self.name + '.pdf')), \
                "Unable to locate %s.pdf." % self.name
            output = PdfFileWriter()
            pdfOne = PdfFileReader(os.path.join(self.dir, 'cbv.pdf'))
            pdfTwo = PdfFileReader(os.path.join(self.dir, self.name + '.pdf'))
            # Add the CBV page
            output.addPage(pdfOne.getPage(0))
            # Add the original DVS page
            output.addPage(pdfTwo.getPage(pdfTwo.numPages - 1))
            # Write the final PDF
            outputStream = open(os.path.join(self.dir, self._mission.DVSFile(
                self.ID, self.season, self.cadence)), "wb")
            output.write(outputStream)
            outputStream.close()
            os.remove(os.path.join(self.dir, 'cbv.pdf'))

            # Make the FITS file
            MakeFITS(self)

        except:

            self.exception_handler(self.debug)

    def publish_csv(self, **kwargs):
        '''


        '''

        try:

            # HACK: Force these params for publication
            self.cbv_win = 999
            self.cbv_order = 3
            self.cbv_num = 1

            # Get the CBVs
            self._mission.GetTargetCBVs(self)

            # Write to file!
            outfile = os.path.join(self.dir, self._mission.CSVFile(self.ID))
            header = self._mission.CSVHEADER % self.ID
            mask = np.zeros_like(self.cadn)
            for i in range(len(mask)):
                if i in self.nanmask:
                    mask[i] = 1
                elif i in self.badmask:
                    mask[i] = 2
                elif i in self.outmask:
                    mask[i] = 3
            data = np.vstack([self.time, self.cadn, self.fcor,
                              self.flux, self.fraw, mask]).T
            np.savetxt(outfile, data,
                       fmt='%.6f,%d,%.6f,%.6f,%.6f,%d', header=header)

        except:

            self.exception_handler(self.debug)


class rPLD(Detrender):
    '''
    The regular PLD model. Nothing fancy.

    '''

    pass


class nPLD(Detrender):
    '''
    The "neighboring stars" *PLD* model. This model uses the
    *PLD* vectors of neighboring stars to help in the de-trending and can lead
    to increased performance over the regular :py:class:`rPLD` model,
    particularly for dimmer stars.

    '''

    def setup(self, **kwargs):
        '''
        This is called during production de-trending, prior to
        calling the :py:obj:`Detrender.run()` method.

        :param tuple cdpp_range: If :py:obj:`parent_model` is set, \
               neighbors are selected only if \
               their de-trended CDPPs fall within this range. Default `None`
        :param tuple mag_range: Only select neighbors whose magnitudes are \
               within this range. Default (11., 13.)
        :param int neighbors: The number of neighboring stars to use in \
               the de-trending. The higher this number, the more signals \
               there are and hence the more de-trending information there is. \
               However, the neighboring star signals are regularized together \
               with the target's signals, so adding too many neighbors will \
               inevitably reduce the contribution of the target's own \
               signals, which may reduce performance. Default `10`
        :param str parent_model: By default, :py:class:`nPLD` is run in \
               stand-alone mode. The neighbor signals are computed directly \
               from their TPFs, so there is no need to have run *PLD* on them \
               beforehand. However, if :py:obj:`parent_model` \
               is set, :py:class:`nPLD` will use information from the \
               :py:obj:`parent_model` model of each neighboring star when \
               de-trending. This is particularly useful for identifying \
               outliers in the neighbor signals and preventing them from \
               polluting the current target. Setting :py:obj:`parent_model` \
               to :py:class:`rPLD`, for instance, will use the \
               outlier information in the :py:class:`rPLD` model of the \
               neighbors (this must have been run ahead of time). \
               Note, however, that tests with *K2* data show that including \
               outliers in the neighbor signals actually \
               *improves* the performance, since many of these outliers \
               are associated with events such as thruster firings and are \
               present in all light curves, and therefore *help* in the \
               de-trending. Default `None`

        ..note :: Optionally, the :py:obj:`neighbors` may be specified \
                  directly as a list of target IDs to use. \
                  In this case, users may also provide a list of \
                  :py:class:`everest.utils.DataContainer` instances \
                  corresponding to each of the neighbors in the \
                  :py:obj:`neighbors_data` kwarg.
        '''

        # Get neighbors
        self.parent_model = kwargs.get('parent_model', None)
        neighbors = kwargs.get('neighbors', 10)
        neighbors_data = kwargs.get('neighbors_data', None)
        if hasattr(neighbors, '__len__'):
            self.neighbors = neighbors
        else:
            num_neighbors = neighbors
            self.neighbors = \
                self._mission.GetNeighbors(self.ID,
                                           season=self.season,
                                           cadence=self.cadence,
                                           model=self.parent_model,
                                           neighbors=num_neighbors,
                                           mag_range=kwargs.get(
                                               'mag_range', (11., 13.)),
                                           cdpp_range=kwargs.get(
                                               'cdpp_range', None),
                                           aperture_name=self.aperture_name)
            if len(self.neighbors):
                if len(self.neighbors) < num_neighbors:
                    log.warn("%d neighbors requested, but only %d found." %
                             (num_neighbors, len(self.neighbors)))
            elif num_neighbors > 0:
                log.warn("No neighbors found! Running standard PLD...")

        for n, neighbor in enumerate(self.neighbors):
            log.info("Loading data for neighboring target %d..." % neighbor)
            if neighbors_data is not None:
                data = neighbors_data[n]
                data.mask = np.array(
                    list(set(np.concatenate([data.badmask, data.nanmask]))),
                    dtype=int)
                data.fraw = np.sum(data.fpix, axis=1)
            elif self.parent_model is not None and self.cadence == 'lc':
                # We load the `parent` model. The advantage here is
                # that outliers have properly been identified and masked.
                # I haven't tested this on short
                # cadence data, so I'm going to just forbid it...
                data = eval(self.parent_model)(
                    neighbor, mission=self.mission, is_parent=True)
            else:
                # We load the data straight from the TPF. Much quicker,
                # since no model must be run in advance. Downside is we
                # don't know where the outliers are. But based
                # on tests with K2 data, the de-trending is actually
                # *better* if the outliers are
                # included! These are mostly thruster fire events and other
                #  artifacts common to
                # all the stars, so it makes sense that we might want
                # to keep them in the design matrix.
                data = self._mission.GetData(neighbor, season=self.season,
                                             clobber=self.clobber_tpf,
                                             cadence=self.cadence,
                                             aperture_name=self.aperture_name,
                                             saturated_aperture_name=
                                             self.saturated_aperture_name,
                                             max_pixels=self.max_pixels,
                                             saturation_tolerance=
                                             self.saturation_tolerance,
                                             get_hires=False, get_nearby=False)
                if data is None:
                    raise Exception(
                        "Unable to retrieve data for neighboring target.")
                data.mask = np.array(
                    list(set(np.concatenate([data.badmask, data.nanmask]))),
                    dtype=int)
                data.fraw = np.sum(data.fpix, axis=1)

            # Compute the linear PLD vectors and interpolate over
            # outliers, NaNs and bad timestamps
            X1 = data.fpix / data.fraw.reshape(-1, 1)
            X1 = Interpolate(data.time, data.mask, X1)
            if self.X1N is None:
                self.X1N = np.array(X1)
            else:
                self.X1N = np.hstack([self.X1N, X1])
            del X1
            del data


class iPLD(Detrender):
    '''
    The iterative PLD model.

    ..warning :: Deprecated and not thoroughly tested.

    '''

    def setup(self, **kwargs):
        '''
        This is called during production de-trending, prior to
        calling the :py:obj:`Detrender.run()` method.

        :param str parent_model: The name of the model to operate on. \
               Default `nPLD`

        '''

        # Load the parent model
        self.parent_model = kwargs.get('parent_model', 'nPLD')
        if not self.load_model(self.parent_model):
            raise Exception('Unable to load parent model.')

        # Save static copies of the de-trended flux,
        # the outlier mask and the lambda array
        self._norm = np.array(self.flux)
        self.recmask = np.array(self.mask)
        self.reclam = np.array(self.lam)

        # Now reset the model params
        self.optimize_gp = False
        nseg = len(self.breakpoints)
        self.lam_idx = -1
        self.lam = [
            [1e5] + [None for i in range(self.pld_order - 1)]
            for b in range(nseg)]
        self.cdpp_arr = np.array([np.nan for b in range(nseg)])
        self.cdppr_arr = np.array([np.nan for b in range(nseg)])
        self.cdppv_arr = np.array([np.nan for b in range(nseg)])
        self.cdpp = np.nan
        self.cdppr = np.nan
        self.cdppv = np.nan
        self.cdppg = np.nan
        self.model = np.zeros_like(self.time)
        self.loaded = True


class pPLD(Detrender):
    '''
    A neighboring PLD extension that uses Powell's method to find the
    cross-validation parameter :py:obj:`lambda`.

    '''

    def setup(self, **kwargs):
        '''
        This is called during production de-trending, prior to
        calling the :py:obj:`Detrender.run()` method.

        :param inter piter: The number of iterations in the minimizer. \
               Default 3
        :param int pmaxf: The maximum number of function evaluations per \
               iteration. Default 300
        :param float ppert: The fractional amplitude of the perturbation on \
               the initial guess. Default 0.1

        '''

        # Check for saved model
        clobber = self.clobber
        self.clobber = False
        if not self.load_model('nPLD'):
            raise Exception("Can't find `nPLD` model for target.")
        self.clobber = clobber

        # Powell iterations
        self.piter = kwargs.get('piter', 3)
        self.pmaxf = kwargs.get('pmaxf', 300)
        self.ppert = kwargs.get('ppert', 0.1)

    def run(self):
        '''
        Runs the de-trending.

        '''

        try:

            # Plot original
            self.plot_aperture([self.dvs.top_right() for i in range(4)])
            self.plot_lc(self.dvs.left(), info_right='nPLD', color='k')

            # Cross-validate
            self.cross_validate(self.dvs.right())
            self.compute()
            self.cdpp_arr = self.get_cdpp_arr()
            self.cdpp = self.get_cdpp()

            # Plot new
            self.plot_lc(self.dvs.left(), info_right='Powell', color='k')

            # Save
            self.plot_final(self.dvs.top_left())
            self.plot_info(self.dvs)
            self.save_model()

        except:

            self.exception_handler(self.debug)

    def cross_validate(self, ax):
        '''
        Performs the cross-validation step.

        '''

        # The CDPP to beat
        cdpp_opt = self.get_cdpp_arr()

        # Loop over all chunks
        for b, brkpt in enumerate(self.breakpoints):

            log.info("Cross-validating chunk %d/%d..." %
                     (b + 1, len(self.breakpoints)))

            # Mask for current chunk
            m = self.get_masked_chunk(b)

            # Mask transits and outliers
            time = self.time[m]
            flux = self.fraw[m]
            ferr = self.fraw_err[m]
            med = np.nanmedian(self.fraw)

            # Setup the GP
            gp = GP(self.kernel, self.kernel_params, white=False)
            gp.compute(time, ferr)

            # The masks
            masks = list(Chunks(np.arange(0, len(time)),
                                len(time) // self.cdivs))

            # The pre-computed matrices
            pre_v = [self.cv_precompute(mask, b) for mask in masks]

            # Initialize with the nPLD solution
            log_lam_opt = np.log10(self.lam[b])
            scatter_opt = self.validation_scatter(
                log_lam_opt, b, masks, pre_v, gp, flux, time, med)
            log.info("Iter 0/%d: " % (self.piter) +
                     "logL = (%s), s = %.3f" %
                     (", ".join(["%.3f" % l for l in log_lam_opt]),
                      scatter_opt))

            # Do `piter` iterations
            for p in range(self.piter):

                # Perturb the initial condition a bit
                log_lam = np.array(
                    np.log10(self.lam[b])) * \
                    (1 + self.ppert * np.random.randn(len(self.lam[b])))
                scatter = self.validation_scatter(
                    log_lam, b, masks, pre_v, gp, flux, time, med)
                log.info("Initializing at: " +
                         "logL = (%s), s = %.3f" %
                         (", ".join(["%.3f" % l for l in log_lam]), scatter))

                # Call the minimizer
                log_lam, scatter, _, _, _, _ = \
                    fmin_powell(self.validation_scatter, log_lam,
                                args=(b, masks, pre_v, gp, flux, time, med),
                                maxfun=self.pmaxf, disp=False,
                                full_output=True)

                # Did it improve the CDPP?
                tmp = np.array(self.lam[b])
                self.lam[b] = 10 ** log_lam
                self.compute()
                cdpp = self.get_cdpp_arr()[b]
                self.lam[b] = tmp
                if cdpp < cdpp_opt[b]:
                    cdpp_opt[b] = cdpp
                    log_lam_opt = log_lam

                # Log it
                log.info("Iter %d/%d: " % (p + 1, self.piter) +
                         "logL = (%s), s = %.3f" %
                         (", ".join(["%.3f" % l for l in log_lam]), scatter))

            # The best solution
            log.info("Found minimum: logL = (%s), s = %.3f" %
                     (", ".join(["%.3f" % l for l in log_lam_opt]),
                      scatter_opt))
            self.lam[b] = 10 ** log_lam_opt

        # We're just going to plot lambda as a function of chunk number
        bs = np.arange(len(self.breakpoints))
        color = ['k', 'b', 'r', 'g', 'y']
        for n in range(self.pld_order):
            ax[0].plot(bs + 1, [np.log10(self.lam[b][n])
                                for b in bs], '.', color=color[n])
            ax[0].plot(bs + 1, [np.log10(self.lam[b][n])
                                for b in bs], '-', color=color[n], alpha=0.25)
        ax[0].set_ylabel(r'$\log\Lambda$', fontsize=5)
        ax[0].margins(0.1, 0.1)
        ax[0].set_xticks(np.arange(1, len(self.breakpoints) + 1))
        ax[0].set_xticklabels([])

        # Now plot the CDPP
        cdpp_arr = self.get_cdpp_arr()
        ax[1].plot(bs + 1, cdpp_arr, 'b.')
        ax[1].plot(bs + 1, cdpp_arr, 'b-', alpha=0.25)
        ax[1].margins(0.1, 0.1)
        ax[1].set_ylabel(r'Scatter (ppm)', fontsize=5)
        ax[1].set_xlabel(r'Chunk', fontsize=5)
        ax[1].set_xticks(np.arange(1, len(self.breakpoints) + 1))

    def validation_scatter(self, log_lam, b, masks, pre_v, gp, flux,
                           time, med):
        '''
        Computes the scatter in the validation set.

        '''

        # Update the lambda matrix
        self.lam[b] = 10 ** log_lam

        # Validation set scatter
        scatter = [None for i in range(len(masks))]
        for i in range(len(masks)):
            model = self.cv_compute(b, *pre_v[i])
            try:
                gpm, _ = gp.predict(flux - model - med, time[masks[i]])
            except ValueError:
                # Sometimes the model can have NaNs if
                # `lambda` is a crazy value
                return 1.e30
            fdet = (flux - model)[masks[i]] - gpm
            scatter[i] = 1.e6 * (1.4826 * np.nanmedian(np.abs(fdet / med -
                                 np.nanmedian(fdet / med))) /
                                 np.sqrt(len(masks[i])))

        return np.max(scatter)
