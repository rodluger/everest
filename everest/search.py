#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
search.py
---------

Given an `everest` instance, performs a transit search
and returns an array of delta chi-squared values. This
approach is similar to that in Foreman-Mackey et al. (2015).

'''

from __future__ import division, print_function, absolute_import, \
     unicode_literals
import numpy as np
from .mathutils import SavGol
from .gp import GetCovariance
from .transit import TransitShape
from scipy.linalg import cho_solve, cho_factor
try:
    from tqdm import tqdm
    def prange(x): return tqdm(range(x))
except ImportError:
    def prange(x): return range(x)
import logging
log = logging.getLogger(__name__)


def Search(star, pos_tol=2.5, neg_tol=50., **ps_kwargs):
    '''
    NOTE: `pos_tol` is the positive (i.e., above the median)
    outlier tolerance in standard deviations.
    NOTE: `neg_tol` is the negative (i.e., below the median)
    outlier tolerance in standard deviations.

    '''

    # Smooth the light curve
    t = np.delete(star.time, np.concatenate([star.nanmask, star.badmask]))
    f = np.delete(star.flux, np.concatenate([star.nanmask, star.badmask]))
    f = SavGol(f)
    med = np.nanmedian(f)

    # Kill positive outliers
    MAD = 1.4826 * np.nanmedian(np.abs(f - med))
    pos_inds = np.where((f > med + pos_tol * MAD))[0]
    pos_inds = np.array([np.argmax(star.time == t[i]) for i in pos_inds])

    # Kill negative outliers
    MAD = 1.4826 * np.nanmedian(np.abs(f - med))
    neg_inds = np.where((f < med - neg_tol * MAD))[0]
    neg_inds = np.array([np.argmax(star.time == t[i]) for i in neg_inds])

    # Replace the star.outmask array
    star.outmask = np.concatenate([neg_inds, pos_inds])
    star.transitmask = np.array([], dtype=int)

    # Delta chi squared
    TIME = np.array([])
    DEPTH = np.array([])
    VARDEPTH = np.array([])
    DELCHISQ = np.array([])
    for b, brkpt in enumerate(star.breakpoints):

        # Log
        log.info('Running chunk %d/%d...' % (b + 1, len(star.breakpoints)))

        # Masks for current chunk
        m = star.get_masked_chunk(b, pad=False)

        # This block of the masked covariance matrix
        K = GetCovariance(star.kernel, star.kernel_params,
                          star.time[m], star.fraw_err[m])

        # The masked X.L.X^T term
        A = np.zeros((len(m), len(m)))
        for n in range(star.pld_order):
            XM = star.X(n, m)
            A += star.lam[b][n] * np.dot(XM, XM.T)
        K += A
        CDK = cho_factor(K)

        # Baseline
        med = np.nanmedian(star.fraw[m])
        lnL0 = -0.5 * np.dot(star.fraw[m], cho_solve(CDK, star.fraw[m]))
        dt = np.median(np.diff(star.time[m]))

        # Create a uniform time array and get indices of missing cadences
        tol = np.nanmedian(np.diff(star.time[m])) / 5.
        tunif = np.arange(star.time[m][0], star.time[m][-1] + tol, dt)
        tnogaps = np.array(tunif)
        gaps = []
        j = 0
        for i, t in enumerate(tunif):
            if np.abs(star.time[m][j] - t) < tol:
                tnogaps[i] = star.time[m][j]
                j += 1
                if j == len(star.time[m]):
                    break
            else:
                gaps.append(i)
        gaps = np.array(gaps, dtype=int)

        # Compute the normalized transit model for a single transit
        transit_model = TransitShape(**ps_kwargs)

        # Now roll the transit model across each cadence
        dchisq = np.zeros(len(tnogaps))
        d = np.zeros(len(tnogaps))
        vard = np.zeros(len(tnogaps))
        for i in prange(len(tnogaps)):
            trn = transit_model(tnogaps, tnogaps[i])
            trn = np.delete(trn, gaps)
            trn *= med
            vard[i] = 1. / np.dot(trn, cho_solve(CDK, trn))
            if not np.isfinite(vard[i]):
                vard[i] = np.nan
                d[i] = np.nan
                dchisq[i] = np.nan
                continue
            d[i] = vard[i] * np.dot(trn, cho_solve(CDK, star.fraw[m]))
            r = star.fraw[m] - trn * d[i]
            lnL = -0.5 * np.dot(r, cho_solve(CDK, r))
            dchisq[i] = -2 * (lnL0 - lnL)

        TIME = np.append(TIME, tnogaps)
        DEPTH = np.append(DEPTH, d)
        VARDEPTH = np.append(VARDEPTH, vard)
        DELCHISQ = np.append(DELCHISQ, dchisq)

    return TIME, DEPTH, VARDEPTH, DELCHISQ
