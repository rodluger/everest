#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`sysrem.py` - CBV routines
----------------------------------

Routines for computing the co-trending basis vectors (CBVs)
for each `K2` campaign using the :py:obj:`SysRem` algorithm.

'''

from __future__ import division, print_function, absolute_import, \
     unicode_literals
from ...config import EVEREST_DAT
from ...utils import InitLog
from .utils import GetK2Campaign, Campaign, Channels
import os
import numpy as np
import matplotlib.pyplot as pl
from scipy.signal import savgol_filter
import logging
log = logging.getLogger(__name__)


def GetChunk(time, breakpoints, b, mask=[]):
    '''
    Returns the indices corresponding to a given light curve chunk.

    :param int b: The index of the chunk to return

    '''

    M = np.delete(np.arange(len(time)), mask, axis=0)
    if b > 0:
        res = M[(M > breakpoints[b - 1]) & (M <= breakpoints[b])]
    else:
        res = M[M <= breakpoints[b]]
    return res


def GetStars(campaign, module, model='nPLD', **kwargs):
    '''
    Returns de-trended light curves for all stars on a given module in
    a given campaign.

    '''

    # Get the channel numbers
    channels = Channels(module)
    assert channels is not None, "No channels available on this module."

    # Get the EPIC numbers
    all = GetK2Campaign(campaign)
    stars = np.array([s[0] for s in all if s[2] in channels and
                      os.path.exists(
        os.path.join(EVEREST_DAT, 'k2', 'c%02d' % int(campaign),
                     ('%09d' % s[0])[:4] + '00000',
                     ('%09d' % s[0])[4:], model + '.npz'))], dtype=int)
    N = len(stars)
    assert N > 0, "No light curves found for campaign %d, module %d." % (
        campaign, module)

    # Loop over all stars and store the fluxes in a list
    fluxes = []
    errors = []
    kpars = []

    for n in range(N):

        # De-trended light curve file name
        nf = os.path.join(EVEREST_DAT, 'k2', 'c%02d' % int(campaign),
                          ('%09d' % stars[n])[:4] + '00000',
                          ('%09d' % stars[n])[4:], model + '.npz')

        # Get the data
        data = np.load(nf)
        t = data['time']
        if n == 0:
            time = t
            breakpoints = data['breakpoints']

        # Get de-trended light curve
        y = data['fraw'] - data['model']
        err = data['fraw_err']

        # De-weight outliers and bad timestamps
        m = np.array(list(set(np.concatenate([data['outmask'], data['badmask'],
                                              data['nanmask'],
                                              data['transitmask']]))),
                     dtype=int)

        # Interpolate over the outliers
        y = np.interp(t, np.delete(t, m), np.delete(y, m))
        err = np.interp(t, np.delete(t, m), np.delete(err, m))

        # Append to our running lists
        fluxes.append(y)
        errors.append(err)
        kpars.append(data['kernel_params'])

    return time, breakpoints, np.array(fluxes), \
           np.array(errors), np.array(kpars)


def SysRem(time, flux, err, ncbv=5, niter=50, sv_win=999,
           sv_order=3, **kwargs):
    '''
    Applies :py:obj:`SysRem` to a given set of light curves.

    :param array_like time: The time array for all of the light curves
    :param array_like flux: A 2D array of the fluxes for each of the light \
           curves, shape `(nfluxes, ntime)`
    :param array_like err: A 2D array of the flux errors for each of the \
           light curves, shape `(nfluxes, ntime)`
    :param int ncbv: The number of signals to recover. Default 5
    :param int niter: The number of :py:obj:`SysRem` iterations to perform. \
           Default 50
    :param int sv_win: The Savitsky-Golay filter window size. Default 999
    :param int sv_order: The Savitsky-Golay filter order. Default 3

    '''

    nflx, tlen = flux.shape

    # Get normalized fluxes
    med = np.nanmedian(flux, axis=1).reshape(-1, 1)
    y = flux - med

    # Compute the inverse of the variances
    invvar = 1. / err ** 2

    # The CBVs for this set of fluxes
    cbvs = np.zeros((ncbv, tlen))

    # Recover `ncbv` components
    for n in range(ncbv):

        # Initialize the weights and regressors
        c = np.zeros(nflx)
        a = np.ones(tlen)
        f = y * invvar

        # Perform `niter` iterations
        for i in range(niter):

            # Compute the `c` vector (the weights)
            c = np.dot(f, a) / np.dot(invvar, a ** 2)

            # Compute the `a` vector (the regressors)
            a = np.dot(c, f) / np.dot(c ** 2, invvar)

        # Remove this component from all light curves
        y -= np.outer(c, a)

        # Save this regressor after smoothing it a bit
        if sv_win >= len(a):
            sv_win = len(a) - 1
            if sv_win % 2 == 0:
                sv_win -= 1
        cbvs[n] = savgol_filter(a - np.nanmedian(a), sv_win, sv_order)

    return cbvs


def GetCBVs(campaign, model='nPLD', clobber=False, **kwargs):
    '''
    Computes the CBVs for a given campaign.

    :param int campaign: The campaign number
    :param str model: The name of the :py:obj:`everest` model. Default `nPLD`
    :param bool clobber: Overwrite existing files? Default `False`

    '''

    # Initialize logging?
    if len(logging.getLogger().handlers) == 0:
        InitLog(file_name=None, screen_level=logging.DEBUG)
    log.info('Computing CBVs for campaign %d...' % (campaign))

    # Output path
    path = os.path.join(EVEREST_DAT, 'k2', 'cbv', 'c%02d' % campaign)
    if not os.path.exists(path):
        os.makedirs(path)

    # Get the design matrix
    xfile = os.path.join(path, 'X.npz')
    if clobber or not os.path.exists(xfile):

        log.info('Obtaining light curves...')
        time = None
        for module in range(2, 25):

            # Get the light curves
            log.info('Loading module %d/%d...' % (module, 24))
            lcfile = os.path.join(path, '%d.npz' % module)
            if clobber or not os.path.exists(lcfile):
                try:
                    time, breakpoints, fluxes, errors, kpars = GetStars(
                        campaign, module, model=model, **kwargs)
                except AssertionError:
                    continue
                np.savez(lcfile, time=time, breakpoints=breakpoints,
                         fluxes=fluxes, errors=errors, kpars=kpars)

            # Load the light curves
            lcs = np.load(lcfile)
            if time is None:
                time = lcs['time']
                breakpoints = lcs['breakpoints']
                fluxes = lcs['fluxes']
                errors = lcs['errors']
                kpars = lcs['kpars']
            else:
                fluxes = np.vstack([fluxes, lcs['fluxes']])
                errors = np.vstack([errors, lcs['errors']])
                kpars = np.vstack([kpars, lcs['kpars']])

        # Compute the design matrix
        log.info('Running SysRem...')
        X = np.ones((len(time), 1 + kwargs.get('ncbv', 5)))

        # Loop over the segments
        new_fluxes = np.zeros_like(fluxes)
        for b in range(len(breakpoints)):

            # Get the current segment's indices
            inds = GetChunk(time, breakpoints, b)

            # Update the error arrays with the white GP component
            for j in range(len(errors)):
                errors[j] = np.sqrt(errors[j] ** 2 + kpars[j][0] ** 2)

            # Get de-trended fluxes
            X[inds, 1:] = SysRem(time[inds], fluxes[:, inds],
                                 errors[:, inds], **kwargs).T

        # Save
        np.savez(xfile, X=X, time=time, breakpoints=breakpoints)

    else:

        # Load from disk
        data = np.load(xfile)
        X = data['X'][()]
        time = data['time'][()]
        breakpoints = data['breakpoints'][()]

    # Plot
    plotfile = os.path.join(path, 'X.pdf')
    if clobber or not os.path.exists(plotfile):
        fig, ax = pl.subplots(2, 3, figsize=(12, 8))
        fig.subplots_adjust(left=0.05, right=0.95)
        ax = ax.flatten()
        for axis in ax:
            axis.set_xticks([])
            axis.set_yticks([])
        for b in range(len(breakpoints)):
            inds = GetChunk(time, breakpoints, b)
            for n in range(min(6, X.shape[1])):
                ax[n].plot(time[inds], X[inds, n])
                ax[n].set_title(n, fontsize=14)
        fig.savefig(plotfile, bbox_inches='tight')

    return X
