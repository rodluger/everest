#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`mathutils.py` - Math utils
-----------------------------------

Miscellaneous math utilities used throughout the code.

'''

from __future__ import division, print_function, absolute_import, \
     unicode_literals
import numpy as np
from scipy.signal import medfilt
from scipy.signal import savgol_filter
from scipy.special import comb
import logging
log = logging.getLogger(__name__)


def Interpolate(time, mask, y):
    '''
    Masks certain elements in the array `y` and linearly
    interpolates over them, returning an array `y'` of the
    same length.

    :param array_like time: The time array
    :param array_like mask: The indices to be interpolated over
    :param array_like y: The dependent array

    '''

    # Ensure `y` doesn't get modified in place
    yy = np.array(y)
    t_ = np.delete(time, mask)
    y_ = np.delete(y, mask, axis=0)
    if len(yy.shape) == 1:
        yy[mask] = np.interp(time[mask], t_, y_)
    elif len(yy.shape) == 2:
        for n in range(yy.shape[1]):
            yy[mask, n] = np.interp(time[mask], t_, y_[:, n])
    else:
        raise Exception("Array ``y`` must be either 1- or 2-d.")
    return yy


def MedianFilter(x, kernel_size=5):
    '''
    A silly wrapper around :py:func:`scipy.signal.medfilt`.

    '''

    if kernel_size % 2 == 0:
        kernel_size += 1
    return medfilt(x, kernel_size=kernel_size)


def Chunks(l, n, all=False):
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
                yield l[i:i + n]
            else:
                if not all:
                    yield l[i:]
                break


def Smooth(x, window_len=100, window='hanning'):
    '''
    Smooth data by convolving on a given timescale.

    :param ndarray x: The data array
    :param int window_len: The size of the smoothing window. Default `100`
    :param str window: The window type. Default `hanning`


    '''

    if window_len == 0:
        return np.zeros_like(x)
    s = np.r_[2 * x[0] - x[window_len - 1::-1],
              x, 2 * x[-1] - x[-1:-window_len:-1]]
    if window == 'flat':
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')
    y = np.convolve(w / w.sum(), s, mode='same')
    return y[window_len:-window_len + 1]


def Scatter(y, win=13, remove_outliers=False):
    '''
    Return the scatter in ppm based on the median running standard deviation
    for a window size of :py:obj:`win` = 13 cadences (for K2, this
    is ~6.5 hours, as in VJ14).

    :param ndarray y: The array whose CDPP is to be computed
    :param int win: The window size in cadences. Default `13`
    :param bool remove_outliers: Clip outliers at 5 sigma before computing \
           the CDPP? Default `False`

    '''

    if remove_outliers:
        # Remove 5-sigma outliers from data
        # smoothed on a 1 day timescale
        if len(y) >= 50:
            ys = y - Smooth(y, 50)
        else:
            ys = y
        M = np.nanmedian(ys)
        MAD = 1.4826 * np.nanmedian(np.abs(ys - M))
        out = []
        for i, _ in enumerate(y):
            if (ys[i] > M + 5 * MAD) or (ys[i] < M - 5 * MAD):
                out.append(i)
        out = np.array(out, dtype=int)
        y = np.delete(y, out)
    if len(y):
        return 1.e6 * np.nanmedian([np.std(yi) / np.sqrt(win)
                                    for yi in Chunks(y, win, all=True)])
    else:
        return np.nan


def SavGol(y, win=49):
    '''
    Subtracts a second order Savitsky-Golay filter with window size `win`
    and returns the result. This acts as a high pass filter.

    '''

    if len(y) >= win:
        return y - savgol_filter(y, win, 2) + np.nanmedian(y)
    else:
        return y


def NumRegressors(npix, pld_order, cross_terms=True):
    '''
    Return the number of regressors for `npix` pixels
    and PLD order `pld_order`.

    :param bool cross_terms: Include pixel cross-terms? Default :py:obj:`True`

    '''

    res = 0
    for k in range(1, pld_order + 1):
        if cross_terms:
            res += comb(npix + k - 1, k)
        else:
            res += npix
    return int(res)


def Downbin(x, newsize, axis=0, operation='mean'):
    '''
    Downbins an array to a smaller size.

    :param array_like x: The array to down-bin
    :param int newsize: The new size of the axis along which to down-bin
    :param int axis: The axis to operate on. Default 0
    :param str operation: The operation to perform when down-binning. \
           Default `mean`
    '''

    assert newsize < x.shape[axis], \
        "The new size of the array must be smaller than the current size."
    oldsize = x.shape[axis]
    newshape = list(x.shape)
    newshape[axis] = newsize
    newshape.insert(axis + 1, oldsize // newsize)
    trim = oldsize % newsize
    if trim:
        xtrim = x[:-trim]
    else:
        xtrim = x

    if operation == 'mean':
        xbin = np.nanmean(xtrim.reshape(newshape), axis=axis + 1)
    elif operation == 'sum':
        xbin = np.nansum(xtrim.reshape(newshape), axis=axis + 1)
    elif operation == 'quadsum':
        xbin = np.sqrt(np.nansum(xtrim.reshape(newshape) ** 2, axis=axis + 1))
    elif operation == 'median':
        xbin = np.nanmedian(xtrim.reshape(newshape), axis=axis + 1)
    else:
        raise ValueError("`operation` must be either `mean`, " +
                         "`sum`, `quadsum`, or `median`.")

    return xbin
