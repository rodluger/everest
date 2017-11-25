#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`overfit.py` - Overfitting metrics
------------------------------------------

EXPERIMENTAL!

'''

from __future__ import division, print_function, absolute_import, \
     unicode_literals
from .transit import TransitShape
from .utils import prange
from .masksolve import MaskSolve
import numpy as np
import matplotlib.pyplot as pl
import logging
log = logging.getLogger(__name__)


class Overfit(object):
    '''
    Generic overfitting metric class.

    '''

    def __init__(self, ID, IDSTRING, t, y, X, XL, XLX, C, K, mask=[],
                 tau=None, w=9):
        '''

        '''

        # Store the values
        self.ID = ID
        self.IDSTRING = IDSTRING
        self.t = t
        self.y = y
        self.X = X
        self.XL = XL
        self.XLX = XLX
        self.C = C
        self.K = K
        self.Kinv = np.linalg.inv(K)
        self.mask = mask
        self.w = w

        # The transit model
        if tau is None:
            self.tau = TransitShape(dur=0.1)
        else:
            self.tau = tau

        # Compute
        self._compute_masked()
        self._compute_unmasked()

    def plot_unmasked(self):
        '''

        '''

        # Set up
        fig = pl.figure(figsize=(8.5, 11))
        fig.subplots_adjust(hspace=0.25, left=0.15, right=0.9)
        axes = [pl.subplot2grid((3, 5), (0, 0), colspan=4, rowspan=1),
                pl.subplot2grid((3, 5), (1, 0), colspan=4, rowspan=1),
                pl.subplot2grid((3, 5), (2, 0), colspan=4, rowspan=1)]
        axesh = [pl.subplot2grid((3, 5), (0, 4), colspan=1, rowspan=1),
                 pl.subplot2grid((3, 5), (1, 4), colspan=1, rowspan=1),
                 pl.subplot2grid((3, 5), (2, 4), colspan=1, rowspan=1)]

        axt = pl.axes([0.15, 0.9, 0.75, 0.075])
        axt.axis('off')
        axt.annotate("%s %d" % (self.IDSTRING, self.ID),
                     xy=(0.5, 0.5), xycoords='axes fraction',
                     ha='center', va='center', fontsize=18)
        axt.annotate(r"Unmasked overfitting metric",
                     xy=(0.5, 0.2), xycoords='axes fraction',
                     ha='center', va='center', fontsize=12, color='k')

        # Loop over three depths
        for depth, ax, axh in zip([0.01, 0.001, 0.0001], axes, axesh):

            # The unmasked overfitting metric
            metric = 1 - (self._O2 + self._O3 / depth) / self._O1
            mm = np.delete(metric, self.mask)
            med = np.nanmedian(mm)
            mad = np.nanmedian(np.abs(mm - med))

            # Plot the metric as a function of time
            ax.plot(np.delete(self.t, self.mask), np.delete(
                metric, self.mask), 'k.', alpha=0.5, ms=2)
            ax.plot(np.delete(self.t, self.mask), np.delete(
                metric, self.mask), 'k-', alpha=0.1, lw=0.5)
            ax.plot(self.t[self.mask], metric[self.mask],
                    'r.', alpha=0.25, ms=2)
            ylim = (-0.2, 1.0)
            ax.margins(0, None)
            ax.axhline(0, color='k', lw=1, alpha=0.5)
            ax.set_ylim(*ylim)
            ax.set_xlabel('Time (days)', fontsize=14)
            ax.set_ylabel('Overfitting (d = %s)' % str(depth), fontsize=14)

            # Plot the histogram
            m1 = np.delete(metric, self.mask)
            m2 = metric
            axh.hist(m1, bins=30, range=ylim, orientation="horizontal",
                     histtype="step", fill=False, color='k')
            axh.hist(m2, bins=30, range=ylim, orientation="horizontal",
                     histtype="step", fill=False, color='r', alpha=0.3)
            axh.axhline(med, color='b', ls='-', lw=1)
            axh.axhspan(med - mad, med + mad, color='b', alpha=0.1)
            axh.axhline(0, color='k', lw=1, alpha=0.5)
            axh.yaxis.tick_right()
            axh.set_ylim(*ax.get_ylim())
            axh.set_xticklabels([])

            bbox = dict(fc="w", ec="1", alpha=0.5)
            info = r"$\mathrm{med}=%.3f$" % med + \
                "\n" + r"$\mathrm{mad}=%.3f$" % mad
            axh.annotate(info, xy=(0.1, 0.925), xycoords='axes fraction',
                         ha="left", va="top", bbox=bbox, color="b")

        return fig

    def plot_masked(self, depth=0.001):
        '''

        '''

        # Overfitting metric
        metric = self._O5 / depth

        # Set up
        fig = pl.figure(figsize=(10, 4))
        fig.subplots_adjust(bottom=0.2)
        ax = pl.subplot2grid((1, 5), (0, 0), colspan=4, rowspan=1)
        axh = pl.subplot2grid((1, 5), (0, 4), colspan=1, rowspan=1)
        ax.set_xlabel('Time (days)', fontsize=14)
        ax.set_ylabel('Masked overfitting', fontsize=14)

        # Plot the metric as a function of time
        ax.plot(np.delete(self.time, self.mask), np.delete(
            metric, self.mask), 'k.', alpha=0.5, ms=2)
        ax.plot(np.delete(self.time, self.mask), np.delete(
            metric, self.mask), 'k-', alpha=0.1, lw=0.5)
        ax.plot(self.time[self.mask], metric[self.mask],
                'r.', alpha=0.25, ms=2)

        # Bound 99% of data
        y = np.delete(metric, self.mask)
        N = int(0.99 * len(y))
        hi, lo = y[np.argsort(y)][[N, -N]]
        fsort = y[np.argsort(y)]
        pad = (hi - lo) * 0.1
        ylim = (lo - pad, hi + pad)
        ax.set_ylim(ylim)

        # Plot the histogram
        m1 = np.delete(metric, self.mask)
        m2 = metric
        axh.hist(m1, bins=30, range=ylim, orientation="horizontal",
                 histtype="step", fill=False, color='k')
        axh.hist(m2, bins=30, range=ylim, orientation="horizontal",
                 histtype="step", fill=False, color='r', alpha=0.3)
        axh.yaxis.tick_right()
        axh.set_ylim(*ax.get_ylim())
        axh.set_xticklabels([])

        # Return
        return fig, ax, axh

    def plot_joint(self, depth=0.001):
        '''

        '''

        # TODO
        pass

    def plot_marginal(self, depth=0.001, std=0.001, mu=0.001):
        '''

        '''

        # TODO
        pass

    def _compute_unmasked(self):
        '''
        Computes the unmasked overfitting metric.

        '''

        # Some masked arrays/matrices
        mt = np.delete(self.t, self.mask)
        my = np.delete(self.y, self.mask)
        MC = np.delete(np.delete(self.C, self.mask, axis=0), self.mask, axis=1)

        # Compute the full `model` and the de-trended flux `f`
        model = np.dot(self.XLX, np.linalg.solve(MC, my))
        model -= np.nanmedian(model)
        f = self.y - model

        # The regression matrix
        log.info("Unmasked light curve: pre-computing some stuff...")
        R = np.linalg.solve(MC, self.XLX.T).T

        # Arrays we'll use for the overfitting calculations
        self._O1 = np.zeros_like(self.t)
        self._O2 = np.zeros_like(self.t)
        self._O3 = np.zeros_like(self.t)
        self._O4 = np.zeros_like(self.t)

        # Loop over all cadences
        log.info("Unmasked light curve: looping through the dataset...")
        for k in prange(len(self.t)):

            # Evaluate the transit model...
            TAU = self.tau(self.t, t0=self.t[k])
            i = np.where(TAU < 0)[0]
            TAU = TAU.reshape(-1, 1)

            # ...and the masked transit model
            MTAU = self.tau(mt, t0=self.t[k])
            mi = np.where(MTAU < 0)[0]
            MTAU = MTAU.reshape(-1, 1)

            # Sparse algebra
            A = np.dot(np.dot(TAU[i].T, self.Kinv[i, :][:, i]), TAU[i])
            B = np.dot(TAU[i].T, self.Kinv[i, :])
            C = TAU - np.dot(R[:, mi], MTAU[mi])

            # Compute the overfitting metric
            self._O1[k] = A
            self._O2[k] = np.dot(B, C)
            self._O3[k] = np.dot(B, f)
            self._O4[k] = np.dot(B, self.y)

    def _compute_masked(self):
        '''
        Computes the masked overfitting metric.

        '''

        # Remove the masked elements for simplicity
        t = np.delete(self.t, self.mask)
        y = np.delete(self.y, self.mask)
        C = np.delete(np.delete(self.C, self.mask, axis=0), self.mask, axis=1)
        X = np.delete(self.X, self.mask, axis=0)
        XL = np.delete(self.XL, self.mask, axis=0)
        K = np.delete(np.delete(self.K, self.mask, axis=0), self.mask, axis=1)
        Kinv = np.linalg.inv(K)

        # Initialize
        N = len(y)
        self._O5 = np.zeros(N) * np.nan

        # Get the solution to the masked problem

        # DEBUG!!!
        import os
        if os.path.exists('A.npz'):
            A = np.load('A.npz')['A']
        else:
            A = MaskSolve(C, y, w=self.w)
            np.savez('A.npz', A=A)

        # The initial condition
        mask = np.arange(0, self.w)
        P = np.dot(XL, np.delete(X, mask, axis=0).T)

        # Loop over all starting points for the mask
        for n in prange(N - self.w + 1):

            # The mask for this iteration goes from `n` to `n + w - 1`
            mask = np.arange(n, self.w + n)

            # The center point of the mask is therefore
            j = n + (self.w + 1) // 2 - 1

            # Update the P matrix
            # DEBUG P[:,?] = np.dot(self.X[?,:], self.XL.T)
            P = np.dot(XL, np.delete(X, mask, axis=0).T)
            
            # The linear solution to this step
            m = np.dot(P, A[n])

            # Evaluate the sparse transit model
            TAU = self.tau(t, t0=t[j])
            i = np.where(TAU < 0)[0]
            TAU = TAU[i].reshape(-1, 1)

            # Dot the transit model in
            den = np.dot(np.dot(TAU.T, Kinv[i, :][:, i]), TAU)
            num = np.dot(TAU.T, Kinv[i, :])

            # Compute the overfitting metric
            # Divide this number by a depth
            # to get the overfitting for that
            # particular depth.
            self._O5[j] = -np.dot(num, y - m) / den


class SavedOverfitObject(Overfit):
    '''

    '''

    def __init__(self, ID, IDSTRING, t, mask, O1, O2, O3, O4, O5):
        '''

        '''

        self.ID = ID
        self.IDSTRING = IDSTRING
        self._O1 = O1
        self._O2 = O2
        self._O3 = O3
        self._O4 = O4
        self._O5 = O5
        self.t = t
        self.mask = mask
