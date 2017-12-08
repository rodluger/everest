#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`inject.py` - Transit injection
---------------------------------------

The :py:class:`Inject` class is the model that handles transit injection
and recovery.

'''

from __future__ import division, print_function, absolute_import, \
     unicode_literals
from .detrender import *
from .transit import Transit
from .dvs import DVS
import os
import sys
import numpy as np
import matplotlib.pyplot as pl
import traceback
import logging
log = logging.getLogger(__name__)

__all__ = ['Inject']


def Inject(ID, inj_model='nPLD', t0=None, per=None, dur=0.1, depth=0.001,
           mask=False, trn_win=5, poly_order=3, make_fits=False, **kwargs):
    '''
    Run one of the :py:obj:`everest` models with injected transits and attempt
    to recover the transit depth at the end with a simple linear regression
    with a polynomial baseline. The depth is stored in the
    :py:obj:`inject` attribute of the model (a dictionary) as
    :py:obj:`rec_depth`. A control injection is also performed, in which the
    transits are injected into the de-trended data; the recovered depth in
    the control run is stored in :py:obj:`inject`
    as :py:obj:`rec_depth_control`.

    :param int ID: The target id
    :param str inj_model: The name of the :py:obj:`everest` model to run. \
           Default `"nPLD"`
    :param float t0: The transit ephemeris in days. Default is to draw from \
           the uniform distributon [0., :py:obj:`per`)
    :param float per: The injected planet period in days. Default is to draw \
           from the uniform distribution [2, 10]
    :param float dur: The transit duration in days. Must be in the range \
           [0.05, 0.5]. Default 0.1
    :param float depth: The fractional transit depth. Default 0.001
    :param bool mask: Explicitly mask the in-transit cadences when computing \
           the PLD model? Default :py:obj:`False`
    :param float trn_win: The size of the transit window in units of the \
           transit duration
    :param int poly_order: The order of the polynomial used to fit the \
           continuum

    '''

    # Randomize the planet params
    if per is None:
        a = 3.
        b = 10.
        per = a + (b - a) * np.random.random()
    if t0 is None:
        t0 = per * np.random.random()

    # Get the actual class
    _model = eval(inj_model)
    inject = {'t0': t0, 'per': per, 'dur': dur, 'depth': depth, 'mask': mask,
              'poly_order': poly_order, 'trn_win': trn_win}

    # Define the injection class
    class Injection(_model):
        '''
        The :py:obj:`Injection` class is a special subclass of a
        user-selected :py:obj:`everest` model.
        See :py:func:`Inject` for more details.

        '''

        def __init__(self, *args, **kwargs):
            '''

            '''

            self.inject = kwargs.pop('inject', None)
            self.parent_class = kwargs.pop('parent_class', None)
            self.kwargs = kwargs
            super(Injection, self).__init__(*args, **kwargs)

        @property
        def name(self):
            '''

            '''

            if self.inject['mask']:
                maskchar = 'M'
            else:
                maskchar = 'U'
            return '%s_Inject_%s%g' % (self.parent_class,
                                       maskchar, self.inject['depth'])

        def load_tpf(self):
            '''
            Loads the target pixel files and injects transits at the pixel level.

            '''

            # Load the TPF
            super(Injection, self).load_tpf()
            log.info("Injecting transits...")

            # Inject the transits into the regular data
            transit_model = Transit(
                self.time, t0=self.inject['t0'], per=self.inject['per'],
                dur=self.inject['dur'], depth=self.inject['depth'])
            for i in range(self.fpix.shape[1]):
                self.fpix[:, i] *= transit_model
            self.fraw = np.sum(self.fpix, axis=1)
            if self.inject['mask']:
                self.transitmask = np.array(list(set(np.concatenate(
                    [self.transitmask, np.where(transit_model < 1.)[0]]))),
                    dtype=int)

            # Update the PLD normalization
            self.get_norm()

        def recover_depth(self):
            '''
            Recovers the injected transit depth from the long
            cadence data with a simple LLS solver.
            The results are all stored in the :py:obj:`inject`
            attribute of the model.

            '''

            # Control run
            transit_model = Transit(
                self.time, t0=self.inject['t0'], per=self.inject['per'],
                dur=self.inject['dur'], depth=self.inject['depth'])
            kwargs = dict(self.kwargs)
            kwargs.update({'clobber': False})
            control = eval(self.parent_class)(
                self.ID, is_parent=True, **kwargs)
            control.fraw *= transit_model

            # Get params
            log.info("Recovering transit depth...")
            t0 = self.inject['t0']
            per = self.inject['per']
            dur = self.inject['dur']
            depth = self.inject['depth']
            trn_win = self.inject['trn_win']
            poly_order = self.inject['poly_order']

            for run, tag in zip([self, control], ['', '_control']):

                # Compute the model
                mask = np.array(
                    list(set(np.concatenate([run.badmask, run.nanmask]))),
                    dtype=int)
                flux = np.delete(run.flux / np.nanmedian(run.flux), mask)
                time = np.delete(run.time, mask)
                transit_model = (Transit(time, t0=t0, per=per,
                                         dur=dur, depth=depth) - 1) / depth

                # Count the transits
                t0 += np.ceil((time[0] - dur - t0) / per) * per
                ttimes0 = np.arange(t0, time[-1] + dur, per)
                tinds = []
                for tt in ttimes0:
                    # Get indices for this chunk
                    inds = np.where(np.abs(time - tt) < trn_win * dur / 2.)[0]
                    # Ensure there's a transit in this chunk, and that
                    # there are enough points for the polynomial fit
                    if np.any(transit_model[inds] < 0.) and \
                            len(inds) > poly_order:
                        tinds.append(inds)

                # Our design matrix
                sz = (poly_order + 1) * len(tinds)
                X = np.empty((0, 1 + sz), dtype=float)
                Y = np.array([], dtype=float)
                T = np.array([], dtype=float)

                # Loop over all transits
                for i, inds in enumerate(tinds):
                    # Get the transit model
                    trnvec = transit_model[inds].reshape(-1, 1)
                    # Normalize the time array
                    t = time[inds]
                    t = (t - t[0]) / (t[-1] - t[0])
                    # Cumulative arrays
                    T = np.append(T, time[inds])
                    Y = np.append(Y, flux[inds])
                    # Polynomial vector
                    polyvec = np.array(
                        [t ** o for o in range(0, poly_order + 1)]).T
                    # Update the design matrix with this chunk
                    lzeros = np.zeros((len(t), i * (poly_order + 1)))
                    rzeros = np.zeros(
                        (len(t), sz - (i + 1) * (poly_order + 1)))
                    chunk = np.hstack((trnvec, lzeros, polyvec, rzeros))
                    X = np.vstack((X, chunk))

                # Get the relative depth
                A = np.dot(X.T, X)
                B = np.dot(X.T, Y)
                C = np.linalg.solve(A, B)
                rec_depth = C[0]

                # Get the uncertainties
                sig = 1.4826 * \
                    np.nanmedian(np.abs(flux - np.nanmedian(flux))
                                 ) / np.nanmedian(flux)
                cov = sig ** 2 * np.linalg.solve(A, np.eye(A.shape[0]))
                err = np.sqrt(np.diag(cov))
                rec_depth_err = err[0]

                # Store the results
                self.inject.update(
                    {'rec_depth%s' % tag: rec_depth,
                     'rec_depth_err%s' % tag: rec_depth_err})

                # Store the detrended, folded data
                D = (Y - np.dot(C[1:], X[:, 1:].T) +
                     np.nanmedian(Y)) / np.nanmedian(Y)
                T = (T - t0 - per / 2.) % per - per / 2.
                self.inject.update(
                    {'fold_time%s' % tag: T, 'fold_flux%s' % tag: D})

        def plot_final(self, ax):
            '''
            Plots the injection recovery results.

            '''

            from mpl_toolkits.axes_grid.inset_locator import inset_axes
            ax.axis('off')
            ax1 = inset_axes(ax, width="47%", height="100%", loc=6)
            ax2 = inset_axes(ax, width="47%", height="100%", loc=7)

            # Plot the recovered folded transits
            ax1.plot(self.inject['fold_time'],
                     self.inject['fold_flux'], 'k.', alpha=0.3)
            x = np.linspace(np.min(self.inject['fold_time']), np.max(
                self.inject['fold_time']), 500)
            try:
                y = Transit(
                    x, t0=0., per=self.inject['per'], dur=self.inject['dur'],
                    depth=self.inject['rec_depth'])
            except:
                # Log the error, and carry on
                exctype, value, tb = sys.exc_info()
                for line in traceback.format_exception_only(exctype, value):
                    l = line.replace('\n', '')
                    log.error(l)
                y = np.ones_like(x) * np.nan
            ax1.plot(x, y, 'r-')
            ax1.annotate('INJECTED', xy=(0.98, 0.025),
                         xycoords='axes fraction',
                         ha='right', va='bottom', fontsize=10, alpha=0.5,
                         fontweight='bold')
            ax1.annotate('True depth:\nRecovered depth:',
                         xy=(0.02, 0.025),
                         xycoords='axes fraction',
                         ha='left', va='bottom', fontsize=6, color='r')
            ax1.annotate('%.6f\n%.6f' % (self.inject['depth'],
                         self.inject['rec_depth']),
                         xy=(0.4, 0.025),
                         xycoords='axes fraction',
                         ha='left', va='bottom', fontsize=6, color='r')
            ax1.margins(0, None)
            ax1.ticklabel_format(useOffset=False)

            # Plot the recovered folded transits (control)
            ax2.plot(self.inject['fold_time_control'],
                     self.inject['fold_flux_control'], 'k.', alpha=0.3)
            x = np.linspace(np.min(self.inject['fold_time_control']), np.max(
                self.inject['fold_time_control']), 500)
            try:
                y = Transit(
                    x, t0=0., per=self.inject['per'], dur=self.inject['dur'],
                    depth=self.inject['rec_depth_control'])
            except:
                # Log the error, and carry on
                exctype, value, tb = sys.exc_info()
                for line in traceback.format_exception_only(exctype, value):
                    l = line.replace('\n', '')
                    log.error(l)
                y = np.ones_like(x) * np.nan
            ax2.plot(x, y, 'r-')
            ax2.annotate('CONTROL', xy=(0.98, 0.025), xycoords='axes fraction',
                         ha='right', va='bottom', fontsize=10, alpha=0.5,
                         fontweight='bold')
            ax2.annotate('True depth:\nRecovered depth:',
                         xy=(0.02, 0.025),
                         xycoords='axes fraction',
                         ha='left', va='bottom', fontsize=6, color='r')
            ax2.annotate('%.6f\n%.6f' % (self.inject['depth'],
                         self.inject['rec_depth_control']),
                         xy=(0.4, 0.025),
                         xycoords='axes fraction',
                         ha='left', va='bottom', fontsize=6, color='r')
            ax2.margins(0, None)
            ax2.ticklabel_format(useOffset=False)
            N = int(0.995 * len(self.inject['fold_flux_control']))
            hi, lo = self.inject['fold_flux_control'][np.argsort(
                self.inject['fold_flux_control'])][[N, -N]]
            fsort = self.inject['fold_flux_control'][np.argsort(
                self.inject['fold_flux_control'])]
            pad = (hi - lo) * 0.2
            ylim = (lo - 2 * pad, hi + pad)
            ax2.set_ylim(ylim)
            ax1.set_ylim(ylim)

            ax2.set_yticklabels([])
            for tick in ax1.get_xticklabels() + ax1.get_yticklabels() + \
                    ax2.get_xticklabels():
                tick.set_fontsize(5)

        def finalize(self):
            '''
            Calls the depth recovery routine at the end
            of the de-trending step.

            '''

            super(Injection, self).finalize()
            self.recover_depth()

    return Injection(ID, inject=inject, parent_class=inj_model,
                     make_fits=make_fits, **kwargs)
