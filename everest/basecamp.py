#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`basecamp.py` - The Everest base class
----------------------------------------------

The :py:obj:`everest` engine. All :py:obj:`everest` models
inherit from :py:class:`Basecamp`.

'''

from __future__ import division, print_function, absolute_import, \
    unicode_literals
from . import missions
from .utils import AP_SATURATED_PIXEL, prange
from .mathutils import SavGol
from .masksolve import MaskSolve
from .gp import GetCovariance
from .search import Search
from .transit import TransitModel, TransitShape
from .dvs import OVERFIT
from scipy.linalg import block_diag, cholesky, cho_factor, cho_solve
import os
import numpy as np
import matplotlib.pyplot as pl
from scipy.ndimage import zoom
from itertools import combinations_with_replacement as multichoose
import logging
import platform
import subprocess
log = logging.getLogger(__name__)

__all__ = ['Basecamp', 'Overfitting']


class Overfitting(object):
    """Stores information on the overfitting metrics for a light curve."""

    def __init__(self, O1, O2, O3, O4, O5, pdf):
        """Store values."""
        self._O1 = O1
        self._O2 = O2
        self._O3 = O3
        self._O4 = O4
        self._O5 = O5
        self.pdf = pdf

    def masked(self, depth=0.01):
        """Return the masked overfitting metric for a given transit depth."""
        return np.hstack(self._O5) / depth

    def unmasked(self, depth=0.01):
        """Return the unmasked overfitting metric for a given transit depth."""
        return 1 - (np.hstack(self._O2) +
                    np.hstack(self._O3) / depth) / np.hstack(self._O1)

    def show(self):
        """Show the overfitting PDF summary."""
        try:
            if platform.system().lower().startswith('darwin'):
                subprocess.call(['open', self.pdf])
            elif os.name == 'nt':
                os.startfile(self.pdf)
            elif os.name == 'posix':
                subprocess.call(['xdg-open', self.pdf])
            else:
                raise IOError("")
        except IOError:
            log.info("Unable to open the pdf. Try opening it manually:")
            log.info(self.pdf)

class Basecamp(object):
    '''

    '''

    @property
    def _mission(self):
        '''

        '''

        return getattr(missions, self.mission)

    @_mission.setter
    def _mission(self, value):
        '''

        '''

        raise NotImplementedError("Can't set this property.")

    @property
    def dir(self):
        '''
        Returns the directory where the raw data and output for the target is
        stored.

        '''

        return self._mission.TargetDirectory(self.ID, self.season)

    @dir.setter
    def dir(self, value):
        '''

        '''

        raise NotImplementedError("Can't set this property.")

    @property
    def logfile(self):
        '''
        Returns the full path to the log file for the current run.

        '''

        return os.path.join(self.dir, '%s.log' % self.name)

    @logfile.setter
    def logfile(self, value):
        '''

        '''

        raise NotImplementedError("Can't set this property.")

    @property
    def season(self):
        """
        Return the current observing season.

        For *K2*, this is the observing campaign, while for *Kepler*,
        it is the current quarter.

        """
        try:
            self._season
        except AttributeError:
            self._season = self._mission.Season(self.ID)
            if hasattr(self._season, '__len__'):
                raise AttributeError(
                    "Please choose a campaign/season for this target: %s." %
                    self._season)
        return self._season

    @season.setter
    def season(self, value):
        '''

        '''

        raise NotImplementedError("Can't set this property.")

    @property
    def flux(self):
        '''
        The corrected/de-trended flux. This is computed by subtracting
        the linear model from the raw SAP flux.

        '''

        return self.fraw - self.model

    @flux.setter
    def flux(self, value):
        '''

        '''

        raise NotImplementedError("Can't set this property.")

    @property
    def fcor(self):
        '''
        The CBV-corrected de-trended flux.

        '''

        if self.XCBV is None:
            return None
        else:
            return self.flux - self._mission.FitCBVs(self)

    @fcor.setter
    def fcor(self, value):
        '''

        '''

        raise NotImplementedError("Can't set this property.")

    @property
    def norm(self):
        '''
        The PLD normalization. Typically, this is just the simple aperture
        photometry flux (i.e., the sum of all the pixels in the aperture).

        '''

        return self._norm

    @norm.setter
    def norm(self, value):
        '''

        '''

        raise NotImplementedError("Can't set this property.")

    @property
    def cdpps(self):
        '''
        The string version of the current value of the CDPP in *ppm*. This
        displays the CDPP for each segment of the light curve individually
        (if breakpoints are present).

        '''

        return " / ".join(["%.2f ppm" % c for c in self.cdpp_arr]) + \
               (" (%.2f ppm)" % self.cdpp)

    @cdpps.setter
    def cdpps(self, value):
        '''

        '''

        raise NotImplementedError("Can't set this property.")

    @property
    def mask(self):
        '''
        The array of indices to be masked. This is the union of the sets of
        outliers, bad (flagged) cadences, transit cadences, and :py:obj:`NaN`
        cadences.

        '''

        return np.array(list(set(np.concatenate([self.outmask, self.badmask,
                        self.transitmask, self.nanmask]))), dtype=int)

    @mask.setter
    def mask(self, value):
        '''

        '''

        raise NotImplementedError("Can't set this property.")

    @property
    def weights(self):
        '''
        The PLD weights vector. The model may be computed by dotting the design
        matrix :py:attr:`X` with this vector. Note that these are computed just
        for plotting purpoeses -- the actual weights are never explicitly
        computed during the de-trending, since it can be rather slow.

        '''

        if self._weights is None:
            self.get_weights()
        return self._weights

    @weights.setter
    def weights(self, value):
        '''

        '''

        raise NotImplementedError("Can't set this property.")

    @property
    def transit_model(self):
        '''

        '''

        try:
            self._transit_model
        except AttributeError:
            self._transit_model = None
        return self._transit_model

    @transit_model.setter
    def transit_model(self, val):
        '''

        '''

        if val is None:
            self._transit_model = None
            self.transit_depth = None
        else:
            val = np.atleast_1d(val)
            for tm in val:
                assert type(tm) is TransitModel, \
                       "Kwarg `transit_model` must be an instance or " + \
                       "a list of instances of `everest.TransitModel`."
            self._transit_model = val
            self.transit_depth = None

    def get_norm(self):
        '''
        Computes the PLD normalization. In the base class, this is just
        the sum of all the pixel fluxes.

        '''

        self._norm = self.fraw

    def X(self, i, j=slice(None, None, None)):
        '''
        Computes the design matrix at the given *PLD* order and the given
        indices. The columns are the *PLD* vectors for the target at the
        corresponding order, computed as the product of the fractional pixel
        flux of all sets of :py:obj:`n` pixels, where :py:obj:`n` is the *PLD*
        order.

        '''

        X1 = self.fpix[j] / self.norm[j].reshape(-1, 1)
        X = np.product(list(multichoose(X1.T, i + 1)), axis=1).T
        if self.X1N is not None:
            return np.hstack([X, self.X1N[j] ** (i + 1)])
        else:
            return X

    def plot_info(self, dvs):
        '''
        Plots miscellaneous de-trending information on the data
        validation summary figure.

        :param dvs: A :py:class:`dvs.DVS` figure instance

        '''

        axl, axc, axr = dvs.title()
        axc.annotate("%s %d" % (self._mission.IDSTRING, self.ID),
                     xy=(0.5, 0.5), xycoords='axes fraction',
                     ha='center', va='center', fontsize=18)

        axc.annotate(r"%.2f ppm $\rightarrow$ %.2f ppm" %
                     (self.cdppr, self.cdpp),
                     xy=(0.5, 0.2), xycoords='axes fraction',
                     ha='center', va='center', fontsize=8, color='k',
                     fontstyle='italic')

        axl.annotate("%s %s%02d: %s" %
                     (self.mission.upper(),
                      self._mission.SEASONCHAR, self.season, self.name),
                     xy=(0.5, 0.5), xycoords='axes fraction',
                     ha='center', va='center', fontsize=12,
                     color='k')

        axl.annotate(self.aperture_name if len(self.neighbors) == 0
                     else "%s, %d neighbors" %
                     (self.aperture_name, len(self.neighbors)),
                     xy=(0.5, 0.2), xycoords='axes fraction',
                     ha='center', va='center', fontsize=8, color='k',
                     fontstyle='italic')

        axr.annotate("%s %.3f" % (self._mission.MAGSTRING, self.mag),
                     xy=(0.5, 0.5), xycoords='axes fraction',
                     ha='center', va='center', fontsize=12,
                     color='k')

        if not np.isnan(self.cdppg) and self.cdppg > 0:
            axr.annotate(r"GP %.3f ppm" % (self.cdppg),
                         xy=(0.5, 0.2), xycoords='axes fraction',
                         ha='center', va='center', fontsize=8, color='k',
                         fontstyle='italic')

    def compute(self):
        '''
        Compute the model for the current value of lambda.

        '''

        # Is there a transit model?
        if self.transit_model is not None:
            return self.compute_joint()

        log.info('Computing the model...')

        # Loop over all chunks
        model = [None for b in self.breakpoints]
        for b, brkpt in enumerate(self.breakpoints):

            # Masks for current chunk
            m = self.get_masked_chunk(b)
            c = self.get_chunk(b)

            # This block of the masked covariance matrix
            mK = GetCovariance(self.kernel, self.kernel_params,
                               self.time[m], self.fraw_err[m])

            # Get median
            med = np.nanmedian(self.fraw[m])

            # Normalize the flux
            f = self.fraw[m] - med

            # The X^2 matrices
            A = np.zeros((len(m), len(m)))
            B = np.zeros((len(c), len(m)))

            # Loop over all orders
            for n in range(self.pld_order):

                # Only compute up to the current PLD order
                if (self.lam_idx >= n) and (self.lam[b][n] is not None):
                    XM = self.X(n, m)
                    XC = self.X(n, c)
                    A += self.lam[b][n] * np.dot(XM, XM.T)
                    B += self.lam[b][n] * np.dot(XC, XM.T)
                    del XM, XC

            # Compute the model
            W = np.linalg.solve(mK + A, f)
            model[b] = np.dot(B, W)

        # Free up some memory
        del A, B, W

        # Join the chunks after applying the correct offset
        if len(model) > 1:

            # First chunk
            self.model = model[0][:-self.bpad]

            # Center chunks
            for m in model[1:-1]:
                # Join the chunks at the first non-outlier cadence
                i = 1
                while len(self.model) - i in self.mask:
                    i += 1
                offset = self.model[-i] - m[self.bpad - i]
                self.model = np.concatenate(
                    [self.model, m[self.bpad:-self.bpad] + offset])

            # Last chunk
            i = 1
            while len(self.model) - i in self.mask:
                i += 1
            offset = self.model[-i] - model[-1][self.bpad - i]
            self.model = np.concatenate(
                [self.model, model[-1][self.bpad:] + offset])

        else:

            self.model = model[0]

        # Subtract the global median
        self.model -= np.nanmedian(self.model)

        # Get the CDPP and reset the weights
        self.cdpp_arr = self.get_cdpp_arr()
        self.cdpp = self.get_cdpp()
        self._weights = None

    def compute_joint(self):
        '''
        Compute the model in a single step, allowing for a light curve-wide
        transit model. This is a bit more expensive to compute.

        '''

        # Init
        log.info('Computing the joint model...')
        A = [None for b in self.breakpoints]
        B = [None for b in self.breakpoints]

        # We need to make sure that we're not masking the transits we are
        # trying to fit!
        # NOTE: If there happens to be an index that *SHOULD* be masked during
        # a transit (cosmic ray, detector anomaly), update `self.badmask`
        # to include that index.
        # Bad data points are *never* used in the regression.
        if self.transit_model is not None:
            outmask = np.array(self.outmask)
            transitmask = np.array(self.transitmask)
            transit_inds = np.where(
                np.sum([tm(self.time) for tm in self.transit_model],
                       axis=0) < 0)[0]
            self.outmask = np.array(
                [i for i in self.outmask if i not in transit_inds])
            self.transitmask = np.array(
                [i for i in self.transitmask if i not in transit_inds])

        # Loop over all chunks
        for b, brkpt in enumerate(self.breakpoints):

            # Masks for current chunk
            m = self.get_masked_chunk(b, pad=False)
            c = self.get_chunk(b, pad=False)

            # The X^2 matrices
            A[b] = np.zeros((len(m), len(m)))
            B[b] = np.zeros((len(c), len(m)))

            # Loop over all orders
            for n in range(self.pld_order):

                # Only compute up to the current PLD order
                if (self.lam_idx >= n) and (self.lam[b][n] is not None):
                    XM = self.X(n, m)
                    XC = self.X(n, c)
                    A[b] += self.lam[b][n] * np.dot(XM, XM.T)
                    B[b] += self.lam[b][n] * np.dot(XC, XM.T)
                    del XM, XC

        # Merge chunks. BIGA and BIGB are sparse, but unfortunately
        # scipy.sparse doesn't handle sparse matrix inversion all that
        # well when the *result* is not itself sparse. So we're sticking
        # with regular np.linalg.
        BIGA = block_diag(*A)
        del A
        BIGB = block_diag(*B)
        del B

        # Compute the full covariance matrix
        mK = GetCovariance(self.kernel, self.kernel_params, self.apply_mask(
             self.time), self.apply_mask(self.fraw_err))

        # The normalized, masked flux array
        f = self.apply_mask(self.fraw)
        med = np.nanmedian(f)
        f -= med

        # Are we computing a joint transit model?
        if self.transit_model is not None:

            # Get the unmasked indices
            m = self.apply_mask()

            # Subtract off the mean total transit model
            mean_transit_model = med * \
                np.sum([tm.depth * tm(self.time[m])
                        for tm in self.transit_model], axis=0)
            f -= mean_transit_model

            # Now add each transit model to the matrix of regressors
            for tm in self.transit_model:
                XM = tm(self.time[m]).reshape(-1, 1)
                XC = tm(self.time).reshape(-1, 1)
                BIGA += med ** 2 * tm.var_depth * np.dot(XM, XM.T)
                BIGB += med ** 2 * tm.var_depth * np.dot(XC, XM.T)
                del XM, XC

            # Dot the inverse of the covariance matrix
            W = np.linalg.solve(mK + BIGA, f)
            self.model = np.dot(BIGB, W)

            # Compute the transit weights and maximum likelihood transit model
            w_trn = med ** 2 * np.concatenate([tm.var_depth * np.dot(
                    tm(self.time[m]).reshape(1, -1), W)
                    for tm in self.transit_model])
            self.transit_depth = np.array(
                 [med * tm.depth + w_trn[i] for i, tm in
                  enumerate(self.transit_model)]) / med

            # Remove the transit prediction from the model
            self.model -= np.dot(np.hstack([tm(self.time).reshape(-1, 1)
                                 for tm in self.transit_model]),
                                 w_trn)

        else:

            # No transit model to worry about
            W = np.linalg.solve(mK + BIGA, f)
            self.model = np.dot(BIGB, W)

        # Subtract the global median
        self.model -= np.nanmedian(self.model)

        # Restore the mask
        if self.transit_model is not None:
            self.outmask = outmask
            self.transitmask = transitmask

        # Get the CDPP and reset the weights
        self.cdpp_arr = self.get_cdpp_arr()
        self.cdpp = self.get_cdpp()
        self._weights = None

    def apply_mask(self, x=None):
        '''
        Returns the outlier mask, an array of indices corresponding to the
        non-outliers.

        :param numpy.ndarray x: If specified, returns the masked version of \
               :py:obj:`x` instead. Default :py:obj:`None`

        '''

        if x is None:
            return np.delete(np.arange(len(self.time)), self.mask)
        else:
            return np.delete(x, self.mask, axis=0)

    def get_chunk(self, b, x=None, pad=True):
        '''
        Returns the indices corresponding to a given light curve chunk.

        :param int b: The index of the chunk to return
        :param numpy.ndarray x: If specified, applies the mask to array \
               :py:obj:`x`. Default :py:obj:`None`

        '''

        M = np.arange(len(self.time))
        if b > 0:
            res = M[(M > self.breakpoints[b - 1] - int(pad) * self.bpad)
                    & (M <= self.breakpoints[b] + int(pad) * self.bpad)]
        else:
            res = M[M <= self.breakpoints[b] + int(pad) * self.bpad]
        if x is None:
            return res
        else:
            return x[res]

    def get_masked_chunk(self, b, x=None, pad=True):
        '''
        Same as :py:meth:`get_chunk`, but first removes the outlier indices.
        :param int b: The index of the chunk to return
        :param numpy.ndarray x: If specified, applies the mask to \
               array :py:obj:`x`. Default :py:obj:`None`

        '''

        M = self.apply_mask(np.arange(len(self.time)))
        if b > 0:
            res = M[(M > self.breakpoints[b - 1] - int(pad) * self.bpad)
                    & (M <= self.breakpoints[b] + int(pad) * self.bpad)]
        else:
            res = M[M <= self.breakpoints[b] + int(pad) * self.bpad]
        if x is None:
            return res
        else:
            return x[res]

    def get_weights(self):
        '''
        Computes the PLD weights vector :py:obj:`w`.

        ..warning :: Deprecated and not thoroughly tested.

        '''

        log.info("Computing PLD weights...")

        # Loop over all chunks
        weights = [None for i in range(len(self.breakpoints))]
        for b, brkpt in enumerate(self.breakpoints):

            # Masks for current chunk
            m = self.get_masked_chunk(b)
            c = self.get_chunk(b)

            # This block of the masked covariance matrix
            _mK = GetCovariance(self.kernel, self.kernel_params,
                                self.time[m], self.fraw_err[m])

            # This chunk of the normalized flux
            f = self.fraw[m] - np.nanmedian(self.fraw)

            # Loop over all orders
            _A = [None for i in range(self.pld_order)]
            for n in range(self.pld_order):
                if self.lam_idx >= n:
                    X = self.X(n, m)
                    _A[n] = np.dot(X, X.T)
                    del X

            # Compute the weights
            A = np.sum([l * a for l, a in zip(self.lam[b], _A)
                        if l is not None], axis=0)
            W = np.linalg.solve(_mK + A, f)
            weights[b] = [l * np.dot(self.X(n, m).T, W)
                          for n, l in enumerate(self.lam[b]) if l is not None]

        self._weights = weights

    def get_cdpp_arr(self, flux=None):
        '''
        Returns the CDPP value in *ppm* for each of the
        chunks in the light curve.

        '''

        if flux is None:
            flux = self.flux
        return np.array([self._mission.CDPP(flux[self.get_masked_chunk(b)],
                        cadence=self.cadence)
                        for b, _ in enumerate(self.breakpoints)])

    def get_cdpp(self, flux=None):
        '''
        Returns the scalar CDPP for the light curve.

        '''

        if flux is None:
            flux = self.flux
        return self._mission.CDPP(self.apply_mask(flux), cadence=self.cadence)

    def plot_aperture(self, axes, labelsize=8):
        '''
        Plots the aperture and the pixel images at the beginning, middle,
        and end of the time series. Also plots a high resolution image of
        the target, if available.

        '''

        log.info('Plotting the aperture...')

        # Get colormap
        plasma = pl.get_cmap('plasma')
        plasma.set_bad(alpha=0)

        # Get aperture contour
        def PadWithZeros(vector, pad_width, iaxis, kwargs):
            vector[:pad_width[0]] = 0
            vector[-pad_width[1]:] = 0
            return vector
        ny, nx = self.pixel_images[0].shape
        contour = np.zeros((ny, nx))
        contour[np.where(self.aperture)] = 1
        contour = np.lib.pad(contour, 1, PadWithZeros)
        highres = zoom(contour, 100, order=0, mode='nearest')
        extent = np.array([-1, nx, -1, ny])

        # Plot first, mid, and last TPF image
        title = ['start', 'mid', 'end']
        for i, image in enumerate(self.pixel_images):
            ax = axes[i]
            ax.imshow(image, aspect='auto',
                      interpolation='nearest', cmap=plasma)
            ax.contour(highres, levels=[0.5], extent=extent,
                       origin='lower', colors='r', linewidths=1)

            # Check for saturated columns
            for x in range(self.aperture.shape[0]):
                for y in range(self.aperture.shape[1]):
                    if self.aperture[x][y] == AP_SATURATED_PIXEL:
                        ax.fill([y - 0.5, y + 0.5, y + 0.5, y - 0.5],
                                [x - 0.5, x - 0.5, x + 0.5, x + 0.5],
                                fill=False, hatch='xxxxx', color='r', lw=0)

            ax.axis('off')
            ax.set_xlim(-0.7, nx - 0.3)
            ax.set_ylim(-0.7, ny - 0.3)
            ax.annotate(title[i], xy=(0.5, 0.975), xycoords='axes fraction',
                        ha='center', va='top', size=labelsize, color='w')
            if i == 1:
                for source in self.nearby:
                    ax.annotate('%.1f' % source['mag'],
                                xy=(source['x'] - source['x0'],
                                    source['y'] - source['y0']),
                                ha='center', va='center', size=labelsize - 2,
                                color='w', fontweight='bold')

        # Plot hi res image
        if self.hires is not None:
            ax = axes[-1]
            ax.imshow(self.hires, aspect='auto',
                      extent=(-0.5, nx - 0.5, -0.5, ny - 0.5),
                      interpolation='bicubic', cmap=plasma)
            ax.contour(highres, levels=[0.5], extent=extent,
                       origin='lower', colors='r', linewidths=1)
            ax.axis('off')
            ax.set_xlim(-0.7, nx - 0.3)
            ax.set_ylim(-0.7, ny - 0.3)
            ax.annotate('hires', xy=(0.5, 0.975), xycoords='axes fraction',
                        ha='center', va='top', size=labelsize, color='w')
        else:
            ax = axes[-1]
            ax.axis('off')

    def search(self, pos_tol=2.5, neg_tol=50., clobber=False,
               name='search', **kwargs):
        '''

        '''

        log.info("Searching for transits...")
        fname = os.path.join(self.dir, self.name + '_%s.npz' % name)
        pname = os.path.join(self.dir, self.name + '_%s.pdf' % name)

        # Compute
        if not os.path.exists(fname) or clobber:
            time, depth, vardepth, delchisq = Search(
                self, pos_tol=pos_tol, neg_tol=neg_tol, **kwargs)
            data = np.vstack([time, depth, vardepth, delchisq]).T
            header = "TIME, DEPTH, VARDEPTH, DELTACHISQ"
            np.savetxt(fname, data, fmt=str('%.10e'), header=header)
        else:
            time, depth, vardepth, delchisq = np.loadtxt(
                fname, unpack=True, skiprows=1)

        # Plot
        if not os.path.exists(pname) or clobber:
            fig, ax = pl.subplots(1, figsize=(10, 4))
            ax.plot(time, delchisq, lw=1)
            ax.set_ylabel(r'$\Delta \chi^2$', fontsize=18)
            ax.set_xlabel('Time (days)', fontsize=18)
            ax.set_xlim(time[0], time[-1])
            fig.savefig(pname, bbox_inches='tight')
            pl.close()

        return time, depth, vardepth, delchisq

    def overfit(self, tau=None, plot=True, clobber=False, w=9, **kwargs):
        r"""
        Compute the masked & unmasked overfitting metrics for the light curve.

        This routine injects a transit model given by `tau` at every cadence
        in the light curve and recovers the transit depth when (1) leaving
        the transit unmasked and (2) masking the transit prior to performing
        regression.

        :param tau: A function or callable that accepts two arguments, \
               `time` and `t0`, and returns an array corresponding to a \
               zero-mean, unit depth transit model centered at \
               `t0` and evaluated at `time`. \
               The easiest way to provide this is to use an instance of \
               :py:class:`everest.transit.TransitShape`. Default is \
               :py:class:`everest.transit.TransitShape(dur=0.1)`, a transit \
               with solar-like limb darkening and a duratio of 0.1 days.
        :param bool plot: Plot the results as a PDF? Default :py:obj:`True`
        :param bool clobber: Overwrite the results if present? Default \
               :py:obj:`False`
        :param int w: The size of the masking window in cadences for \
               computing the masked overfitting metric. Default `9` \
               (about 4.5 hours for `K2` long cadence).

        :returns: An instance of `everest.basecamp.Overfitting`.
        """
        fname = os.path.join(self.dir, self.name + '_overfit.npz')
        figname = os.path.join(self.dir, self.name)

        # Compute
        if not os.path.exists(fname) or clobber:

            # Baseline
            med = np.nanmedian(self.fraw)

            # Default transit model
            if tau is None:
                tau = TransitShape(dur=0.1)

            # The overfitting metrics
            O1 = [None for brkpt in self.breakpoints]
            O2 = [None for brkpt in self.breakpoints]
            O3 = [None for brkpt in self.breakpoints]
            O4 = [None for brkpt in self.breakpoints]
            O5 = [None for brkpt in self.breakpoints]

            # Loop over all chunks
            for b, brkpt in enumerate(self.breakpoints):

                # Masks for current chunk
                m = self.get_masked_chunk(b, pad=False)
                time = self.time[m]
                ferr = self.fraw_err[m] / med
                y = self.fraw[m] / med - 1

                # The metrics we're computing here
                O1[b] = np.zeros(len(y)) * np.nan
                O2[b] = np.zeros(len(y)) * np.nan
                O3[b] = np.zeros(len(y)) * np.nan
                O4[b] = np.zeros(len(y)) * np.nan
                O5[b] = np.zeros(len(y)) * np.nan

                # Compute the astrophysical covariance and its inverse
                log.info("Computing the covariance...")
                if self.kernel == 'Basic':
                    wh, am, ta = self.kernel_params
                    wh /= med
                    am /= med
                    kernel_params = [wh, am, ta]
                elif self.kernel == 'QuasiPeriodic':
                    wh, am, ga, pe = self.kernel_params
                    wh /= med
                    am /= med
                    kernel_params = [wh, am, ga, pe]
                K = GetCovariance(self.kernel, kernel_params, time, ferr)
                Kinv = cho_solve((cholesky(K), False), np.eye(len(time)))

                # Loop over all orders
                log.info("Computing some large matrices...")
                X = [None for n in range(self.pld_order)]
                XL = [None for n in range(self.pld_order)]
                XLX = [None for n in range(self.pld_order)]
                for n in range(self.pld_order):
                    if (self.lam_idx >= n) and (self.lam[b][n] is not None):
                        X[n] = self.X(n, m, **kwargs)
                        XL[n] = (self.lam[b][n] / med ** 2) * X[n]
                        XLX[n] = np.dot(XL[n], X[n].T)
                X = np.hstack(X)
                XL = np.hstack(XL)
                XLX = np.sum(XLX, axis=0)

                # The full covariance
                C = XLX + K

                # The unmasked linear problem
                log.info("Solving the unmasked linear problem...")
                m = np.dot(XLX, np.linalg.solve(C, y))
                m -= np.nanmedian(m)
                f = y - m
                R = np.linalg.solve(C, XLX.T).T

                # The masked linear problem
                log.info("Solving the masked linear problem...")
                A = MaskSolve(C, y, w=w)

                # Now loop through and compute the metric
                log.info("Computing the overfitting metrics...")
                for n in prange(len(y)):

                    #
                    # *** Unmasked overfitting metric ***
                    #

                    # Evaluate the sparse transit model
                    TAU = tau(time, t0=time[n])
                    i = np.where(TAU < 0)[0]
                    TAU = TAU.reshape(-1, 1)

                    # Fast sparse algebra
                    AA = np.dot(np.dot(TAU[i].T, Kinv[i, :][:, i]), TAU[i])
                    BB = np.dot(TAU[i].T, Kinv[i, :])
                    CC = TAU - np.dot(R[:, i], TAU[i])
                    O1[b][n] = AA
                    O2[b][n] = np.dot(BB, CC)
                    O3[b][n] = np.dot(BB, f)
                    O4[b][n] = np.dot(BB, y)

                    #
                    # *** Masked overfitting metric ***
                    #

                    # The current mask and mask centerpoint
                    mask = np.arange(n, n + w)
                    j = n + (w + 1) // 2 - 1
                    if j >= len(y) - w:
                        continue

                    # The regularized design matrix
                    # This is the same as
                    # XLmX[:, n - 1] = \
                    #   np.dot(XL, np.delete(X, mask, axis=0).T)[:, n - 1]
                    if n == 0:
                        XLmX = np.dot(XL, np.delete(X, mask, axis=0).T)
                    else:
                        XLmX[:, n - 1] = np.dot(XL, X[n - 1, :].T)

                    # The linear solution to this step
                    m = np.dot(XLmX, A[n])

                    # Evaluate the sparse transit model
                    TAU = tau(time, t0=time[j])
                    i = np.where(TAU < 0)[0]
                    TAU = TAU[i].reshape(-1, 1)

                    # Dot the transit model in
                    den = np.dot(np.dot(TAU.T, Kinv[i, :][:, i]), TAU)
                    num = np.dot(TAU.T, Kinv[i, :])

                    # Compute the overfitting metric
                    # Divide this number by a depth
                    # to get the overfitting for that
                    # particular depth.
                    O5[b][j] = -np.dot(num, y - m) / den

            # Save!
            np.savez(fname, O1=O1, O2=O2, O3=O3, O4=O4, O5=O5)

        else:

            data = np.load(fname)
            O1 = data['O1']
            O2 = data['O2']
            O3 = data['O3']
            O4 = data['O4']
            O5 = data['O5']

        # Plot
        if plot and (clobber or not os.path.exists(figname + '_overfit.pdf')):

            log.info("Plotting the overfitting metrics...")

            # Masked time array
            time = self.apply_mask(self.time)

            # Plot the final corrected light curve
            ovr = OVERFIT()
            self.plot_info(ovr)

            # Loop over the two metrics
            for kind, axes, axesh in zip(['unmasked', 'masked'],
                                         [ovr.axes1, ovr.axes2],
                                         [ovr.axes1h, ovr.axes2h]):

                # Loop over three depths
                for depth, ax, axh in zip([0.01, 0.001, 0.0001], axes, axesh):

                    # Get the metric
                    if kind == 'unmasked':
                        metric = 1 - (np.hstack(O2) +
                                      np.hstack(O3) / depth) / np.hstack(O1)
                        color = 'r'
                    elif kind == 'masked':
                        metric = np.hstack(O5) / depth
                        color = 'b'
                    else:
                        raise ValueError("Invalid metric.")

                    # Median and median absolute deviation
                    med = np.nanmedian(metric)
                    mad = np.nanmedian(np.abs(metric - med))

                    # Plot the metric as a function of time
                    ax.plot(time, metric, 'k.', alpha=0.5, ms=2)
                    ax.plot(time, metric, 'k-', alpha=0.1, lw=0.5)
                    ylim = (-0.2, 1.0)
                    ax.margins(0, None)
                    ax.axhline(0, color='k', lw=1, alpha=0.5)
                    ax.set_ylim(*ylim)
                    if kind == 'masked' and depth == 0.0001:
                        ax.set_xlabel('Time (days)', fontsize=14)
                    else:
                        ax.set_xticklabels([])

                    # Plot the histogram
                    rng = (max(ylim[0], np.nanmin(metric)),
                           min(ylim[1], np.nanmax(metric)))
                    axh.hist(metric, bins=30, range=rng,
                             orientation="horizontal",
                             histtype="step", fill=False, color='k')
                    axh.axhline(med, color=color, ls='-', lw=1)
                    axh.axhspan(med - mad, med + mad, color=color, alpha=0.1)
                    axh.axhline(0, color='k', lw=1, alpha=0.5)
                    axh.yaxis.tick_right()
                    axh.set_ylim(*ax.get_ylim())
                    axh.set_xticklabels([])

                    bbox = dict(fc="w", ec="1", alpha=0.5)
                    info = r"$\mathrm{med}=%.3f$" % med + \
                        "\n" + r"$\mathrm{mad}=%.3f$" % mad
                    axh.annotate(info, xy=(0.1, 0.925),
                                 xycoords='axes fraction',
                                 ha="left", va="top", bbox=bbox, color=color)

                    bbox = dict(fc="w", ec="1", alpha=0.95)
                    ax.annotate("%s overfitting metric" % kind,
                                xy=(1-0.035, 0.92),
                                xycoords='axes fraction',
                                ha='right', va='top',
                                bbox=bbox, color=color)

            pl.figtext(0.025, 0.77, "depth = 0.01", rotation=90,
                       ha='left', va='center', fontsize=18)
            pl.figtext(0.025, 0.48, "depth = 0.001", rotation=90,
                       ha='left', va='center', fontsize=18)
            pl.figtext(0.025, 0.19, "depth = 0.0001", rotation=90,
                       ha='left', va='center', fontsize=18)

            ovr.fig.savefig(figname + '_overfit.pdf')
            log.info("Saved plot to %s_overfit.pdf" % figname)
            pl.close()

        return Overfitting(O1, O2, O3, O4, O5, figname + '_overfit.pdf')

    def lnlike(self, model, refactor=False, pos_tol=2.5, neg_tol=50.,
               full_output=False):
        r"""
        Return the likelihood of the astrophysical model `model`.

        Returns the likelihood of `model` marginalized over the PLD model.

        :param ndarray model: A vector of the same shape as `self.time` \
               corresponding to the astrophysical model.
        :param bool refactor: Re-compute the Cholesky decomposition? This \
               typically does not need to be done, except when the PLD \
               model changes. Default :py:obj:`False`.
        :param float pos_tol: the positive (i.e., above the median) \
               outlier tolerance in standard deviations.
        :param float neg_tol: the negative (i.e., below the median) \
               outlier tolerance in standard deviations.
        :param bool full_output: If :py:obj:`True`, returns the maximum \
               likelihood model amplitude and the variance on the amplitude \
               in addition to the log-likelihood. In the case of a transit \
               model, these are the transit depth and depth variance. Default \
               :py:obj:`False`.
        """
        lnl = 0

        # Re-factorize the Cholesky decomposition?
        try:
            self._ll_info
        except AttributeError:
            refactor = True
        if refactor:

            # Smooth the light curve and reset the outlier mask
            t = np.delete(self.time,
                          np.concatenate([self.nanmask, self.badmask]))
            f = np.delete(self.flux,
                          np.concatenate([self.nanmask, self.badmask]))
            f = SavGol(f)
            med = np.nanmedian(f)
            MAD = 1.4826 * np.nanmedian(np.abs(f - med))
            pos_inds = np.where((f > med + pos_tol * MAD))[0]
            pos_inds = np.array([np.argmax(self.time == t[i])
                                 for i in pos_inds])
            MAD = 1.4826 * np.nanmedian(np.abs(f - med))
            neg_inds = np.where((f < med - neg_tol * MAD))[0]
            neg_inds = np.array([np.argmax(self.time == t[i])
                                 for i in neg_inds])
            outmask = np.array(self.outmask)
            transitmask = np.array(self.transitmask)
            self.outmask = np.concatenate([neg_inds, pos_inds])
            self.transitmask = np.array([], dtype=int)

            # Now re-factorize the Cholesky decomposition
            self._ll_info = [None for b in self.breakpoints]
            for b, brkpt in enumerate(self.breakpoints):

                # Masks for current chunk
                m = self.get_masked_chunk(b, pad=False)

                # This block of the masked covariance matrix
                K = GetCovariance(self.kernel, self.kernel_params,
                                  self.time[m], self.fraw_err[m])

                # The masked X.L.X^T term
                A = np.zeros((len(m), len(m)))
                for n in range(self.pld_order):
                    XM = self.X(n, m)
                    A += self.lam[b][n] * np.dot(XM, XM.T)
                K += A
                self._ll_info[b] = [cho_factor(K), m]

            # Reset the outlier masks
            self.outmask = outmask
            self.transitmask = transitmask

        # Compute the likelihood for each chunk
        amp = [None for b in self.breakpoints]
        var = [None for b in self.breakpoints]
        for b, brkpt in enumerate(self.breakpoints):
            # Get the inverse covariance and the mask
            CDK = self._ll_info[b][0]
            m = self._ll_info[b][1]
            # Compute the maximum likelihood model amplitude
            # (for transits, this is the transit depth)
            var[b] = 1. / np.dot(model[m], cho_solve(CDK, model[m]))
            amp[b] = var[b] * np.dot(model[m], cho_solve(CDK, self.fraw[m]))
            # Compute the residual
            r = self.fraw[m] - amp[b] * model[m]
            # Finally, compute the likelihood
            lnl += -0.5 * np.dot(r, cho_solve(CDK, r))

        if full_output:
            # We need to multiply the Gaussians for all chunks to get the
            # amplitude and amplitude variance for the entire dataset
            vari = var[0]
            ampi = amp[0]
            for v, a in zip(var[1:], amp[1:]):
                ampi = (ampi * v + a * vari) / (vari + v)
                vari = vari * v / (vari + v)
            med = np.nanmedian(self.fraw)
            return lnl, ampi / med, vari / med ** 2
        else:
            return lnl
