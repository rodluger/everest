#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`transit.py` - Transit models
-------------------------------------
These are routines used to generate a transit model, primarily for
transit injection/recovery tests. These are wrappers around
:py:func:`pysyzygy.Transit`, with the added feature that
the transit :py:obj:`depth` and the transit :py:obj:`duration` can be specified
as input variables (as opposed to the planet-star radius ratio
and the stellar density, which :py:mod:`pysyzygy` expects).

'''

from __future__ import division, print_function, absolute_import, \
     unicode_literals
import numpy as np
import matplotlib.pyplot as pl
try:
    import pysyzygy as ps
except:
    ps = None
from scipy.optimize import fmin
import logging
log = logging.getLogger(__name__)


class TransitModel(object):
    r"""
    A transit model for use in regression. Note that the transit model
    is linear in the transit depth, so the depth (or planet/star radius
    ratio) is not a free parameter here!

    :param str name: The name of the planet, such as "b".
    :param bool single: Is this a single transit? Default :py:obj:`False` \
        (transit is assumed to be periodic at period given by `per`).
    :param float sig_RpRs: The standard deviation of the prior on the \
        planet/star radius ratio. Used only when the transit model is \
        included in the design matrix. Default `0.001`
    :param kwargs: Any of the :py:mod:`pysyzygy.Transit` keyword arguments. \
        These include:

    - **b** or **bcirc** - The (circular) impact parameter. Default `0.`
    - **MpMs** - The planet-star mass ratio. Default `0.`
    - **per** - The planet orbital period in days. Default `10.`
    - **rhos** or **aRs** - The stellar density in `g/cm^3` or the semi-major \
        axis-stellar radius ratio. Default is `rhos = 1.4`, the density of \
        the Sun
    - **ecc** and **w** or **esw** and **ecw** - The eccentricity and the \
        longitude of pericenter in radians, or the two eccentricity vectors. \
        Default is `ecc = 0.` and `w = 0.`
    - **t0** or **times** - The time of first transit, or the time of each \
        of the transits (in case they are not periodic) in days. Default \
        is `t0 = 0.`
    - **u1** and **u2** or **q1** and **q2** - The quadratic limb darkening \
        parameters (u1, u2) or the modified quadratic limb darkening \
        parameters (q1, q2) from \
        `Kipping (2013) <http://dx.doi.org/10.1093/mnras/stt1435>`_. \
        Default is `u1 = 0.40` and `u2 = 0.26`
    - **exptime** - The exposure time in days for binning the model. \
        Default `ps.KEPLONGEXP`
    - **fullorbit** - Compute the orbital parameters for the entire orbit? \
        Only useful if you're interested in the full arrays of orbital \
        parameters. Default `False`
    - **maxpts** - Maximum number of points in the model. Increase \
        this if you're getting errors. Default `10,000`
    - **exppts** - The number of exposure points per cadence when \
        binning the model. Default `50`
    - **keptol** - The tolerance of the Kepler solver. Default `1.e-15`
    - **maxkepiter** - Maximum number of iterations in the Kepler solver. \
        Default `100`

    """

    def __init__(self, name, single=False, sig_RpRs=0.001, **kwargs):
        """Instantiate the pysyzygy model."""
        # The planet/transit model ID
        assert type(name) is str, "Arg `name` must be a string."
        self.name = name

        # Skip if pysyzygy not installed
        if ps is None:
            return

        # The transit model
        self._transit = ps.Transit(**kwargs)

        # Compute the depth
        times = kwargs.get('times', None)
        if times is not None:
            self.t0 = times[0]
        else:
            self.t0 = kwargs.get('t0', 0.)
        self.per = kwargs.get('per', 10.)
        self.single = single
        self.depth = (1. - self._transit([self.t0]))[0]

        # Approximate variance on the depth
        self.var_depth = (2 * sig_RpRs) ** 2

        # Save the kwargs
        self.params = kwargs

    def __call__(self, time):
        """Return the model evaluated at `time`."""
        if ps is None:
            raise Exception("Unable to import `pysyzygy`.")

        model = (self._transit(time) - 1) / self.depth

        # Single transit?
        if self.single:
            model[np.where(np.abs(time - self.t0) > self.per / 5.)[0]] = 0.

        return model


def Get_RpRs(d, **kwargs):
    '''
    Returns the value of the planet radius over the stellar radius
    for a given depth :py:obj:`d`, given
    the :py:class:`everest.pysyzygy` transit :py:obj:`kwargs`.

    '''
    if ps is None:
            raise Exception("Unable to import `pysyzygy`.")

    def Depth(RpRs, **kwargs):
        return 1 - ps.Transit(RpRs=RpRs, **kwargs)([kwargs.get('t0', 0.)])

    def DiffSq(r):
        return 1.e10 * (d - Depth(r, **kwargs)) ** 2

    return fmin(DiffSq, [np.sqrt(d)], disp=False)


def Get_rhos(dur, **kwargs):
    '''
    Returns the value of the stellar density for a given transit
    duration :py:obj:`dur`, given
    the :py:class:`everest.pysyzygy` transit :py:obj:`kwargs`.

    '''
    if ps is None:
            raise Exception("Unable to import `pysyzygy`.")

    assert dur >= 0.01 and dur <= 0.5, "Invalid value for the duration."

    def Dur(rhos, **kwargs):
        t0 = kwargs.get('t0', 0.)
        time = np.linspace(t0 - 0.5, t0 + 0.5, 1000)
        try:
            t = time[np.where(ps.Transit(rhos=rhos, **kwargs)(time) < 1)]
        except:
            return 0.
        return t[-1] - t[0]

    def DiffSq(rhos):
        return (dur - Dur(rhos, **kwargs)) ** 2

    return fmin(DiffSq, [0.2], disp=False)


def Transit(time, t0=0., dur=0.1, per=3.56789, depth=0.001, **kwargs):
    '''
    A `Mandel-Agol <http://adsabs.harvard.edu/abs/2002ApJ...580L.171M>`_
    transit model, but with the depth and the duration as primary
    input variables.

    :param numpy.ndarray time: The time array
    :param float t0: The time of first transit in units of \
           :py:obj:`BJD` - 2454833.
    :param float dur: The transit duration in days. Don't go too crazy on \
           this one -- very small or very large values will break the \
           inverter. Default 0.1
    :param float per: The orbital period in days. Default 3.56789
    :param float depth: The fractional transit depth. Default 0.001
    :param dict kwargs: Any additional keyword arguments, passed directly \
           to :py:func:`pysyzygy.Transit`
    :returns tmod: The transit model evaluated at the same times as the \
                   :py:obj:`time` array

    '''
    if ps is None:
        raise Exception("Unable to import `pysyzygy`.")

    # Note that rhos can affect RpRs, so we should really do this iteratively,
    # but the effect is pretty negligible!
    RpRs = Get_RpRs(depth, t0=t0, per=per, **kwargs)
    rhos = Get_rhos(dur, t0=t0, per=per, **kwargs)
    return ps.Transit(t0=t0, per=per, RpRs=RpRs, rhos=rhos, **kwargs)(time)


class TransitShape(object):
    r"""
    Return a single transit model used for injection/recovery tests.

    A `Mandel-Agol <http://adsabs.harvard.edu/abs/2002ApJ...580L.171M>`_
    transit model, but with the depth and the duration as primary
    input variables. Note that the time(s) of transit is not an
    input parameter, nor is the period.

    :param float depth: The transit depth. Set this to 1. for transit \
           injection/recovery tests. Default `1`
    :param float dur: The transit duration in days. Default `0.1`
    :param kwargs: Any additional keyword arguments, passed directly \
           to :py:func:`pysyzygy.Transit`
    """

    def __init__(self, depth=1, dur=0.1, **kwargs):
        """Initialize the pysyzygy model."""
        if ps is None:
            return
        kwargs.pop('t0', None)
        kwargs.pop('times', None)

        # Update kwargs with correct duration
        kwargs.update({'per': 3.56789})
        kwargs.update({'rhos': Get_rhos(dur, **kwargs)})

        # Transit window size w/ padding
        window = dur * 3
        t = np.linspace(-window / 2, window / 2, 5000)

        # Construct a transit model
        trn = ps.Transit(t0=0., **kwargs)
        transit_model = trn(t)
        transit_model -= 1
        transit_model *= depth / (1 - trn([0.])[0])
        self.x = t
        self.y = transit_model

    def __call__(self, time, t0=0.):
        """Evalaute the transit model."""
        if ps is None:
            raise Exception("Unable to import `pysyzygy`.")
        return np.interp(time, self.x + t0, self.y)
