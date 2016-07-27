#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`transit_plot.py` - A quick light curve
-----------------------------------------------

Plots a simple light curve using the :py:func:`pysyzygy.PlotTransit`
function, which also shows the orbit, the limb-darkened stellar disk,
and a table with all the planet parameters.

.. figure:: ../pysyzygy/img/transit_plot.png
    :width: 600px
    :align: center
    :height: 100px
    :alt: alternate text
    :figclass: align-center

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import pysyzygy as ps
import numpy as np
import matplotlib.pyplot as pl

if __name__ == '__main__':

  # Transit model kwargs
  kwargs = dict(RpRs = 0.1, b = 0.5, 
                ecc = 0.5, w = 0.,
                MpMs = 0., rhos = 1.4,
                per = 5., t0 = 0.,
                q1 = 0.45, q2 = 0.3,
                maxpts = 20000,
                exptime = ps.KEPLONGEXP)

  fig = ps.PlotTransit(**kwargs)
  pl.show()