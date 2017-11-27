#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`dvs.py` - Data Validation Summary
------------------------------------------

Code for handling the "Data Validation Summary" plot.

'''

from __future__ import division, print_function, absolute_import, \
     unicode_literals
import matplotlib.pyplot as pl
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition


class Frame(object):
    '''
    A not-so-elegant object that adds an inset axis at a given
    position within a matplotlib axes instance when called.

    '''

    def __init__(self, fig, ax, pos=[0, 0, 1, 1]):
        '''

        '''

        self.fig = fig
        self.ax = ax
        self.pos = pos

    def __call__(self, pos=None, on=True):
        '''

        '''

        if pos is None:
            pos = self.pos
        res = []
        for axis in np.atleast_1d(self.ax):
            ax = self.fig.add_subplot(111, label=np.random.randn())
            ax.set_axes_locator(InsetPosition(axis, pos))
            for tick in ax.get_xticklabels() + ax.get_yticklabels():
                tick.set_fontsize(5)
            if not on:
                ax.axis('off')
            res.append(ax)
        if len(res) == 1:
            # This is a single axis
            return res[0]
        else:
            # This is a list of axes
            return res


class DVS(object):
    '''
    The "Data Validation Summary" figure container.

    :param int nchunks: The number of light curve segments. Default 2
    :param int pld_order: The PLD order. Default 3

    '''

    def __init__(self, nchunks=2, pld_order=3):
        '''

        '''

        if pld_order <= 3:
            hght = 28
            nrows = 160
        else:
            hght = 32
            nrows = 174 + hght * (pld_order - 3)
        self.fig = pl.figure(figsize=(8.5, 11))
        self.fig.subplots_adjust(
            left=0.025 * (11 / 8.5), right=1 - 0.025 * (11 / 8.5), top=0.975,
            bottom=0.025)

        def GetFrame(y, x, dx, dy):
            return Frame(self.fig, pl.subplot2grid((nrows, 160), (y, x),
                         colspan=dx, rowspan=dy))
        self.title_left = GetFrame(0, 6, 44, 10)
        self.title_center = GetFrame(0, 50, 66, 10)
        self.title_right = GetFrame(0, 116, 44, 10)
        self.body_top_left = GetFrame(12, 6, 102, 26)
        self.body_top_right = [GetFrame(12, 116, 21, 26),
                               GetFrame(12, 139, 21, 26),
                               GetFrame(12 + hght, 116, 21, 26),
                               GetFrame(12 + hght, 139, 21, 26)]
        self.body_left = [GetFrame(12 + hght * n, 6, 102, 26)
                          for n in range(1, 2 + pld_order)]
        if (nchunks == 2) or (nchunks > 3):
            self.body_right = [Frame(self.fig,
                                     [pl.subplot2grid((nrows, 160),
                                      (12 + hght * n, 116),
                                      colspan=44, rowspan=13),
                                      pl.subplot2grid((nrows, 160),
                                      (25 + hght * n, 116), colspan=44,
                                      rowspan=13)])
                               for n in range(2, 2 + pld_order)]
        elif nchunks == 3:
            self.body_right = [Frame(self.fig,
                                     [pl.subplot2grid((nrows, 160),
                                                      (12 + hght * n, 116),
                                                      colspan=44,
                                                      rowspan=9),
                                      pl.subplot2grid((nrows, 160),
                                                      (21 + hght * n, 116),
                                                      colspan=44,
                                                      rowspan=8),
                                      pl.subplot2grid((nrows, 160),
                                                      (29 + hght * n, 116),
                                                      colspan=44,
                                                      rowspan=9)])
                               for n in range(2, 2 + pld_order)]
        else:
            self.body_right = [GetFrame(12 + hght * n, 116, 44, 26)
                               for n in range(2, 2 + pld_order)]
        self.footer_left = GetFrame(nrows - 6, 6, 44, 6)
        self.footer_center = GetFrame(nrows - 6, 50, 66, 6)
        self.footer_right = GetFrame(nrows - 6, 116, 44, 6)
        for ax in self.fig.get_axes():
            ax.axis('off')
        self.tcount = 0
        self.lcount = 0
        self.rcount = 0

    def title(self):
        '''
        Returns the axis instance where the title will be printed

        '''

        return self.title_left(on=False), self.title_center(on=False), \
               self.title_right(on=False)

    def footer(self):
        '''
        Returns the axis instance where the footer will be printed

        '''

        return self.footer_left(on=False), self.footer_center(on=False), \
               self.footer_right(on=False)

    def top_right(self):
        '''
        Returns the axis instance at the top right of the page,
        where the postage stamp and aperture is displayed

        '''

        res = self.body_top_right[self.tcount]()
        self.tcount += 1
        return res

    def top_left(self):
        '''
        Returns the axis instance at the top left of the page,
        where the final de-trended light curve is displayed

        '''

        return self.body_top_left()

    def left(self):
        '''
        Returns the current axis instance on the left side of
        the page where each successive light curve is displayed

        '''

        res = self.body_left[self.lcount]()
        self.lcount += 1
        return res

    def right(self):
        '''
        Returns the current axis instance on the right side of the
        page, where cross-validation information is displayed

        '''

        res = self.body_right[self.rcount]()
        self.rcount += 1
        return res


class CBV(object):
    '''

    '''

    def __init__(self):
        '''

        '''

        self.fig = pl.figure(figsize=(8.5, 11))
        self.fig.subplots_adjust(
            left=0.025 * (11 / 8.5), right=1 - 0.025 * (11 / 8.5),
            top=0.975, bottom=0.025)

        def GetFrame(y, x, dx, dy):
            return Frame(self.fig, pl.subplot2grid((160, 160), (y, x),
                         colspan=dx, rowspan=dy))
        self.title_left = GetFrame(0, 6, 44, 10)
        self.title_center = GetFrame(0, 50, 66, 10)
        self.title_right = GetFrame(0, 116, 44, 10)
        self._body = [GetFrame(12, 6, 148, 42),
                      GetFrame(62, 6, 148, 42),
                      GetFrame(112, 6, 148, 42)]
        for ax in self.fig.get_axes():
            ax.axis('off')
        self.bcount = 0

    def title(self):
        '''
        Returns the axis instance where the title will be printed

        '''

        return self.title_left(on=False), self.title_center(on=False), \
            self.title_right(on=False)

    def body(self):
        '''
        Returns the axis instance where the light curves will be shown

        '''

        res = self._body[self.bcount]()
        self.bcount += 1
        return res

class OVERFIT(object):
    '''

    '''

    def __init__(self):
        '''

        '''

        self.fig = pl.figure(figsize=(8.5, 11))
        self.fig.subplots_adjust(
            left=0.025 * (11 / 8.5), right=1 - 0.025 * (11 / 8.5),
            top=0.975, bottom=0.025, hspace=0.5, wspace=0.5)
        def GetFrame(y, x, dx, dy):
            return Frame(self.fig, pl.subplot2grid((160, 160), (y, x),
                         colspan=dx, rowspan=dy))
        self.title_left = GetFrame(0, 6, 44, 10)
        self.title_center = GetFrame(0, 50, 66, 10)
        self.title_right = GetFrame(0, 116, 44, 10)
        for ax in self.fig.get_axes():
            ax.axis('off')
        kw = dict(colspan=40, rowspan=10)
        kwh = dict(colspan=10, rowspan=10)
        self.axes1 = [pl.subplot2grid((70, 60), (5, 5), **kw),
                      pl.subplot2grid((70, 60), (26, 5), **kw),
                      pl.subplot2grid((70, 60), (47, 5), **kw)]
        self.axes1h = [pl.subplot2grid((70, 60), (5, 45), **kwh),
                       pl.subplot2grid((70, 60), (26, 45), **kwh),
                       pl.subplot2grid((70, 60), (47, 45), **kwh)]
        self.axes2 = [pl.subplot2grid((70, 60), (15, 5), **kw),
                      pl.subplot2grid((70, 60), (36, 5), **kw),
                      pl.subplot2grid((70, 60), (57, 5), **kw)]
        self.axes2h = [pl.subplot2grid((70, 60), (15, 45), **kwh),
                       pl.subplot2grid((70, 60), (36, 45), **kwh),
                       pl.subplot2grid((70, 60), (57, 45), **kwh)]
        for ax in [self.axes1, self.axes1h, self.axes2, self.axes2h]:
            for axis in ax:
                axis.tick_params(direction='in')

    def title(self):
        '''
        Returns the axis instance where the title will be printed

        '''

        return self.title_left(on=False), self.title_center(on=False), \
               self.title_right(on=False)
