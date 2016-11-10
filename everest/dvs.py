#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`dvs.py`
----------------

Code for handling the "Data Validation Summary" plot.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import matplotlib.pyplot as pl
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

class Frame(object):
  '''
  
  '''
  
  def __init__(self, fig, ax, pos = [0, 0, 1, 1]):
    '''
    
    '''
    
    self.fig = fig
    self.ax = ax
    self.pos = pos

  def __call__(self, pos = None, on = True):
    '''
  
    '''
    
    if pos is None:
      pos = self.pos
    res = []
    for axis in np.atleast_1d(self.ax):
      ax = self.fig.add_subplot(111, label = np.random.randn())
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
    
class DVS1(object):
  '''
  
  '''
  
  def __init__(self, breakpoint = False, pld_order = 3):
    '''
    
    '''
    
    if pld_order <= 3:
      hght = 28
      nrows = 160
    else:
      hght = 32
      nrows = 174 + hght * (pld_order - 3)  
    self.fig = pl.figure(figsize = (8.5, 11))
    self.fig.subplots_adjust(left = 0.025 * (11/8.5), right = 1 - 0.025 * (11/8.5), top = 0.975, bottom = 0.025)
    def GetFrame(y, x, dx, dy):
      return Frame(self.fig, pl.subplot2grid((nrows, 160), (y, x), colspan = dx, rowspan = dy))
    self.title_left = GetFrame(0, 6, 44, 10)
    self.title_center = GetFrame(0, 50, 66, 10)
    self.title_right = GetFrame(0, 116, 44, 10)
    self.body_top_left = GetFrame(12, 6, 102, 26)
    self.body_top_right = [GetFrame(12, 116, 21, 26), GetFrame(12, 139, 21, 26), GetFrame(12 + hght, 116, 21, 26), GetFrame(12 + hght, 139, 21, 26)]
    self.body_left = [GetFrame(12 + hght * n, 6, 102, 26) for n in range(1,2 + pld_order)]
    if not breakpoint:
      self.body_right = [GetFrame(12 + hght * n, 116, 44, 26) for n in range(2,2 + pld_order)]
    else:
      self.body_right = [Frame(self.fig,[pl.subplot2grid((nrows, 160), (12 + hght * n, 116), colspan=44, rowspan=13), pl.subplot2grid((nrows, 160), (25 + hght * n, 116), colspan=44, rowspan=13)]) for n in range(2, 2 + pld_order)]                
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
    
    '''
    
    return self.title_left(on = False), self.title_center(on = False), self.title_right(on = False)

  def footer(self):
    '''
    
    '''
    
    return self.footer_left(on = False), self.footer_center(on = False), self.footer_right(on = False)
  
  def top_right(self):
    '''
    
    '''
    
    res = self.body_top_right[self.tcount]()
    self.tcount += 1
    return res
  
  def top_left(self):
    '''
    
    '''
    
    return self.body_top_left()
  
  def left(self):
    '''
    
    '''
    
    res = self.body_left[self.lcount]()
    self.lcount += 1
    return res

  def right(self):
    '''
    
    '''
    
    res = self.body_right[self.rcount]()
    self.rcount += 1
    return res

class DVS2(object):
  '''
  
  '''
  
  def __init__(self, breakpoint = False, pld_order = 3):
    '''
    
    '''
    
    self.fig = pl.figure(figsize = (8.5, 11))
    self.fig.subplots_adjust(left = 0.025 * (11/8.5), right = 1 - 0.025 * (11/8.5), top = 0.975, bottom = 0.025)
    def GetFrame(y, x, dx, dy):
      return Frame(self.fig, pl.subplot2grid((160, 160), (y, x), colspan = dx, rowspan = dy))
    self.title_left = GetFrame(0, 6, 44, 10)
    self.title_center = GetFrame(0, 50, 66, 10)
    self.title_right = GetFrame(0, 116, 44, 10)
    self.body_top_right = [GetFrame(12, 116, 21, 26), GetFrame(12, 139, 21, 26), GetFrame(40, 116, 21, 26), GetFrame(40, 139, 21, 26)]
    self.lc1 = GetFrame(12, 6, 102, 26)
    self.lc2 = GetFrame(12 + 28, 6, 102, 26)
    self.footer_left = GetFrame(154, 6, 44, 6)
    self.footer_center = GetFrame(154, 50, 66, 6)
    self.footer_right = GetFrame(154, 116, 44, 6)
    for ax in self.fig.get_axes():
      ax.axis('off')
    self.tcount = 0
    if breakpoint:
      self.nseg = 2
    else:
      self.nseg = 1
    nord = 1 + 2 * (pld_order - 1)
    self.wax = [GetFrame(76 + 20 * j, 14 + 25 * i, 23, 18) for j in range(self.nseg * 2) for i in range(nord)]
    self.wcax = [GetFrame(76 + 20 * j, 14 + 25 * nord, 4, 18) for j in range(self.nseg * 2)]
    for ax in self.fig.get_axes():
      ax.axis('off')
      
  def weights_grid(self):
    '''
    
    '''
    
    ax = np.array([a() for a in self.wax]).reshape(2 * self.nseg, -1)
    cax = np.array([c() for c in self.wcax])
    
    return ax, cax
    
  def title(self):
    '''
    
    '''
    
    return self.title_left(on = False), self.title_center(on = False), self.title_right(on = False)

  def footer(self):
    '''
    
    '''
    
    return self.footer_left(on = False), self.footer_center(on = False), self.footer_right(on = False)
  
  def top_right(self):
    '''
    
    '''
    
    res = self.body_top_right[self.tcount]()
    self.tcount += 1
    return res
  
  def top_left(self):
    '''
    
    '''
    
    return self.body_top_left()
  
  def left(self):
    '''
    
    '''
    
    res = self.body_left[self.lcount]()
    self.lcount += 1
    return res

  def right(self):
    '''
    
    '''
    
    res = self.body_right[self.rcount]()
    self.rcount += 1
    return res