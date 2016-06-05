#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`ccd.py` - Plotting the CCD image
-----------------------------------------

Classes for constructing and plotting the
`Kepler` detector.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import matplotlib.pyplot as pl
from matplotlib.patches import Rectangle
import numpy as np
dx = 0.05

class Output(object):
  '''
  
  '''
  
  def __init__(self, module, output, channel):
    self.module = module
    self.output = output
    self.channel = channel
    # A list of (column, row) tuples
    self.sources = []       
  
  def __repr__(self):
    return "<Module %d: Output %d (Channel %d)>" % (self.module, self.output, self.channel)
  
  def draw(self, ax, xy):
    ax.add_patch(Rectangle(xy, 0.5 - 2 * dx, 0.5 - 2 * dx, fill = False))
    ax.annotate('%d/%d' % (self.channel, self.output), xy = xy + [0.2, 0.2], 
                va = 'center', ha = 'center', fontsize = 8)
    
    for source in self.sources:
      x = (source[0] - 12) / (1111 - 12)
      y = (source[1] - 20) / (1043 - 20)
      ax.plot(xy[0] + (0.5 - 2 * dx) * x, 
              xy[1] + (0.5 - 2 * dx) * y, 
              'ro', zorder = 99, markersize = 8)

class Module(object):
  '''
  
  '''
  
  def __init__(self, n, c0, orientation):
    if orientation == 0:
      self._outputs = np.array([[Output(n, 4, c0 + 3), Output(n, 3, c0 + 2)],
                                [Output(n, 1, c0), Output(n, 2, c0 + 1)]])
      self._inds = np.array([(1,0), (1,1), (0,1), (0,0)])
    elif orientation == 1:
      self._outputs = np.array([[Output(n, 3, c0 + 2), Output(n, 2, c0 + 1)],
                                [Output(n, 4, c0 + 3), Output(n, 1, c0)]])
      self._inds = np.array([(1,1), (0,1), (0,0), (1,0)])
    elif orientation == 2:
      self._outputs = np.array([[Output(n, 1, c0), Output(n, 4, c0 + 3)],
                                [Output(n, 2, c0 + 1), Output(n, 3, c0 + 2)]])
      self._inds = np.array([(0,0), (1,0), (1,1), (0,1)])  
    elif orientation == 3:
      self._outputs = np.array([[Output(n, 2, c0 + 1), Output(n, 1, c0)],
                                [Output(n, 3, c0 + 2), Output(n, 4, c0 + 3)]])
      self._inds = np.array([(0,1), (0,0), (1,0), (1,1)])
      
    self.n = n
  
  def __getitem__(self, n):
    if n in range(1,5):
      return self._outputs[tuple(self._inds[n - 1])]
    else:
      return None

  def __repr__(self):
    return "<Module %d>" % self.n

  def draw(self, ax, xy):
    ax.add_patch(Rectangle(xy, 1, 1, fill = False))
    ax.annotate('%02d' % self.n, xy = xy + [0.5, 0.5], 
                      va = 'center', ha = 'center',
                      color = 'k', fontsize = 60, alpha = 0.05)
    for n in range(1,5):
      i = self._inds[n - 1]
      offset = i / 2
      for n in range(2):
        if tuple(i)[n] == 0:
          offset[n] += dx + dx/3
        else:
          offset[n] += dx - dx/3
      self._outputs.T[tuple(i)].draw(ax, xy + offset) 

class Empty(Module):
  '''
  
  '''
  
  def __init__(self, n):
    self.n = n
  
  def __getitem__(self, n):
    return None 

  def draw(self, ax, xy):
    pass

class CCD(object):
  '''
  
  '''
  
  def __init__(self):
    self._modules = np.array([[Empty(1), Module(2, 1, 0), Module(3, 5, 0), Module(4, 9, 0), Empty(5)],
                              [Module(6, 13, 1), Module(7, 17, 0), Module(8, 21, 0), Module(9, 25, 2), Module(10, 29, 2)],
                              [Module(11, 33, 1), Module(12, 37, 1), Module(13, 41, 1), Module(14, 45, 2), Module(15, 49, 2)],
                              [Module(16, 53, 1), Module(17, 57, 1), Module(18, 61, 3), Module(19, 65, 3), Module(20, 69, 2)],
                              [Empty(21), Module(22, 73, 3), Module(23, 77, 3), Module(24, 81, 3), Empty(25)]])
    self._inds = np.array([(0,0), (0,1), (0,2), (0,3), (0,4),
                           (1,0), (1,1), (1,2), (1,3), (1,4),
                           (2,0), (2,1), (2,2), (2,3), (2,4),
                           (3,0), (3,1), (3,2), (3,3), (3,4),
                           (4,0), (4,1), (4,2), (4,3), (4,4)])
    self.fig = None
    self.draw()
    
  def __getitem__(self, n):
    if n in range(1,26):
      return self._modules[tuple(self._inds[n - 1])]
    else:
      return None
  
  def channel(self, n):
    if n in range(1, 85):
      modules = [m for m in self._modules.flatten() if type(m) is not Empty]
      i, j = divmod(n, 4)
      if j == 0:
        j = 4
        i -= 1
      return modules[i][j]
    else:
      return None
  
  def draw(self):
    '''
    
    '''
    
    if self.fig is not None:
      pl.close(self.fig)
    self.fig = pl.figure(figsize = (7,7))
    self.fig.subplots_adjust(top = 0.97, bottom = 0.03, left = 0.03, right = 0.97)
    self.ax = self.fig.add_subplot(111, aspect = 'equal')
    self.ax.xaxis.set_visible(False)
    self.ax.yaxis.set_visible(False)
    self.ax.set_xlim(0-dx,5+dx)
    self.ax.set_ylim(5+dx,0-dx)
    
    for n in range(1,26):
      i = self._inds[n - 1]
      self._modules.T[tuple(i)].draw(self.ax, i)
  
  def add_source(self, channel, crval1p, crval2p):
    '''
    
    '''
    
    self.channel(channel).sources.append((crval1p, crval2p))
    self.draw()
    