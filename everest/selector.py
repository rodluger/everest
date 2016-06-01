#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
selector.py
-----------

   .. code-block:: python
      fig, ax = pl.subplots(1)  
      pl.plot([0,100],[0,100], 'r-')
      sel = Selector(fig, ax)
      print(sel.transits)

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import matplotlib
import matplotlib.pyplot as pl
from matplotlib.widgets import RectangleSelector
import numpy as np

# Make this a global variable
rcParams = None

# Keyboard shortcuts
TRAN = 't'
SUBT = 'alt'

def DisableShortcuts():
  '''
  Disable MPL keyboard shortcuts and the plot toolbar.

  '''
  
  global rcParams
  rcParams = dict(pl.rcParams)
  pl.rcParams['keymap.all_axes'] = ''
  pl.rcParams['keymap.back'] = ''
  pl.rcParams['keymap.forward'] = ''
  pl.rcParams['keymap.fullscreen'] = ''
  pl.rcParams['keymap.grid'] = ''
  pl.rcParams['keymap.home'] = ''
  pl.rcParams['keymap.pan'] = ''
  pl.rcParams['keymap.quit'] = ''
  pl.rcParams['keymap.save'] = ''
  pl.rcParams['keymap.xscale'] = ''
  pl.rcParams['keymap.yscale'] = ''
  pl.rcParams['keymap.zoom'] = ''

def EnableShortcuts():
  '''
  Resets pl.rcParams to its original state.
  
  '''
  
  global rcParams
  pl.rcParams.update(rcParams)

class Selector(object):
  '''
  
  '''
  
  def __init__(self, fig, ax, x = np.arange(0, 100), y = np.random.randn(100)):
  
    # Initialize
    DisableShortcuts()
    self.fig = fig
    self.ax = ax    
    self.x = x
    self.y = y
        
    # Initialize arrays
    self._transits = []
    self._plots = []
    self.info = ""
    self.subt = False
    
    # Initialize our selectors
    self.Transits = RectangleSelector(ax, self.oselect,
                                      rectprops = dict(facecolor='red', 
                                                       edgecolor = 'black',
                                                       alpha=0.1, fill=True),
                                      useblit = True,
                                      minspanx = 0,
                                      minspany = 0)
    self.Transits.set_active(False)
  
    pl.connect('key_press_event', self.on_key_press)
    pl.connect('key_release_event', self.on_key_release)
                                         
    self.redraw()
    
  def redraw(self):
  
    # Reset the figure
    for p in self._plots:
      if len(p): p.pop(0).remove()
    self._plots = []
    
    # Non-transit inds
    inds = [i for i, _ in enumerate(self.x) if (i not in self.transits)]
    p = self.ax.plot(self.x[inds], self.y[inds], 'b.', markersize = 3, alpha = 0.5)
    self._plots.append(p)
    color = p[0].get_color()
    
    # Transits
    t = self.transits
    p = self.ax.plot(self.x[t], self.y[t], 'r.', markersize = 3, alpha = 0.5)
    self._plots.append(p)
          
    # Labels
    self.ax.set_xlabel('Time', fontsize = 22)
    self.ax.set_ylabel('Flux', fontsize = 22)
    
    # Operations
    label = None
    if self.Transits.active:
      if self.subt:
        label = 'transits (-)'
      else:
        label = 'transits (+)'
    
    if label is not None:
      a = self.ax.text(0.005, 0.96, label, fontsize=12, transform=self.ax.transAxes, color = 'r', alpha = 0.75)
      self._plots.append([a])
      
    # Refresh
    self.fig.canvas.draw()
  
  def get_inds(self, eclick, erelease):
    '''
    
    '''
        
    # Get coordinates
    x1 = min(eclick.xdata, erelease.xdata)
    x2 = max(eclick.xdata, erelease.xdata)
    y1 = min(eclick.ydata, erelease.ydata)
    y2 = max(eclick.ydata, erelease.ydata)

    if (x1 == x2) and (y1 == y2):   
      # User clicked 
      d = ((self.x - x1)/(max(self.x) - min(self.x))) ** 2 + ((self.y - y1)/(max(self.y) - min(self.y))) ** 2
      return [np.argmin(d)]
        
    else:
      # User selected
      xy = zip(self.x, self.y)
      return [i for i, pt in enumerate(xy) if x1<=pt[0]<=x2 and y1<=pt[1]<=y2] 

  def oselect(self, eclick, erelease):
    '''
    
    '''
    
    inds = self.get_inds(eclick, erelease)
    for i in inds:
      if self.subt:
        if i in self._transits:
          self._transits.remove(i)
      else:
        if i not in self._transits:
          self._transits.append(i)   
    self.redraw()
  
  def on_key_release(self, event):
    
    # Alt
    if event.key == SUBT:
      self.subt = False
      self.redraw()
      
  def on_key_press(self, event):
    
    # Alt
    if event.key == SUBT:
      self.subt = True
      self.redraw()
    
    # Transits
    elif event.key == TRAN:
      self.Transits.set_active(not self.Transits.active)
      self.redraw()

  @property
  def transits(self):
    return np.array(sorted(self._transits), dtype = int)