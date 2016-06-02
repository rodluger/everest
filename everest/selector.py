#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
selector.py
-----------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import matplotlib
import matplotlib.pyplot as pl
from matplotlib.widgets import RectangleSelector
import numpy as np

# Make this a global variable
rcParams = None

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
  
  def __init__(self, fig, ax, fxy, selected = [], 
               key_sel = 't', key_sub = 'alt', 
               key_fxy = 'r'):
  
    # Initialize
    DisableShortcuts()
    self.fig = fig
    self.ax = ax    
    self.fxy = fxy
    self.x, self.y = fxy(selected)
    self.key_sel = key_sel
    self.key_sub = key_sub
    self.key_fxy = key_fxy
        
    # Initialize arrays
    self._selected = selected
    self._plots = []
    self.info = ""
    self.subt = False
    
    # Initialize our selectors
    self.Selector = RectangleSelector(ax, self.oselect,
                                      rectprops = dict(facecolor='red', 
                                                       edgecolor = 'black',
                                                       alpha=0.1, fill=True),
                                      useblit = True,
                                      minspanx = 0,
                                      minspany = 0)
    self.Selector.set_active(False)
  
    pl.connect('key_press_event', self.on_key_press)
    pl.connect('key_release_event', self.on_key_release)
                                         
    self.redraw()
    
  def redraw(self):
  
    # Reset the figure
    for p in self._plots:
      if len(p): p.pop(0).remove()
    self._plots = []
    
    # Non-selected points
    inds = [i for i, _ in enumerate(self.x) if (i not in self.selected)]
    p = self.ax.plot(self.x[inds], self.y[inds], 'b.', markersize = 3, alpha = 0.5)
    self._plots.append(p)
    color = p[0].get_color()
    
    # Selected points
    t = self.selected
    p = self.ax.plot(self.x[t], self.y[t], 'r.', markersize = 3, alpha = 0.5)
    self._plots.append(p)
           
    # Labels
    label = None
    if self.Selector.active:
      if self.subt:
        label = '-'
      else:
        label = '+'
    
    if label is not None:
      a = self.ax.text(0.975, 0.025, label, fontsize=42, transform=self.ax.transAxes, 
                       va = 'bottom', ha = 'right',
                       color = 'r', alpha = 0.25, fontweight = 'bold',
                       fontname = 'monospace')
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
        if i in self._selected:
          self._selected.remove(i)
      else:
        if i not in self._selected:
          self._selected.append(i)   
    self.redraw()
  
  def on_key_release(self, event):
    
    # Alt
    if event.key == self.key_sub:
      self.subt = False
      self.redraw()
      
  def on_key_press(self, event):
    
    # Deselect
    if event.key == self.key_sub:
      self.subt = True
      self.redraw()
    
    # Select
    elif event.key == self.key_sel:
      self.Selector.set_active(not self.Selector.active)
      self.redraw()
    
    # Recalculate
    elif event.key == self.key_fxy:
      self.Selector.set_active(False)
      self.x, self.y = self.fxy(self.selected)
      self.redraw()
      
  @property
  def selected(self):
    return np.array(sorted(self._selected), dtype = int)