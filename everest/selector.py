#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
selector.py
-----------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import matplotlib
import matplotlib.pyplot as pl
from matplotlib.widgets import RectangleSelector, Button
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
    
    # Initialize our widgets
    self.Selector = RectangleSelector(ax, self.on_select,
                                      rectprops = dict(facecolor='red', 
                                                       edgecolor = 'black',
                                                       alpha=0.1, fill=True),
                                      useblit = True,
                                      minspanx = 0,
                                      minspany = 0)
    self.Selector.set_active(False)
                                   
    self.Unselector = RectangleSelector(ax, self.on_unselect,
                                      rectprops = dict(facecolor='red', 
                                                       edgecolor = 'black',
                                                       alpha=0.1, fill=True),
                                      useblit = True,
                                      minspanx = 0,
                                      minspany = 0)
    self.Unselector.set_active(False)
    
    axsel = pl.axes([0.86, 0.12, 0.08, 0.04])
    self.SelectButton = Button(axsel, 'Select')
    self.SelectButton.on_clicked(self.on_select_button)
    
    axunsel = pl.axes([0.77, 0.12, 0.08, 0.04])
    self.UnselectButton = Button(axunsel, 'Unselect')
    self.UnselectButton.on_clicked(self.on_unselect_button)
    
    axrec = pl.axes([0.68, 0.12, 0.08, 0.04])
    self.RecalcButton = Button(axrec, 'Recalc')
    self.RecalcButton.on_clicked(self.on_recalc_button)
                
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

  def on_select(self, eclick, erelease):
    '''
    
    '''
    
    inds = self.get_inds(eclick, erelease)
    for i in inds:
      if i not in self._selected:
        self._selected.append(i)   
    self.redraw()

  def on_unselect(self, eclick, erelease):
    '''
    
    '''
    
    inds = self.get_inds(eclick, erelease)
    for i in inds:
      if i in self._selected:
        self._selected.remove(i) 
    self.redraw()
  
  def on_select_button(self, event):
    '''
    
    '''
    
    self.Unselector.set_active(False)
    self.UnselectButton.label.set_weight('normal')
    
    self.Selector.set_active(not self.Selector.active)
    if self.Selector.active:
      self.SelectButton.label.set_weight('bold')
    else:
      self.SelectButton.label.set_weight('normal')
    self.redraw()
  
  def on_unselect_button(self, event):
    '''
    
    '''
    
    self.Selector.set_active(False)
    self.SelectButton.label.set_weight('normal')
    self.Unselector.set_active(not self.Unselector.active)
    if self.Unselector.active:
      self.UnselectButton.label.set_weight('bold')
    else:
      self.UnselectButton.label.set_weight('normal')
    self.redraw()
  
  def on_recalc_button(self, event):
    '''
    
    '''
    
    self.Selector.set_active(False)
    self.Unselector.set_active(False)
    self.UnselectButton.label.set_weight('normal')
    self.SelectButton.label.set_weight('normal')
    self.x, self.y = self.fxy(self.selected)
    self.redraw()
  
  @property
  def selected(self):
    return np.array(sorted(self._selected), dtype = int)