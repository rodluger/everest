#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`selector.py` - Interactive plotting
--------------------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import matplotlib
import matplotlib.pyplot as pl
from matplotlib.widgets import RectangleSelector, Button, Slider
import numpy as np

def ToggleAxis(ax):
  '''
  
  '''
  
  pos = ax.get_position()
  pos = [-pos.x0, -pos.y0, pos.width, pos.height]
  ax.set_position(pos)
  ax.set_visible(not ax.get_visible())

def HideAxis(ax):
  '''
  
  '''
  
  pos = ax.get_position()
  pos = [-np.abs(pos.x0), -np.abs(pos.y0), pos.width, pos.height]
  ax.set_position(pos)
  ax.set_visible(False)

def ShowAxis(ax):
  '''
  
  '''
  
  pos = ax.get_position()
  pos = [np.abs(pos.x0), np.abs(pos.y0), pos.width, pos.height]
  ax.set_position(pos)
  ax.set_visible(True)

class TransitSelector(object):
  '''
  
  '''
  
  def __init__(self, ti, tf):
    self.ti = ti
    self.tf = tf
    self.set_active(False)

  def set_active(self, active):
    self.active = active
    self._t0 = None
    self._t1 = None
    self.per = None
    self.ttimes = []

  @property
  def t0(self):
    return self._t0
  
  @t0.setter
  def t0(self, v):
    self._t0 = v
    self.ttimes = [v]

  @property
  def t1(self):
    return self._t1
  
  @t1.setter
  def t1(self, v):
    self._t1 = v
    self.per = np.abs(self._t1 - self._t0) 
    self.ttimes = np.arange(self._t0 + 
                            np.ceil((self.ti - 1. - self._t0) / self.per) * self.per, 
                            self.tf + 1., self.per)
    
class Selector(object):
  '''
  
  '''
  
  def __init__(self, fig, ax, fxy, selected = [], 
               key_sel = 't', key_sub = 'alt', 
               key_fxy = 'r'):
  
    # Initialize
    self.fig = fig
    self.ax = ax    
    self.fxy = fxy
    self.x, self.y = fxy(selected)
    self.key_sel = key_sel
    self.key_sub = key_sub
    self.key_fxy = key_fxy
        
    # Initialize arrays
    self._original_selected = list(selected)
    self._selected = list(selected)
    self._plots = []
    self.info = ""
    self.subt = False
    
    # Connect the mouse
    pl.connect('button_press_event', self.on_mouse_click) 
    
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
                                      rectprops = dict(facecolor='blue', 
                                                       edgecolor = 'black',
                                                       alpha=0.1, fill=True),
                                      useblit = True,
                                      minspanx = 0,
                                      minspany = 0)
    self.Unselector.set_active(False)
    
    self.TransitSelector = TransitSelector(self.x[0], self.x[-1])
    self.TransitSelector.set_active(False)
    
    axtrn = pl.axes([0.86, 0.12, 0.08, 0.04])
    self.TransitButton = Button(axtrn, 'Transit')
    self.TransitButton.on_clicked(self.on_transit_button)
    
    axsel = pl.axes([0.77, 0.12, 0.08, 0.04])
    self.SelectButton = Button(axsel, 'Select')
    self.SelectButton.on_clicked(self.on_select_button)
    
    axunsel = pl.axes([0.68, 0.12, 0.08, 0.04])
    self.UnselectButton = Button(axunsel, 'Unselect')
    self.UnselectButton.on_clicked(self.on_unselect_button)
    
    axrec = pl.axes([0.59, 0.12, 0.08, 0.04])
    self.RecalcButton = Button(axrec, 'Detrend')
    self.RecalcButton.on_clicked(self.on_recalc_button)
    
    axres = pl.axes([0.50, 0.12, 0.08, 0.04])
    self.ResetButton = Button(axres, 'Reset')
    self.ResetButton.on_clicked(self.on_reset_button)
    
    axcancel = pl.axes([0.86, 0.17, 0.08, 0.04])
    self.CancelButton = Button(axcancel, 'Cancel')
    self.CancelButton.on_clicked(self.on_cancel_button)
    HideAxis(axcancel)
    
    axok = pl.axes([0.77, 0.17, 0.08, 0.04])
    self.OKButton = Button(axok, 'OK')
    self.OKButton.on_clicked(self.on_ok_button)
    HideAxis(axok)
    
    axslider = pl.axes([0.86, 0.22, 0.08, 0.02])
    self.Slider = Slider(axslider, '', 0., 1., valinit = 0.25, facecolor = 'gray')
    def update(val):
      self.redraw()
    self.Slider.on_changed(update)
    self.Slider.valtext.set_visible(False)
    axslider.set_xlim(0,1)
    axslider.set_xticks([0,1])
    axslider.set_xticklabels(["0.","1."], fontsize = 8)
    axslider.xaxis.tick_top()
    axslider.set_title('Duration (d)', fontsize = 8, y = 1.1)
    HideAxis(axslider)
            
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
    
    # Transits
    if self.TransitSelector.active:
      for t in self.TransitSelector.ttimes:
        p = [self.ax.axvline(t, color = 'r', ls = '-', alpha = 1.0)]
        self._plots.append(p)
        p = [self.ax.axvspan(t - self.Slider.val / 2., t + self.Slider.val / 2.,
                            color = 'r', alpha = 0.25)]
        self._plots.append(p)
                 
    # Refresh
    self.fig.canvas.draw()
  
  def on_mouse_click(self, event):
    '''
    
    '''
    
    if (event.inaxes is not None) and (self.TransitSelector.active):
      s = np.argmin(np.abs(self.x - event.xdata))
            
      if self.TransitSelector.t0 is None:
        self.TransitSelector.t0 = self.x[s]
      elif self.TransitSelector.t1 is None:
        self.TransitSelector.t1 = self.x[s]
        ShowAxis(self.OKButton.ax)
        ShowAxis(self.CancelButton.ax)
        ShowAxis(self.Slider.ax)
      self.redraw()
  
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
  
  def on_ok_button(self, event):
    '''
    
    '''
    if self.OKButton.ax.get_visible():
      for t in self.TransitSelector.ttimes:
        self._selected.extend(np.where(np.abs(self.x - t) < self.Slider.val / 2.)[0])
      self.TransitSelector.set_active(False)
      self.TransitButton.label.set_weight('normal')
      HideAxis(self.OKButton.ax)
      HideAxis(self.CancelButton.ax)
      HideAxis(self.Slider.ax)
      self.redraw()
  
  def on_cancel_button(self, event):
    '''
    
    '''
    if self.CancelButton.ax.get_visible():
      self.TransitSelector.set_active(False)
      self.TransitButton.label.set_weight('normal')
      HideAxis(self.OKButton.ax)
      HideAxis(self.CancelButton.ax)
      HideAxis(self.Slider.ax)
      self.redraw()
  
  def on_transit_button(self, event):
    '''
    
    '''
    
    self.Selector.set_active(False)
    self.Unselector.set_active(False)
    HideAxis(self.OKButton.ax)
    HideAxis(self.CancelButton.ax)
    HideAxis(self.Slider.ax)
    self.SelectButton.label.set_weight('normal')
    self.UnselectButton.label.set_weight('normal')
    self.TransitSelector.set_active(not self.TransitSelector.active)
    if self.TransitSelector.active:
      self.TransitButton.label.set_weight('bold')
    else:
      self.TransitButton.label.set_weight('normal')
    self.redraw()
    
  
  def on_select_button(self, event):
    '''
    
    '''
    
    self.Unselector.set_active(False)
    self.TransitSelector.set_active(False)
    HideAxis(self.OKButton.ax)
    HideAxis(self.CancelButton.ax)
    HideAxis(self.Slider.ax)
    self.UnselectButton.label.set_weight('normal')
    self.TransitButton.label.set_weight('normal')
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
    self.TransitSelector.set_active(False)
    HideAxis(self.OKButton.ax)
    HideAxis(self.CancelButton.ax)
    HideAxis(self.Slider.ax)
    self.SelectButton.label.set_weight('normal')
    self.TransitButton.label.set_weight('normal')
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
    self.TransitSelector.set_active(False)
    HideAxis(self.OKButton.ax)
    HideAxis(self.CancelButton.ax)
    HideAxis(self.Slider.ax)
    self.UnselectButton.label.set_weight('normal')
    self.SelectButton.label.set_weight('normal')
    self.TransitButton.label.set_weight('normal')
    self.x, self.y = self.fxy(self.selected)
    self.redraw()

  def on_reset_button(self, event):
    '''
    
    '''
    
    self._selected = list(self._original_selected)
    self.Selector.set_active(False)
    self.Unselector.set_active(False)
    self.TransitSelector.set_active(False)
    HideAxis(self.OKButton.ax)
    HideAxis(self.CancelButton.ax)
    HideAxis(self.Slider.ax)
    self.UnselectButton.label.set_weight('normal')
    self.SelectButton.label.set_weight('normal')
    self.TransitButton.label.set_weight('normal')
    self.x, self.y = self.fxy(self._selected)
    self.TransitSelector = TransitSelector(self.x[0], self.x[-1])
    self.redraw()
  
  @property
  def selected(self):
    return np.array(sorted(self._selected), dtype = int)