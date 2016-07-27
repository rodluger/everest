#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`detector.py` - Plotting the Kepler Field
-------------------------------------------------

Tools for interactively plotting the `Kepler` field.

'''

from __future__ import print_function
from .. import __version__ as EVEREST_VERSION
from ..config import EVEREST_DAT, EVEREST_SRC
from ..data import Campaign, GetK2Campaign
from k2plr.config import KPLR_ROOT
try:
  import pyfits
except ImportError:
  try:
    import astropy.io.fits as pyfits
  except ImportError:
    raise Exception('Please install the `pyfits` package.')
import os, sys, subprocess
import matplotlib.pyplot as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib.widgets import Slider, Button
from matplotlib.widgets import AxesWidget
import numpy as np

# Color map
coolwarm = pl.cm.coolwarm
coolwarm.set_under(alpha = 0.)
def color(kp, hide):
  x = (kp - 13) / 5. + 0.5
  x[x<0] = 0; x[x>1] = 1
  x[hide] = -1
  return coolwarm(x)

class VerticalSlider(AxesWidget):
    """
    A vertical slider widget. Adapted from
    `this page <http://stackoverflow.com/questions/25934279/add-a-vertical-slider-with-matplotlib>`_.

    The following attributes are defined
      *ax*        : the slider :class:`matplotlib.axes.Axes` instance

      *val*       : the current slider value

      *vline*     : a :class:`matplotlib.lines.Line2D` instance
                     representing the initial value of the slider

      *poly*      : A :class:`matplotlib.patches.Polygon` instance
                     which is the slider knob

      *valfmt*    : the format string for formatting the slider text

      *label*     : a :class:`matplotlib.text.Text` instance
                     for the slider label

      *closedmin* : whether the slider is closed on the minimum

      *closedmax* : whether the slider is closed on the maximum

      *slidermin* : another slider - if not *None*, this slider must be
                     greater than *slidermin*

      *slidermax* : another slider - if not *None*, this slider must be
                     less than *slidermax*

      *dragging*  : allow for mouse dragging on slider

    Call :meth:`on_changed` to connect to the slider event
    """
    def __init__(self, ax, label, valmin, valmax, valinit=0.5, valfmt='%1.2f',
                 closedmin=True, closedmax=True, slidermin=None,
                 slidermax=None, dragging=True, **kwargs):
        """
        Create a slider from *valmin* to *valmax* in axes *ax*

        *valinit*
            The slider initial position

        *label*
            The slider label

        *valfmt*
            Used to format the slider value

        *closedmin* and *closedmax*
            Indicate whether the slider interval is closed

        *slidermin* and *slidermax*
            Used to constrain the value of this slider to the values
            of other sliders.

        additional kwargs are passed on to ``self.poly`` which is the
        :class:`matplotlib.patches.Rectangle` which draws the slider
        knob.  See the :class:`matplotlib.patches.Rectangle` documentation
        valid property names (e.g., *facecolor*, *edgecolor*, *alpha*, ...)
        """
        AxesWidget.__init__(self, ax)

        self.valmin = valmin
        self.valmax = valmax
        self.val = valinit
        self.valinit = valinit
        self.poly = ax.axhspan(valmin, valinit, 0, 1, **kwargs)

        self.vline = ax.axhline(valinit, 0, 1, color='r', lw=1)

        self.valfmt = valfmt
        ax.set_xticks([])
        ax.set_ylim((valmin, valmax))
        ax.set_yticks([])
        ax.set_navigate(False)

        self.connect_event('button_press_event', self._update)
        self.connect_event('button_release_event', self._update)
        if dragging:
            self.connect_event('motion_notify_event', self._update)
        self.label = ax.text(0.5, 1.03, label, transform=ax.transAxes,
                             verticalalignment='center',
                             horizontalalignment='center')

        self.valtext = ax.text(0.5, -0.03, valfmt % valinit,
                               transform=ax.transAxes,
                               verticalalignment='center',
                               horizontalalignment='center')

        self.cnt = 0
        self.observers = {}

        self.closedmin = closedmin
        self.closedmax = closedmax
        self.slidermin = slidermin
        self.slidermax = slidermax
        self.drag_active = False

    def _update(self, event):
        """update the slider position"""
        if self.ignore(event):
            return

        if event.button != 1:
            return

        if event.name == 'button_press_event' and event.inaxes == self.ax:
            self.drag_active = True
            event.canvas.grab_mouse(self.ax)

        if not self.drag_active:
            return

        elif ((event.name == 'button_release_event') or
              (event.name == 'button_press_event' and
               event.inaxes != self.ax)):
            self.drag_active = False
            event.canvas.release_mouse(self.ax)
            return

        val = event.ydata
        if val <= self.valmin:
            if not self.closedmin:
                return
            val = self.valmin
        elif val >= self.valmax:
            if not self.closedmax:
                return
            val = self.valmax

        if self.slidermin is not None and val <= self.slidermin.val:
            if not self.closedmin:
                return
            val = self.slidermin.val

        if self.slidermax is not None and val >= self.slidermax.val:
            if not self.closedmax:
                return
            val = self.slidermax.val

        self.set_val(val)

    def set_val(self, val):
        xy = self.poly.xy
        xy[1] = 0, val
        xy[2] = 1, val
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % val)
        if self.drawon:
            self.ax.figure.canvas.draw()
        self.val = val
        if not self.eventson:
            return
        for cid, func in self.observers.items():
            func(val)

    def on_changed(self, func):
        """
        When the slider value is changed, call *func* with the new
        slider position

        A connection id is returned which can be used to disconnect
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        """remove the observer with connection id *cid*"""
        try:
            del self.observers[cid]
        except KeyError:
            pass

    def reset(self):
        """reset the slider to the initial value if needed"""
        if (self.val != self.valinit):
            self.set_val(self.valinit)

class Detector(object):
  '''
  Plots all stars in a given campaign as they appear in the `Kepler` field
  of view. Stars are colored according to their `Kp` magnitude, and a slider
  at the left allows the user to make cuts based on the contamination metric.
  Double-clicking on a star will download and plot its :py:mod:`everest`
  light curve.
  
  '''
  
  def __init__(self, campaign, epic = None, maxc = 0.25):
        
    # Init
    self._maxc = maxc
    self.maxc = maxc
    self.cancel = True
    self.campaign = campaign
    self.GetStars()
  
    # Create the figure
    self.fig = pl.figure(figsize = (8,8))
    self.fig.subplots_adjust(top = 0.95, bottom = 0.05, left = 0.05, right = 0.95)
    self.ax = self.fig.add_subplot(111)
    self.ax.set_axis_bgcolor('k')
    self.scatter = self.ax.scatter(self.x, self.y, 
                                   s = 50,
                                   c = color(self.kp, (self.c > self.maxc)),
                                   edgecolor = 'none', zorder = 99,
                                   picker = True)
    self.ax.margins(0.1,0.1)
    self.ax.xaxis.set_visible(False)
    self.ax.yaxis.set_visible(False)
    self.fig.canvas.set_window_title('Campaign %d' % self.campaign)
    
    # Colorbar
    self.divider = make_axes_locatable(self.ax)
    self.axcbar = self.divider.append_axes("top", size="5%", pad=0.1)
    self.axcbar.yaxis.set_visible(False)
    self.axcbar.xaxis.tick_top()
    z = np.atleast_2d(np.linspace(0,1,100))
    self.axcbar.imshow(z, cmap = coolwarm, extent = (8,18,0,1))
    self.axcbar.set_aspect(1./2.5)
    self.axcbar.annotate('Kepler Magnitude', xy = (0.5, 0.5), xycoords = 'axes fraction',
                          ha = 'center', va = 'center', fontweight = 'bold', alpha = 0.5)
      
    # Selector
    self.circle = None
    self.fig.canvas.mpl_connect('pick_event', self.onpick)
    
    # Slider 1 (Contamination)
    self.axctm = self.divider.append_axes("left", size="5%", pad=0.1)
    self.slctm = VerticalSlider(self.axctm, '', 0., 1., valinit = self.maxc,
                                facecolor = 'lightcoral')
    self.axctm.set_ylim(0,1)
    self.axctm.set_yticks(np.linspace(0,1,11))
    self.slctm.on_changed(self.update_ctm)
    self.slctm.valtext.set_visible(False)
    self.axctm.annotate('Contamination', xy = (0.5, 0.5), xycoords = 'axes fraction',
                        ha = 'center', va = 'center', fontweight = 'bold', alpha = 0.5,
                        rotation = 90)
    
    # Labels
    self.label = self.ax.annotate('', xy = (0.95, 0.95), ha = 'right', 
                                  va = 'top', xycoords = 'axes fraction',
                                  color = 'w', zorder = 99)
    
    # Initial EPIC?
    if epic is not None:
      i = np.argmax(self.e == epic)
      self.label.set_text('EPIC %d\nK$_p$ = %.2f' % (self.e[i], self.kp[i]))
      self.circle = self.ax.scatter(self.x[i],self.y[i],s=500,c='yellow',alpha=0.5,zorder=99)
      self.fig.canvas.draw_idle()
                                      
    # Show
    pl.show()

  def onpick(self, event):  
    select = False
    for i in event.ind:
      if (self.c[i] > self.maxc):
        continue
      else:
        self.label.set_text('EPIC %d\nK$_p$ = %.2f' % (self.e[i], self.kp[i]))
        if self.circle is not None:
          self.circle.remove()
        self.circle = self.ax.scatter(self.x[i],self.y[i],s=500,c='yellow',alpha=0.5,zorder=99)
        self.fig.canvas.draw_idle()
        select = True
        break
    if select and event.mouseevent.dblclick:
      print("Downloading EPIC %d..." % self.e[i])
      subprocess.Popen(['everest', '%d' % self.e[i]])
    
  def update_ctm(self, val):
    self.maxc = self.slctm.val
    self.scatter.set_color(color(self.kp, (self.c > self.maxc)))
    self.fig.canvas.draw()
        
  def GetStars(self):
    '''
    Get all stars on the detector for this campaign.
    
    '''
      
    # Get info on all neighbors
    info_file = os.path.join(EVEREST_SRC, 'tables', 'C%02d.info' % self.campaign)
    if not os.path.exists(info_file):
      with open(info_file, 'w') as f:
        stars = GetK2Campaign(self.campaign)
        nstars = len(stars)
        for i, star in enumerate(stars):
          print("Running EPIC %d (%d/%d)..." % (star, i + 1, nstars))
          file = os.path.join(KPLR_ROOT, 'data', 'everest', str(star), str(star) + '.npz')
          if os.path.exists(file):
            data = np.load(file)
            ra = data['fitsheader'][0]['RA_OBJ'][1]
            dec = data['fitsheader'][0]['DEC_OBJ'][1]
            foo = [s['kepmag'] for s in data['nearby'] if s['epic'] == star]
            if len(foo):
              kp = foo[0]
            else:
              kp = np.nan
            contamination = data['contamination'][()]
            channel = data['fitsheader'][0]['CHANNEL'][1]
            print("{:>09d} {:>03d} {:>15.10f} {:>15.10f} {:>6.3f} {:>6.3f}".format(star, channel, ra, dec, kp, contamination), file = f)

    self.e, channel, self.x, self.y, self.kp, self.c  = np.loadtxt(info_file, unpack = True)  
    self.e = np.array(self.e, dtype = int)