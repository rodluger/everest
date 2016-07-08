#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
neigbbors.py
------------

'''

from .. import __version__ as EVEREST_VERSION
from ..config import EVEREST_DAT, EVEREST_SRC
from ..data import Campaign, GetK2Campaign
from .ccd import ModuleNumber, Channels
from k2plr.config import KPLR_ROOT
try:
  import pyfits
except ImportError:
  try:
    import astropy.io.fits as pyfits
  except ImportError:
    raise Exception('Please install the `pyfits` package.')
import os, sys
import matplotlib.pyplot as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib.widgets import Slider, Button
from matplotlib.widgets import AxesWidget
import numpy as np

# Color map
coolwarm = pl.cm.coolwarm
coolwarm.set_under(alpha = 0.)
def color(delta_k, hide):
  x = delta_k / 10. + 0.5
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

class Neighbors(object):
    
  def __init__(self, EPIC, maxc = 0.25, maxr = 250, maxm = 5, include_self = False,
               command_line = False):
        
    # Init
    self._maxr = maxr
    self._maxc = maxc
    self._maxm = maxm
    self.EPIC = EPIC
    self.maxc = maxc
    self.maxr = maxr
    self.maxm = maxm
    self.include_self = include_self
    self.cancel = True
    self.command_line = command_line
    self.GetNeighbors()
  
    # Create the figure
    self.fig = pl.figure(figsize = (8,8))
    self.fig.subplots_adjust(top = 0.95, bottom = 0.05, left = 0.05, right = 0.95)
    self.ax = self.fig.add_subplot(111, aspect = 'equal')
    self.ax.set_axis_bgcolor('k')
    self.scatter = self.ax.scatter(self.x, self.y, 
                                   s = 50,
                                   c = color(self.dk, (self.c > self.maxc) | (self.r > self.maxr)),
                                   edgecolor = 'none', zorder = 99,
                                   picker = True)
    self.ax.set_xlim(-self.maxr, self.maxr)
    self.ax.set_ylim(-self.maxr, self.maxr)
    self.ax.xaxis.set_visible(False)
    self.ax.yaxis.set_visible(False)
  
    # Colorbar
    self.divider = make_axes_locatable(self.ax)
    self.axcbar = self.divider.append_axes("top", size="5%", pad=0.1)
    self.axcbar.yaxis.set_visible(False)
    self.axcbar.xaxis.tick_top()
    z = np.atleast_2d(np.linspace(0,1,100))
    self.axcbar.imshow(z, cmap = coolwarm, extent = (-5,5,0,1))
    self.axcbar.set_aspect(1./2.5)
    self.axcbar.set_xticks([-5,-4,-3,-2,-1,0,1,2,3,4,5])
    self.axcbar.set_xticklabels(['-5','-4','-3','-2','-1','0','+1','+2','+3','+4','+5'])
      
    # Selector
    self.circle = pl.Circle((0,0), 0.1 * self.maxr, color = 'lightskyblue', alpha = 0.5, ec = 'none', zorder = 99)
    self.ax.add_artist(self.circle)
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

    # Slider 2 (Radius)
    self.axrad = self.divider.append_axes("bottom", size="5%", pad=0.1)
    self.slrad = Slider(self.axrad, '', 1., 4.25, valinit = np.log10(self.maxr),
                        facecolor = 'lightcoral')
    self.axrad.set_xlim(1.,4.25)
    self.axrad.set_xticks([1., 2., 3., 4.])
    self.axrad.set_xticklabels([10., 100., 1000., 1000.])
    self.slrad.on_changed(self.update_rad)
    self.slrad.valtext.set_visible(False)
    self.axrad.annotate('Radius', xy = (0.5, 0.5), xycoords = 'axes fraction',
                        ha = 'center', va = 'center', fontweight = 'bold', alpha = 0.5)

    # Slider 3 (Delta mag)
    self.axmag = self.divider.append_axes("right", size="5%", pad=0.1)
    self.slmag = VerticalSlider(self.axmag, '', 0., 10., valinit = self.maxm,
                                facecolor = 'lightcoral')
    self.axmag.set_ylim(0,10)
    self.axmag.set_yticks([0, 2.5, 5, 7.5, 10])
    self.axmag.yaxis.tick_right()
    self.slmag.on_changed(self.update_mag)
    self.slmag.valtext.set_visible(False)
    self.axmag.annotate('Delta Mag', xy = (0.5, 0.5), xycoords = 'axes fraction',
                        ha = 'center', va = 'center', fontweight = 'bold', alpha = 0.5,
                        rotation = -90)
    
    # Buttons
    self.axbuttons = self.divider.append_axes("bottom", size="7.5%", pad=0.35)
    self.axbuttons.xaxis.set_visible(False)
    self.axbuttons.yaxis.set_visible(False)
    self.axbuttons.set_frame_on(False)
    self.inset_ok = inset_axes(self.axbuttons, width="25%", height="100%", loc=2)
    self.inset_cancel = inset_axes(self.axbuttons, width="25%", height="100%", loc=9)
    self.inset_reset = inset_axes(self.axbuttons, width="25%", height="100%", loc=1)
    self.ok = Button(self.inset_ok, 'OK')
    self.ok.on_clicked(self.ok_click)
    self.cancel = Button(self.inset_cancel, 'Cancel')
    self.cancel.on_clicked(self.cancel_click)
    self.reset = Button(self.inset_reset, 'Reset')
    self.reset.on_clicked(self.reset_click)
    
    # Labels
    self.label = self.ax.annotate('EPIC %d\nK$_p$ = %.2f' % (self.EPIC, self.kepmag), xy = (0.95, 0.95), ha = 'right', 
                                  va = 'top', xycoords = 'axes fraction',
                                  color = 'w')
    self.label_state = self.ax.annotate('Radius: %d\nContam: %.2f\ndMag: %.2f' % (self.maxr, self.maxc, self.maxm), 
                                        xy = (0.95, 0.05), ha = 'right', 
                                        va = 'bottom', xycoords = 'axes fraction',
                                        color = 'w')
    self.label_number = self.ax.annotate('%d stars' % len(self._targets()), 
                                         xy = (0.05, 0.05), ha = 'left', 
                                         va = 'bottom', xycoords = 'axes fraction',
                                         color = 'w')
                                        
    # Show
    pl.show()

  def onpick(self, event):
    i = event.ind[0]
    self.label.set_text('EPIC %d\nK$_p$ = %.2f' % (self.e[i], self.kepmag + self.dk[i]))
    self.circle.center = (self.x[i], self.y[i])
    self.fig.canvas.draw_idle()
    
  def update_ctm(self, val):
    self.maxc = self.slctm.val
    self.scatter.set_color(color(self.dk, (self.c > self.maxc) | (self.r > self.maxr) | (np.abs(self.dk) > self.maxm)))
    self.label_state.set_text('Radius: %d\nContam: %.2f\ndMag: %.2f' % (self.maxr, self.maxc, self.maxm))
    self.label_number.set_text('%d stars' % len(self._targets()))
    self.fig.canvas.draw()
  
  def update_rad(self, val):
    self.maxr = max(10, 10 ** self.slrad.val)
    self.scatter.set_color(color(self.dk, (self.c > self.maxc) | (self.r > self.maxr) | (np.abs(self.dk) > self.maxm)))
    self.ax.set_xlim(-self.maxr, self.maxr)
    self.ax.set_ylim(-self.maxr, self.maxr)
    self.circle.set_radius(0.1 * self.maxr)
    self.label_state.set_text('Radius: %d\nContam: %.2f\ndMag: %.2f' % (self.maxr, self.maxc, self.maxm))
    self.label_number.set_text('%d stars' % len(self._targets()))
    self.fig.canvas.draw()

  def update_mag(self, val):
    self.maxm = max(0, min(10, self.slmag.val))
    self.scatter.set_color(color(self.dk, (self.c > self.maxc) | (self.r > self.maxr) | (np.abs(self.dk) > self.maxm)))
    self.label_state.set_text('Radius: %d\nContam: %.2f\ndMag: %.2f' % (self.maxr, self.maxc, self.maxm))
    self.label_number.set_text('%d stars' % len(self._targets()))
    self.fig.canvas.draw()
  
  def ok_click(self, event):
    self.cancel = False
    pl.close()
    if self.command_line:
      print(self._targets())
  
  def cancel_click(self, event):
    pl.close()
    sys.exit()
  
  def reset_click(self, event):
    self.slctm.set_val(self._maxc); self.update_ctm(None)
    self.slrad.set_val(np.log10(self._maxr)); self.update_rad(None)
    self.slmag.set_val(self._maxm); self.update_mag(None)
  
  def _targets(self):
    '''
    
    '''
    
    if self.include_self:
      return self.e[(self.c <= self.maxc) & (self.r <= self.maxr) & (np.abs(self.dk) <= self.maxm)]
    else:
      return self.e[(self.c <= self.maxc) & (self.r <= self.maxr) & (np.abs(self.dk) <= self.maxm) & (self.e != self.EPIC)]

  @property
  def targets(self):
    if not self.cancel:
      return self._targets()
    else:
      return None
  
  @property
  def radii(self):
    if not self.cancel:
      return [self.r[np.argmax(self.e == t)] for t in self.targets]
    else:
      return None
      
  @property
  def deltak(self):
    if not self.cancel:
      return [self.dk[np.argmax(self.e == t)] for t in self.targets]
    else:
      return None
        
  def GetNeighbors(self, same_module = False):
    '''
  
    '''
  
    # Get campaign number
    campaign = Campaign(self.EPIC)
    
    # Get info on all neighbors
    info_file = os.path.join(EVEREST_SRC, 'tables', 'C%02d.npz' % campaign)
    if not os.path.exists(info_file):
      info = {}
      for channel in range(1,85):
        info.update({channel: []})
      stars = GetK2Campaign(campaign)
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
          info[channel].append([star, ra, dec, kp, contamination])
      np.savez(info_file, info = info)
    info = np.load(info_file)['info'][()]    
  
    # Get position from FITS file
    file = os.path.join(EVEREST_DAT, 'fits', 'c%02d' % campaign, 
                        ('%09d' % self.EPIC)[:4] + '00000', ('%09d' % self.EPIC)[4:],
                        'hlsp_everest_k2_llc_%d-c%02d_kepler_v%s_lc.fits' % 
                        (self.EPIC, campaign, EVEREST_VERSION))
    channel = pyfits.getval(file, 'CHANNEL', 0)
    ra = pyfits.getval(file, 'RA_OBJ', 0)
    dec = pyfits.getval(file, 'DEC_OBJ', 0)
    crpix1 = pyfits.getval(file, 'CRPIX1', 4)
    crpix2 = pyfits.getval(file, 'CRPIX2', 4)
    crval1 = pyfits.getval(file, 'CRVAL1', 4) 
    crval2 = pyfits.getval(file, 'CRVAL2', 4) 
    cdelt1 = pyfits.getval(file, 'CDELT1', 4)   
    cdelt2 = pyfits.getval(file, 'CDELT1', 4) 
    pc1_1 = pyfits.getval(file, 'PC1_1', 4) 
    pc1_2 = pyfits.getval(file, 'PC1_2', 4)
    pc2_1 = pyfits.getval(file, 'PC2_1', 4)
    pc2_2 = pyfits.getval(file, 'PC2_2', 4)
    pc = np.array([[pc1_1, pc1_2], [pc2_1, pc2_2]])
    pc = np.linalg.inv(pc)
    crpix1p = pyfits.getval(file, 'CRPIX1P', 4)
    crpix2p = pyfits.getval(file, 'CRPIX2P', 4)
    crval1p = pyfits.getval(file, 'CRVAL1P', 4) 
    crval2p = pyfits.getval(file, 'CRVAL2P', 4) 
    cdelt1p = pyfits.getval(file, 'CDELT1P', 4)    
    cdelt2p = pyfits.getval(file, 'CDELT2P', 4) 
  
    if same_module:
      # Get all stars on the same module
      channels = Channels(ModuleNumber(channel))
      info = info[channels[0]] + info[channels[1]] + info[channels[2]] + info[channels[3]]
    else:
      # Get all stars on detector
      all_info = []
      for c in range(1,85):
        all_info.extend(info[c])
      info = all_info
    
    # Get the kepler mag of this target
    kepmag = info[np.argmax(np.array([i[0] == self.EPIC for i in info]))][3]
  
    e = np.zeros(len(info)) * np.nan
    x = np.zeros(len(info)) * np.nan
    y = np.zeros(len(info)) * np.nan
    r = np.zeros(len(info)) * np.nan
    dk = np.zeros(len(info)) * np.nan
    c = np.zeros(len(info)) * np.nan
  
    for i in range(len(info)):
      ra = info[i][1]
      dec = info[i][2]
      kp = info[i][3]
      dra = (ra - crval1) * np.cos(np.radians(dec)) / cdelt1
      ddec = (dec - crval2) / cdelt2
      x[i] = pc[0,0] * dra + pc[0,1] * ddec
      y[i] = pc[1,0] * dra + pc[1,1] * ddec
      r[i] = np.sqrt(x[i] ** 2 + y[i] ** 2)
      dk[i] = kp - kepmag
      c[i] = info[i][4]
      e[i] = info[i][0]
    
    # Make self contam 0 (for plotting)
    c[np.argmax(e == self.EPIC)] = 0
    
    self.e = np.array(e, dtype = int)
    self.x = x
    self.y = y
    self.r = r
    self.dk = dk
    self.c = c
    self.kepmag = kepmag