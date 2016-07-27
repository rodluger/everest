#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`utils.py` - Plotting routines
--------------------------------------

Plots a simple :py:mod:`pysyzygy` light curve, along with the
observer view of the orbit, the limb darkening profile, and a 
table of the transit parameters.

.. figure:: ../pysyzygy/img/transit_plot.png
    :width: 600px
    :align: center
    :height: 100px
    :alt: alternate text
    :figclass: align-center

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .transit import Transit, QUADRATIC, KIPPING, NONLINEAR
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap, colorConverter
import sys

__all__ = ['PlotTransit']

def I(r, limbdark):
  '''
  The standard quadratic limb darkening law.
  
  :param ndarray r: The radius vector
  :param limbdark: A :py:class:`pysyzygy.transit.LIMBDARK` instance containing 
                   the limb darkening law information
  
  :returns: The stellar intensity as a function of `r`
  
  '''
  
  if limbdark.ldmodel == QUADRATIC:
    u1 = limbdark.u1
    u2 = limbdark.u2
    return (1-u1*(1-np.sqrt(1-r**2))-u2*(1-np.sqrt(1-r**2))**2)/(1-u1/3-u2/6)/np.pi
  elif limbdark.ldmodel == KIPPING:
    a = np.sqrt(limbdark.q1)
    b = 2*limbdark.q2
    u1 = a*b
    u2 = a*(1 - b)
    return (1-u1*(1-np.sqrt(1-r**2))-u2*(1-np.sqrt(1-r**2))**2)/(1-u1/3-u2/6)/np.pi
  elif limbdark.ldmodel == NONLINEAR:
    raise Exception('Nonlinear model not yet implemented!')                           # TODO!
  else:
    raise Exception('Invalid limb darkening model.')

def PlotTransit(compact = False, ldplot = True, plottitle = "", 
                xlim = None, binned = True, **kwargs):
  '''
  Plots a light curve described by `kwargs`
  
  :param bool compact: Display the compact version of the plot? Default `False`
  :param bool ldplot: Displat the limb darkening inset? Default `True`
  :param str plottitle: The title of the plot. Default `""`
  :param float xlim: The half-width of the orbit plot in stellar radii. Default is to \
                     auto adjust this
  :param bool binned: Bin the light curve model to the exposure time? Default `True`
  :param kwargs: Any keyword arguments to be passed to :py:func:`pysyzygy.transit.Transit`
  
  :returns fig: The :py:mod:`matplotlib` figure object
  
  '''

  # Plotting
  fig = pl.figure(figsize = (12,8))
  fig.subplots_adjust(hspace=0.3)
  ax1, ax2 = pl.subplot(211), pl.subplot(212)
  if not compact:
    fig.subplots_adjust(right = 0.7)
    
  t0 = kwargs.pop('t0', 0.)
  trn = Transit(**kwargs)
  try:
    trn.Compute()
    notransit = False
  except Exception as e:
    if str(e) == "Object does not transit the star.":
      notransit = True
    else: raise Exception(e)

  time = trn.arrays.time + t0
  
  if not notransit:
    if binned:
      trn.Bin()
      flux = trn.arrays.bflx
    else:
      flux = trn.arrays.flux

    time = np.concatenate(([-1.e5], time, [1.e5]))                                    # Add baseline on each side
    flux = np.concatenate(([1.], flux, [1.]))
    ax1.plot(time, flux, '-', color='DarkBlue')
    rng = np.max(flux) - np.min(flux)
  
    if rng > 0:
      ax1.set_ylim(np.min(flux) - 0.1*rng, np.max(flux) + 0.1*rng)
      left = np.argmax(flux < (1. - 1.e-8))
      right = np.argmax(flux[left:] > (1. - 1.e-8)) + left
      rng = time[right] - time[left]        
      ax1.set_xlim(time[left] - rng, time[right] + rng)
  
  ax1.set_xlabel('Time (Days)', fontweight='bold')
  ax1.set_ylabel('Normalized Flux', fontweight='bold')

  # Adjust these for full-orbit plotting
  maxpts = kwargs.get('maxpts', 10000); kwargs.update({'maxpts': maxpts})
  per = kwargs.get('per', 10.); kwargs.update({'per': per})
  kwargs.update({'fullorbit': True})
  kwargs.update({'exppts': 30})
  kwargs.update({'exptime': 50 * per / maxpts})
  trn = Transit(**kwargs)
  
  try:
    trn.Compute()
  except Exception as e:
    if str(e) == "Object does not transit the star.":
      pass
    else: raise Exception(e)

  # Sky-projected motion
  x = trn.arrays.x
  y = trn.arrays.y
  z = trn.arrays.z
  inc = (np.arccos(trn.transit.bcirc/trn.transit.aRs)*180./np.pi)                     # Orbital inclination
  
  # Mask the star
  for j in range(len(x)):
    if (x[j]**2 + y[j]**2) < 1. and (z[j] > 0):
      x[j] = np.nan
      y[j] = np.nan
  
  # The star
  r = np.linspace(0,1,100)
  Ir = I(r,trn.limbdark)/I(0,trn.limbdark)
  
  for ri,Iri in zip(r[::-1],Ir[::-1]):
    star = pl.Circle((0, 0), ri, color=str(0.95*Iri), alpha=1.)
    ax2.add_artist(star)

  # Inset: Limb darkening
  if ldplot:
    if compact:
      inset1 = pl.axes([0.145, 0.32, .09, .1])
    else:
      inset1 = fig.add_axes([0.725,0.3,0.2,0.15])
    inset1.plot(r,Ir,'k-')
    pl.setp(inset1, xlim=(-0.1,1.1), ylim=(-0.1,1.1), xticks=[0,1], yticks=[0,1])
    for tick in inset1.xaxis.get_major_ticks() + inset1.yaxis.get_major_ticks():
      tick.label.set_fontsize(8)
    inset1.set_ylabel(r'I/I$_0$', fontsize=8, labelpad=-8)
    inset1.set_xlabel(r'r/R$_\star$', fontsize=8, labelpad=-8)
    inset1.set_title('Limb Darkening', fontweight='bold', fontsize=9)
    
  # Inset: Top view of orbit
  if compact:
    inset2 = pl.axes([0.135, 0.115, .1, .1])
  else:
    inset2 = fig.add_axes([0.725,0.1,0.2,0.15])
  pl.setp(inset2, xticks=[], yticks=[])
  trn.transit.bcirc = trn.transit.aRs                                                 # This ensures we are face-on
  try:
    trn.Compute()
  except Exception as e:
    if str(e) == "Object does not transit the star.":
      pass
    else: raise Exception(e)
  xp = trn.arrays.x
  yp = trn.arrays.y
  inset2.plot(xp, yp, '-', color='DarkBlue', alpha=0.5)
  # Draw some invisible dots at the corners to set the window size
  xmin, xmax, ymin, ymax = np.nanmin(xp), np.nanmax(xp), np.nanmin(yp), np.nanmax(yp)
  xrng = xmax - xmin
  yrng = ymax - ymin
  xmin -= 0.1*xrng; xmax += 0.1*xrng;
  ymin -= 0.1*yrng; ymax += 0.1*yrng;
  inset2.scatter([xmin,xmin,xmax,xmax], [ymin,ymax,ymin,ymax], alpha = 0.)
  # Plot the star
  for ri,Iri in zip(r[::-10],Ir[::-10]):
    star = pl.Circle((0, 0), ri, color=str(0.95*Iri), alpha=1.)
    inset2.add_artist(star)
  # Plot the planet
  ycenter = yp[np.where(np.abs(xp) == np.nanmin(np.abs(xp)))][0]
  while ycenter > 0:
    xp[np.where(np.abs(xp) == np.nanmin(np.abs(xp)))] = np.nan
    ycenter = yp[np.where(np.abs(xp) == np.nanmin(np.abs(xp)))][0]
  planet = pl.Circle((0, ycenter), trn.transit.RpRs, color='DarkBlue', alpha=1.)
  inset2.add_artist(planet)
  inset2.set_title('Top View', fontweight='bold', fontsize=9)
  inset2.set_aspect('equal','datalim')
  
  # The orbit itself
  with np.errstate(invalid='ignore'):
    ax2.plot(x, y, '-', color='DarkBlue', lw = 1. if per < 30. else 
                                               max(1. - (per - 30.) / 100., 0.3) )

  # The planet
  with np.errstate(invalid = 'ignore'):
    ycenter = y[np.where(np.abs(x) == np.nanmin(np.abs(x)))][0]
  while ycenter > 0:
    x[np.where(np.abs(x) == np.nanmin(np.abs(x)))] = np.nan
    ycenter = y[np.where(np.abs(x) == np.nanmin(np.abs(x)))][0]
  planet = pl.Circle((0, ycenter), trn.transit.RpRs, color='DarkBlue', alpha=1.)
  ax2.add_artist(planet)
  
  # Force aspect
  if xlim is None:
    xlim = 1.1 * max(np.nanmax(x), np.nanmax(-x))
  ax2.set_ylim(-xlim/3.2,xlim/3.2)
  ax2.set_xlim(-xlim,xlim)
  
  ax2.set_xlabel(r'X (R$_\star$)', fontweight='bold')
  ax2.set_ylabel(r'Y (R$_\star$)', fontweight='bold')
  ax1.set_title(plottitle, fontsize=12)
  
  if not compact:
    rect = 0.725,0.55,0.2,0.35
    ax3 = fig.add_axes(rect)
    ax3.xaxis.set_visible(False)
    ax3.yaxis.set_visible(False)

    # Table of parameters
    ltable = [ r'$P:$',
               r'$e:$',
               r'$i:$',
               r'$\omega:$',
               r'$\rho_\star:$',
               r'$M_p:$',
               r'$R_p:$',
               r'$q_1:$',
               r'$q_2:$']
    rtable = [ r'$%.4f\ \mathrm{days}$' % trn.transit.per,
               r'$%.5f$' % trn.transit.ecc,
               r'$%.4f^\circ$' % inc,
               r'$%.3f^\circ$' % (trn.transit.w*180./np.pi),
               r'$%.5f\ \mathrm{g/cm^3}$' % trn.transit.rhos,
               r'$%.5f\ M_\star$' % trn.transit.MpMs,
               r'$%.5f\ R_\star$' % trn.transit.RpRs,
               r'$%.5f$' % trn.limbdark.q1,
               r'$%.5f$' % trn.limbdark.q2]
    yt = 0.875
    for l,r in zip(ltable, rtable):
      ax3.annotate(l, xy=(0.25, yt), xycoords="axes fraction", ha='right', fontsize=16)
      ax3.annotate(r, xy=(0.35, yt), xycoords="axes fraction", fontsize=16)
      yt -= 0.1

  return fig