#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`feature_test.py` - Feature testing
-------------------------------------------

A simple test that plots a bunch of different light curves.
The output should be similar to this figure:

.. figure:: ../pysyzygy/img/feature_test.png
    :width: 600px
    :align: center
    :height: 100px
    :alt: alternate text
    :figclass: align-center

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import pysyzygy as ps
import numpy as np
import matplotlib.pyplot as pl

if __name__ == '__main__':

  fig, ax = pl.subplots(5, 5, figsize = (24,16))
  fig.subplots_adjust(wspace = 0, hspace = 0, left = 0.01, right = 0.99, bottom = 0.01, top = 0.99)
  ax = ax.flatten()
  for axis in ax:
    axis.xaxis.set_ticklabels([])
    axis.yaxis.set_ticklabels([])

  # Different radii
  time = np.linspace(-0.5,0.5,1000)
  for RpRs in [0.1, 0.09, 0.08]:
    trn = ps.Transit(RpRs = RpRs)
    ax[0].plot(time, trn(time), label = 'RpRs = %.2f' % RpRs)
  ax[0].legend(loc = 'lower left', fontsize = 8)
  ax[0].margins(0.,0.2)
  ax[0].annotate('DIFFERENT RADII', xy = (0.05, 0.9), xycoords = 'axes fraction')
  
  # Different periods
  time = np.linspace(-0.5,0.5,1000)
  for per in [10., 20., 30.]:
    trn = ps.Transit(per = per)
    ax[1].plot(time, trn(time), label = 'per = %.0f' % per)
  ax[1].legend(loc = 'lower left', fontsize = 8)
  ax[1].margins(0.,0.2)  
  ax[1].annotate('DIFFERENT PERIODS', xy = (0.05, 0.9), xycoords = 'axes fraction')

  # Different impact params
  time = np.linspace(-0.5,0.5,1000)
  for b in [0., 0.25, 0.5, 0.75, 0.99]:
    trn = ps.Transit(b = b)
    ax[2].plot(time, trn(time), label = 'b = %.2f' % b)
  ax[2].legend(loc = 'lower left', fontsize = 8)
  ax[2].margins(0.,0.2) 
  ax[2].annotate('DIFFERENT IMPACT', xy = (0.05, 0.9), xycoords = 'axes fraction') 

  # Different densities
  time = np.linspace(-0.5,0.5,1000)
  for rhos in [1.4, 10.0, 0.1]:
    trn = ps.Transit(rhos = rhos)
    ax[3].plot(time, trn(time), label = 'rhos = %.2f' % rhos)
  ax[3].legend(loc = 'lower left', fontsize = 8)
  ax[3].margins(0.,0.2) 
  ax[3].annotate('DIFFERENT STELLAR DENSITIES', xy = (0.05, 0.9), xycoords = 'axes fraction') 

  # Different semi-major axes
  time = np.linspace(-0.5,0.5,1000)
  for aRs in [19.5, 5., 50.]:
    trn = ps.Transit(aRs = aRs)
    ax[4].plot(time, trn(time), label = 'aRs = %.2f' % aRs)
  ax[4].legend(loc = 'lower left', fontsize = 8)
  ax[4].margins(0.,0.2)
  ax[4].annotate('DIFFERENT SEMI-MAJOR AXES', xy = (0.05, 0.9), xycoords = 'axes fraction')  

  # Different masses
  time = np.linspace(-0.5,0.5,1000)
  for MpMs in [0., 0.1, 1.]:
    trn = ps.Transit(MpMs = MpMs)
    ax[5].plot(time, trn(time), label = 'MpMs = %.2f' % MpMs)
  ax[5].legend(loc = 'lower left', fontsize = 8)
  ax[5].margins(0.,0.2)  
  ax[5].annotate('DIFFERENT MASSES', xy = (0.05, 0.9), xycoords = 'axes fraction')

  # Different eccentricities
  time = np.linspace(-0.5,0.5,1000)
  for ecc in [0., 0.25, 0.5, 0.75, 0.9]:
    trn = ps.Transit(ecc = ecc)
    ax[6].plot(time, trn(time), label = 'ecc = %.2f' % ecc)
  ax[6].legend(loc = 'lower left', fontsize = 8)
  ax[6].margins(0.,0.2)
  ax[6].annotate('DIFFERENT ECCENTRICITIES', xy = (0.05, 0.9), xycoords = 'axes fraction')  

  # Different omegas
  time = np.linspace(-0.5,0.5,1000)
  for w in [0., np.pi/8, np.pi/4, np.pi/3, np.pi/2]:
    trn = ps.Transit(aRs = 5, ecc = 0.75, w = w)
    ax[7].plot(time, trn(time), label = 'w = %.2f' % w)
  ax[7].legend(loc = 'lower left', fontsize = 8)
  ax[7].margins(0.,0.2) 
  ax[7].annotate('DIFFERENT OMEGAS', xy = (0.05, 0.9), xycoords = 'axes fraction')
  
  # Different esw/ecw
  time = np.linspace(-0.5,0.5,1000)
  for esw, ecw in zip([0., 0.25, -0.25, 0.5], [0., -0.25, 0.5, 0.25]):
    trn = ps.Transit(aRs = 10, esw = esw, ecw = ecw)
    ax[8].plot(time, trn(time), label = 'esw = %.2f, ecw = %.2f' % (esw, ecw))
  ax[8].legend(loc = 'lower left', fontsize = 8)
  ax[8].margins(0.,0.2) 
  ax[8].annotate('DIFFERENT esinw/ecosw', xy = (0.05, 0.9), xycoords = 'axes fraction')  

  # Different limb darkening 1
  time = np.linspace(-0.5,0.5,1000)
  for u1, u2 in zip([0., 0.25, 0.5, 0.75], [0., 0.25, -0.5, -0.25]):
    trn = ps.Transit(aRs = 5, u1 = u1, u2 = u2)
    ax[9].plot(time, trn(time), label = 'u1 = %.2f, u2 = %.2f' % (u1, u2))
  ax[9].legend(loc = 'upper center', fontsize = 8)
  ax[9].margins(0.,0.2)  
  ax[9].annotate('DIFFERENT LD', xy = (0.05, 0.9), xycoords = 'axes fraction')
 
  # Different limb darkening 2
  time = np.linspace(-0.5,0.5,1000)
  for q1, q2 in zip([0., 0.25, 0.5, 0.75], [0., 0.5, 0.25, 0.1]):
    trn = ps.Transit(aRs = 5, q1 = q1, q2 = q2)
    ax[10].plot(time, trn(time), label = 'q1 = %.2f, q2 = %.2f' % (q1, q2))
  ax[10].legend(loc = 'upper center', fontsize = 8)
  ax[10].margins(0.,0.2)  
  ax[10].annotate('DIFFERENT LD', xy = (0.05, 0.9), xycoords = 'axes fraction')

  # Different t0
  time = np.linspace(-0.5,0.5,1000)
  for t0 in [-0.1, 0., 0.1]:
    trn = ps.Transit(t0 = t0)
    ax[11].plot(time, trn(time), label = 't0 = %.2f' % t0)
  ax[11].legend(loc = 'lower left', fontsize = 8)
  ax[11].margins(0.,0.2) 
  ax[11].annotate('DIFFERENT t0', xy = (0.05, 0.9), xycoords = 'axes fraction') 

  # Different times
  time = np.linspace(-0.5,0.5,1000)
  for times in [[-0.3, 0.3], [-0.4, 0., 0.4], [0.]]:
    trn = ps.Transit(times = times)
    ax[12].plot(time, trn(time), label = 'times = %s' % times)
  ax[12].legend(loc = 'lower left', fontsize = 8)
  ax[12].margins(0.,0.2)
  ax[12].annotate('DIFFERENT TRANSIT TIMES', xy = (0.05, 0.9), xycoords = 'axes fraction')

  # Different exposure times
  time = np.linspace(-0.5,0.5,1000)
  for exptime in [(1765.5/86400.), 10 * (1765.5/86400.), 0.2 * (1765.5/86400.)]:
    trn = ps.Transit(aRs = 5, exptime = exptime)
    ax[13].plot(time, trn(time), label = 'exptime = %.3f' % exptime)
  ax[13].legend(loc = 'lower left', fontsize = 8)
  ax[13].margins(0.,0.2)  
  ax[13].annotate('DIFFERENT EXP TIMES', xy = (0.05, 0.9), xycoords = 'axes fraction')
 
  # Different number of exposure points
  time = np.linspace(-0.5,0.5,1000)
  for exppts in [4, 10, 100]:
    trn = ps.Transit(aRs = 5, exptime = 0.1, exppts = exppts)
    ax[14].plot(time, trn(time), label = 'exppts = %d' % exppts)
  ax[14].legend(loc = 'lower left', fontsize = 8)
  ax[14].margins(0.,0.2)   
  ax[14].annotate('DIFFERENT NUMBER OF EXP PTS', xy = (0.05, 0.9), xycoords = 'axes fraction')

  # Different binning method
  time = np.linspace(-0.5,0.5,1000)
  for binmethod in [ps.RIEMANN, ps.TRAPEZOID]:
    trn = ps.Transit(aRs = 5, binmethod = binmethod, exppts = 4, exptime = 0.1)
    ax[15].plot(time, trn(time), label = 'binmethod = %d' % binmethod)
  ax[15].legend(loc = 'lower left', fontsize = 8)
  ax[15].margins(0.,0.2)   
  ax[15].annotate('DIFFERENT BIN METHOD', xy = (0.05, 0.9), xycoords = 'axes fraction')

  # Different integration method (curves should be identical)
  time = np.linspace(-0.5,0.5,1000)
  for intmethod in [ps.SMARTINT, ps.SLOWINT]:
    trn = ps.Transit(aRs = 5, intmethod = intmethod, exppts = 4, exptime = 0.1)
    ax[16].plot(time, trn(time), label = 'intmethod = %d' % intmethod)
  ax[16].legend(loc = 'lower left', fontsize = 8)
  ax[16].margins(0.,0.2)   
  ax[16].annotate('SHOULD BE IDENTICAL', xy = (0.05, 0.9), xycoords = 'axes fraction')

  # Different kepler solver (curves should be identical)
  time = np.linspace(-0.5,0.5,1000)
  for kepsolver in [ps.MDFAST, ps.NEWTON]:
    trn = ps.Transit(aRs = 5, ecc = 0.75, w = np.pi/10, kepsolver = kepsolver)
    ax[17].plot(time, trn(time), label = 'kepsolver = %d' % kepsolver)
  ax[17].legend(loc = 'lower left', fontsize = 8)
  ax[17].margins(0.,0.2)   
  ax[17].annotate('SHOULD BE IDENTICAL', xy = (0.05, 0.9), xycoords = 'axes fraction')

  # Binned/unbinned flux
  time = np.linspace(-0.5,0.5,1000)
  for param in ['binned', 'unbinned']:
    trn = ps.Transit(exptime = 0.1)
    ax[18].plot(time, trn(time, param = param), label = 'param = %s' % param)
  ax[18].legend(loc = 'lower left', fontsize = 8)
  ax[18].margins(0.,0.2)
  ax[18].annotate('BINNED/UNBINNED', xy = (0.05, 0.9), xycoords = 'axes fraction') 

  # Mean anomaly / eccentric anomaly / true anomaly
  time = np.linspace(-6,6,1000)
  trn = ps.Transit(per = 5., ecc = 0.4, fullorbit = True, maxpts = 20000)
  for param in ['M', 'E', 'f']:
    ax[19].plot(time, trn(time, param = param), label = 'param = %s' % param)
  ax[19].legend(loc = 'lower left', fontsize = 8)
  ax[19].margins(0.,0.2) 
  ax[19].annotate('ANOMALIES', xy = (0.05, 0.9), xycoords = 'axes fraction')

  # radius / impact parameter
  time = np.linspace(-6,6,1000)
  trn = ps.Transit(per = 5., ecc = 0.4, fullorbit = True, maxpts = 20000)
  for param in ['r', 'b']:
    ax[20].plot(time, trn(time, param = param), label = 'param = %s' % param)
  ax[20].legend(loc = 'lower left', fontsize = 8)
  ax[20].margins(0.,0.2) 
  ax[20].annotate('r and b', xy = (0.05, 0.9), xycoords = 'axes fraction')

  # xyz
  time = np.linspace(-5,5,1000)
  trn = ps.Transit(ecc = 0.5, w = np.pi/8, b = 0.5, per = 5., fullorbit = True, maxpts = 20000)
  ax[21].plot(trn(time, 'x'), trn(time, 'y'), label = 'XY')
  ax[21].plot(trn(time, 'x'), trn(time, 'z'), label = 'XZ')
  ax[21].plot(trn(time, 'y'), trn(time, 'z'), label = 'YZ')
  ax[21].legend(loc = 'lower left', fontsize = 8)
  ax[21].margins(0.,0.2) 
  ax[21].annotate('ORBIT', xy = (0.05, 0.9), xycoords = 'axes fraction')

  # Different omega, top view orbit
  time = np.linspace(-5,5,1000)
  for w in np.arange(4) * np.pi / 4:
    trn = ps.Transit(w = w, ecc = 0.5, b = 0., per = 5., fullorbit = True, maxpts = 20000)
    ax[22].plot(trn(time, 'x'), trn(time, 'z'), label = r'w = %.0f$^\circ$' % (w * 180./np.pi))
  ax[22].axvline(0, alpha = 0.5)
  ax[22].axhline(0, alpha = 0.5)
  ax[22].legend(loc = 'lower left', fontsize = 8)
  ax[22].margins(0.,0.2) 
  ax[22].annotate('TOP VIEW', xy = (0.05, 0.9), xycoords = 'axes fraction')

  # Show the plot
  pl.show()
