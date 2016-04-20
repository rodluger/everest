#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
plot.py
-------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import matplotlib.pyplot as pl
import numpy as np
import os
import logging
log = logging.getLogger(__name__)

def Plot(data):
  '''
  
  '''
  
  EPIC = data['EPIC']
  outdir = data['outdir'][()]
  
  # Plot the apertures
  log.info('Plotting the apertures...')
  if not os.path.exists(os.path.join(outdir, 'aper.png')):
    fig, ax = PlotApertures(EPIC, data)
    fig.savefig(os.path.join(outdir, 'aper.png'))
    pl.close()

  # Plot the outliers
  log.info('Plotting the outliers...')
  if not os.path.exists(os.path.join(outdir, 'outliers.png')):
    fig, ax = PlotOutliers(EPIC, data)
    fig.savefig(os.path.join(outdir, 'outliers.png'))
    pl.close()  

  # Plot the apertures
  log.info('Plotting the GP...')
  if not os.path.exists(os.path.join(outdir, 'acor.png')):
    fig, ax = PlotGP(EPIC, data)
    fig.savefig(os.path.join(outdir, 'acor.png'))
    pl.close()

  # Plot the scatter curves
  log.info('Plotting the scatter curve...')
  if not os.path.exists(os.path.join(outdir, 'scatter.png')):
    try:
      fig, ax = PlotScatter(EPIC, data)
      fig.savefig(os.path.join(outdir, 'scatter.png'))
      pl.close()
    except:
      log.error('An error occurred while plotting the scatter curve.')

  # Plot the detrended data
  log.info('Plotting the detrended data...')
  if not os.path.exists(os.path.join(outdir, 'detrended.png')):
    fig, ax = PlotDetrended(EPIC, data)
    fig.savefig(os.path.join(outdir, 'detrended.png'))
    pl.close()
  
  # Plot the folded data
  log.info('Plotting the folded data...')
  try:
    PlotFolded(EPIC, data, outdir)
  except:
    log.error('An error occurred while plotting the folded data.')

  # Convert the images to a single PDF and remove the .png files
  images = []
  for img in ['detrended.png', 'folded_01.png', 'folded_02.png', 'folded_03.png', 
              'folded_04.png', 'folded_05.png', 'folded_06.png', 'folded_07.png', 
              'folded_EB01.png', 'folded_EB02.png', 'aper.png', 
              'scatter.png', 'acor.png', 'outliers.png']:
    if os.path.exists(os.path.join(outdir, img)):
      images.append(os.path.join(outdir, img))
  if len(images):
    try:
      subprocess.call(['convert'] + images + [os.path.join(outdir, '%s.pdf' % EPIC)])
      for image in images:
        os.remove(image)
    except:
      log.warn('Unable to generate PDF.')