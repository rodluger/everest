#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`kepler.py` - Main mission routines
-------------------------------------------

Implements several routines specific to the `Kepler` mission.

.. todo:: This mission is not yet supported.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import logging
log = logging.getLogger(__name__)

__all__ = ['Setup', 'Season', 'Breakpoints', 'GetData', 'GetNeighbors', 
           'Statistics', 'TargetDirectory', 'HasShortCadence', 
           'InjectionStatistics', 'HDUCards', 'FITSFile', 'FITSUrl', 'CDPP',
           'GetTargetCBVs', 'FitCBVs', 'PlanetStatistics']

def Setup():
  '''
  Called when the code is installed.
  
  '''
  
  raise NotImplementedError('This mission is not yet supported.')

def Season(ID, **kwargs):
  '''
  Returns the season number for a given target.
  
  '''
  
  raise NotImplementedError('This mission is not yet supported.')

def Breakpoints(ID, cadence = 'lc', **kwargs):  
  '''
  
  Returns the location of the breakpoints for a given target.

  :param int ID: The target ID number
  :param str cadence: The light curve cadence. Default `lc`
  
  .. note :: The number corresponding to a given breakpoint is the number \
            of cadences *since the beginning of the campaign*.
  
  '''
  
  raise NotImplementedError('This mission is not yet supported.')

def CDPP(flux, mask = [], cadence = 'lc'):
  '''
  Compute the CDPP metric for a given timeseries.
  
  :param array_like flux: The flux array to compute the CDPP for
  :param array_like mask: The indices to be masked
  :param str cadence: The light curve cadence. Default `lc`
  
  '''
  
  raise NotImplementedError('This mission is not yet supported.')
    
def GetData(ID, season = None, cadence = 'lc', clobber = False, delete_raw = False, 
            aperture_name = None, saturated_aperture_name = None,
            max_pixels = None, download_only = False, saturation_tolerance = None, 
            bad_bits = None, **kwargs):
  '''
  Returns a :py:obj:`DataContainer` instance with the raw data for the target.
  
  :param int ID: The target ID number
  :param int season: The observing season. Default :py:obj:`None`
  :param str cadence: The light curve cadence. Default `lc`
  :param bool clobber: Overwrite existing files? Default :py:obj:`False`
  :param bool delete_raw: Delete the FITS TPF after processing it? Default :py:obj:`False`
  :param str aperture_name: The name of the aperture to use. Select `custom` to call \
         :py:func:`GetCustomAperture`. Default :py:obj:`None`
  :param str saturated_aperture_name: The name of the aperture to use if the target is \
         saturated. Default :py:obj:`None`
  :param int max_pixels: Maximum number of pixels in the TPF. Default :py:obj:`None`
  :param bool download_only: Download raw TPF and return? Default :py:obj:`False`
  :param float saturation_tolerance: Target is considered saturated if flux is within \
         this fraction of the pixel well depth. Default :py:obj:`None`
  :param array_like bad_bits: Flagged :py:obj`QUALITY` bits to consider outliers when \
         computing the model. Default :py:obj:`None`
  
  '''
  
  raise NotImplementedError('This mission is not yet supported.')

def GetNeighbors(ID, model = None, neighbors = None, mag_range = None, 
                 cdpp_range = None, aperture_name = None, 
                 cadence = 'lc', **kwargs):
  '''
  Return `neighbors` random bright stars on the same module as `EPIC`.
  
  :param int ID: The target ID number
  :param str model: The :py:obj:`everest` model name. Only used when imposing CDPP bounds. Default :py:obj:`None`
  :param int neighbors: Number of neighbors to return. Default None
  :param str aperture_name: The name of the aperture to use. Select `custom` to call \
         :py:func:`GetCustomAperture`. Default :py:obj:`None`
  :param str cadence: The light curve cadence. Default `lc`
  :param tuple mag_range: (`low`, `high`) values for the Kepler magnitude. Default :py:obj:`None`
  :param tuple cdpp_range: (`low`, `high`) values for the de-trended CDPP. Default :py:obj:`None`
  
  '''
  
  raise NotImplementedError('This mission is not yet supported.')

def PlanetStatistics(*args, **kwargs):
  '''
  Computes and plots the CDPP statistics comparison between `model` and
  `compare_to` for all known planets.

  '''
    
  raise NotImplementedError('This mission is not yet supported.')

def ShortCadenceStatistics(*args, **kwargs):
  '''
  Computes and plots the CDPP statistics comparison between short cadence
  and long cadence de-trended light curves
  
  '''
  
  raise NotImplementedError('This mission is not yet supported.')
    
def Statistics(*args, **kwargs):
  '''
  Computes and plots the CDPP statistics comparison between `model` and `compare_to` 
  for all long cadence light curves in a given campaign

  '''
  
  raise NotImplementedError('This mission is not yet supported.')

def HasShortCadence(ID, season = None):
  '''
  Returns `True` if short cadence data is available for this target.
  
  :param int ID: The target ID number
  :param int season: The season number. Default :py:obj:`None`
  
  '''
  
  raise NotImplementedError('This mission is not yet supported.')

def InjectionStatistics(*args, **kwargs):
  '''
  Computes and plots the statistics for injection/recovery tests.
  
  '''
  
  raise NotImplementedError('This mission is not yet supported.')

def HDUCards(headers, hdu = 0):
  '''
  Generates HDU cards for inclusion in the de-trended light curve FITS file.
  Used internally.
  
  '''
  
  raise NotImplementedError('This mission is not yet supported.') 

def TargetDirectory(ID, season, relative = False, **kwargs):
  '''
  Returns the location of the :py:mod:`everest` data on disk
  for a given target.
  
  :param ID: The target ID
  :param int season: The target season number
  :param bool relative: Relative path? Default :py:obj:`False`
  
  '''
  
  raise NotImplementedError('This mission is not yet supported.')

def FITSFile(ID, season, cadence = 'lc'):
  '''
  Returns the name of the FITS file for a given target.
  
  :param ID: The target ID
  :param int season: The target season number
  :param str cadence: The cadence type. Default `lc`
  
  '''
  
  raise NotImplementedError('This mission is not yet supported.')
    
def FITSUrl(ID, season):
  '''
  Returns the online path to the FITS file for a given target.
  
  :param ID: The target ID
  :param int season: The target season number
  
  '''
  
  raise NotImplementedError('This mission is not yet supported.')

def GetTargetCBVs(model):
  '''
  Returns the design matrix of CBVs for the given target.
  
  :param model: An instance of the :py:obj:`everest` model for the target
  
  '''
  
  raise NotImplementedError('This mission is not yet supported.')
  
def FitCBVs(model):
  '''
  Fits the CBV design matrix to the de-trended flux of a given target. This
  is called internally whenever the user accesses the :py:attr:`fcor`
  attribute.
  
  :param model: An instance of the :py:obj:`everest` model for the target
  
  '''
  
  raise NotImplementedError('This mission is not yet supported.')