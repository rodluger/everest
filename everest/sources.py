#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`sources.py` - Nearby sources
-------------------------------------

Adapted from some of the `PyKE tools <http://keplergo.arc.nasa.gov/PyKE.shtml>`_.
Queries the `MAST` archive for sources within or near the aperture of a given
target and returns their magnitudes and positions on the detector.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np
import urllib
import re
import kplr
import logging
log = logging.getLogger(__name__)

class Source(object):
  '''
  A generic source object.
  
  '''
  
  def __init__(self, **kwargs):
    self.__dict__.update(kwargs)
    
  def __repr__(self):
    return "<Source: EPIC %d>" % self.epic
    
def MASTRADec(ra, dec, darcsec):
  '''
  Detector location retrieval based upon RA and Dec.
  Adapted from 
  `PyKE <http://keplergo.arc.nasa.gov/PyKE.shtml>`_.
  
  '''
  
  # coordinate limits
  darcsec /= 3600.0
  ra1 = ra - darcsec / np.cos(dec * np.pi / 180)
  ra2 = ra + darcsec / np.cos(dec * np.pi / 180)
  dec1 = dec - darcsec
  dec2 = dec + darcsec
 
  # build mast query
  url  = 'http://archive.stsci.edu/k2/epic/search.php?'
  url += 'action=Search'
  url += '&k2_ra=' + str(ra1) + '..' + str(ra2)
  url += '&k2_dec=' + str(dec1) + '..' + str(dec2)
  url += '&max_records=10000'
  url += '&selectedColumnsCsv=id,k2_ra,k2_dec,kp'
  url += '&outputformat=CSV'

  # retrieve results from MAST
  try:
    lines = urllib.request.urlopen(url)
  except:
    log.warn('Unable to retrieve source data from MAST.')
    lines = ''

  # collate nearby sources
  epicid = []
  kepmag = []
  ra = []
  dec = []
  for line in lines:
  
    line = line.strip().decode('ascii')
    
    if (len(line) > 0 and 'EPIC' not in line and 'integer' not in line and
                          'no rows found' not in line):
        
      out = line.split(',')
      r, d = sex2dec(out[1], out[2])
      epicid.append(int(out[0]))
      kepmag.append(float(out[3]))
      ra.append(r)
      dec.append(d)
  
  epicid = np.array(epicid)
  kepmag = np.array(kepmag)
  ra = np.array(ra)
  dec = np.array(dec)

  return epicid, ra, dec, kepmag


def sex2dec(ra, dec):
  '''
  Convert sexadecimal hours to decimal degrees. Adapted from 
  `PyKE <http://keplergo.arc.nasa.gov/PyKE.shtml>`_.
  
  :param float ra: The right ascension
  :param float dec: The declination
  
  :returns: The same values, but in decimal degrees
  
  '''

  ra = re.sub('\s+','|', ra.strip())
  ra = re.sub(':','|', ra.strip())
  ra = re.sub(';','|', ra.strip())
  ra = re.sub(',','|', ra.strip())
  ra = re.sub('-','|', ra.strip())
  ra = ra.split('|')
  outra = (float(ra[0]) + float(ra[1]) / 60. + float(ra[2]) / 3600.) * 15.0

  dec = re.sub('\s+','|', dec.strip())
  dec = re.sub(':','|', dec.strip())
  dec = re.sub(';','|', dec.strip())
  dec = re.sub(',','|', dec.strip())
  dec = dec.split('|')

  if float(dec[0]) > 0.0:
      outdec = float(dec[0]) + float(dec[1]) / 60. + float(dec[2]) / 3600.
  else:
      outdec = float(dec[0]) - float(dec[1]) / 60. - float(dec[2]) / 3600.

  return outra, outdec

def GetSources(EPIC):
  '''
  Grabs the EPIC coordinates from the TPF and searches MAST
  for other EPIC targets within the same aperture.
  
  :param int EPIC: The 9-digit `EPIC` number of the target
  
  :returns: A list of :py:class:`Source` instances containing \
            other `EPIC` targets within or close to this target's aperture
  '''
  
  client = kplr.API()
  star = client.k2_star(EPIC)
  tpf = star.get_target_pixel_files()[0]
  with tpf.open() as f:
    crpix1 = f[2].header['CRPIX1']
    crpix2 = f[2].header['CRPIX2']
    crval1 = f[2].header['CRVAL1']  
    crval2 = f[2].header['CRVAL2'] 
    cdelt1 = f[2].header['CDELT1']   
    cdelt2 = f[2].header['CDELT2']
    pc1_1 = f[2].header['PC1_1']
    pc1_2 = f[2].header['PC1_2']
    pc2_1 = f[2].header['PC2_1']
    pc2_2 = f[2].header['PC2_2']
    pc = np.array([[pc1_1, pc1_2], [pc2_1, pc2_2]])
    pc = np.linalg.inv(pc)
    crpix1p = f[2].header['CRPIX1P']
    crpix2p = f[2].header['CRPIX2P']
    crval1p = f[2].header['CRVAL1P']  
    crval2p = f[2].header['CRVAL2P'] 
    cdelt1p = f[2].header['CDELT1P']   
    cdelt2p = f[2].header['CDELT2P']
    darcsec = 4 * max(f[2].data.shape)
    
  epicid, ra, dec, kepmag = MASTRADec(star.k2_ra, star.k2_dec, darcsec)
  sources = []
  for i, epic in enumerate(epicid):
    dra = (ra[i] - crval1) * np.cos(np.radians(dec[i])) / cdelt1
    ddec = (dec[i] - crval2) / cdelt2
    sx = pc[0,0] * dra + pc[0,1] * ddec + crpix1 + crval1p - 1.0
    sy = pc[1,0] * dra + pc[1,1] * ddec + crpix2 + crval2p - 1.0
    
    sources.append(Source(epic = epic, x = sx, y = sy, kepmag = kepmag[i],
                          x0 = crval1p, y0 = crval2p))
    
  return sources