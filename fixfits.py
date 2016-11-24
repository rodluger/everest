from everest.missions.k2 import TargetDirectory, FITSFile
from everest.missions.k2.aux import GetK2Campaign
import os
try:
  import pyfits
except ImportError:
  try:
    import astropy.io.fits as pyfits
  except ImportError:
    raise Exception('Please install the `pyfits` package.')
campaign = 2
stars = GetK2Campaign(campaign, epics_only = True)
for i, ID in enumerate(stars):
  print('Processing EPIC %d (%d/%d)...' % (ID, i + 1, len(stars)))
  file = os.path.join(TargetDirectory(ID, campaign), FITSFile(ID, campaign))
  if os.path.exists(file):
    with pyfits.open(file, 'update') as f:
      if len(f) == 6:
        h = f[1].header
        if h.get('BRKPT1') is None:
          continue
        h.set('BRKPT01', h.get('BRKPT1'), 'Light curve breakpoint', 'BRKPT1')
        h.remove('BRKPT1')
        h.set('BRKPT02', h.get('BRKPT2'), 'Light curve breakpoint', 'BRKPT2')
        h.remove('BRKPT2')
        h.set('LAMB0101', h.get('LAMBDA11'), 'Cross-validation parameter', 'LAMBDA11')
        h.remove('LAMBDA11')
        h.set('LAMB0102', h.get('LAMBDA12'), 'Cross-validation parameter', 'LAMBDA12')
        h.remove('LAMBDA12')
        h.set('LAMB0103', h.get('LAMBDA13'), 'Cross-validation parameter', 'LAMBDA13')
        h.remove('LAMBDA13')
        h.set('LAMB0201', h.get('LAMBDA21'), 'Cross-validation parameter', 'LAMBDA21')
        h.remove('LAMBDA21')
        h.set('LAMB0202', h.get('LAMBDA22'), 'Cross-validation parameter', 'LAMBDA22')
        h.remove('LAMBDA22')
        h.set('LAMB0203', h.get('LAMBDA23'), 'Cross-validation parameter', 'LAMBDA23')
        h.remove('LAMBDA23')
        continue
    print('Removing file...')
    os.remove(file)