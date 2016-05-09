#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
inject_control.py
-----------------

Runs **control** injection/recovery tests on the sample of ~2,000
random EPIC targets returned by ``GetK2InjectionTestStars()``.
In this control group, the injection/recovery is performed on the
**de-trended** light curve. There will be no bias (by construction)
but the width of the recovered depth distributions is useful when
comparing to the actual injection runs. These runs generate the
red histograms in the injection/recovery plot in the paper.


'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
EVEREST_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, EVEREST_ROOT)
from everest.data import GetK2InjectionTestStars
from everest.transit import Transit
from everest.compute import GetTransitDepth
import numpy as np

count = 0
bad = []
for EPIC in GetK2InjectionTestStars():
  
  # Get the de-trended light curve (no injections)
  for campaign in range(10):
    outdir = os.path.join(EVEREST_ROOT, 'output', 'C%02d' % campaign, '%s' % EPIC, 'default')
    if os.path.exists(os.path.join(outdir, 'data.npz')):
      break
  if campaign == 9:
    # No data for this target
    continue
  
  if os.path.exists(os.path.join(outdir, '%d.ctrl6.inj' % EPIC)):
    # We've done this one already
    continue
  
  try:
    data = np.load(os.path.join(outdir, 'data.npz'))
    time = data['time']
    fpld = data['fpld']
  except:
    bad.append(str(EPIC))
    continue
  
  count += 1
  n = 1
  for depth in [0.01, 0.001, 0.0001]:
    for mask in [True, False]:

      # The injection info
      inject = dict(t0 = 0., 
                    per = 3 + 7 * np.random.random(),   # Random periods between 3 and 10 days
                    depth = depth, 
                    dur = 0.1,
                    mask = mask)

      # Recover the depth
      inject.update({'everest': GetTransitDepth(time, fpld * Transit(time, **inject), inject, buf = 5, order = 2)})

      # Save recovery info to disk
      with open(os.path.join(outdir, '%d.ctrl%d.inj' % (EPIC, n)), 'w') as injfile:
        if inject['mask']:
          print('Masked:            True', file = injfile)
        else:
          print('Masked:            False', file = injfile)
        print('Injected period:   %.2f d' % inject['per'], file = injfile)
        print('Injected duration: %.2f d' % inject['dur'], file = injfile)
        print('Injected depth:    %.5f' % inject['depth'], file = injfile)
        print('Everest depth:     %.8f +/- %.8f' % (inject['everest'][0], 
                                  inject['everest'][1]), file = injfile)
                                  
      # Increment tag
      n += 1
      
print('%d targets processed.' % count)
print('ERRORS:', ', '.join(bad))