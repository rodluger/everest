#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
inject_control.py
-----------------

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
for EPIC in GetK2InjectionTestStars():
  
  # Get the de-trended light curve (no injections)
  for campaign in range(10):
    outdir = os.path.join(EVEREST_ROOT, 'output', 'C%02d' % campaign, '%s' % EPIC, 'default')
    if os.path.exists(os.path.join(outdir, 'data.npz')):
      count += 1
      break
  if campaign == 9:
    # No data for this target
    continue
  
  data = np.load(os.path.join(outdir, 'data.npz'))
  time = data['time']
  fpld = data['fpld']
  
  for depth in [0.01, 0.001, 0.0001]:
    for mask in [True, False]:

      # The injection info
      inject = dict(t0 = 0., 
                    per = 3.56789, 
                    depth = depth, 
                    dur = 0.1,
                    mask = mask)

      # Recover the depth
      inject.update({'everest': GetTransitDepth(time, fpld * Transit(time, **inject), inject, buf = 5, order = 2)})

      # Save recovery info to disk
      with open(os.path.join(outdir, '%d.ctrl.inj' % EPIC), 'w') as injfile:
        if inject['mask']:
          print('Masked:            True', file = injfile)
        else:
          print('Masked:            False', file = injfile)
        print('Injected period:   %.2f d' % inject['per'], file = injfile)
        print('Injected duration: %.2f d' % inject['dur'], file = injfile)
        print('Injected depth:    %.5f' % inject['depth'], file = injfile)
        print('Everest depth:     %.8f +/- %.8f' % (inject['everest'][0], 
                                  inject['everest'][1]), file = injfile)

print('%d targets processed.' % count)