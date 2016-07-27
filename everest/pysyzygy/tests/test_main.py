#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test_main.py
------------

'''


import numpy as np
from pysyzygy.transit import Transit

def test_main():
  '''
  
  '''
  
  time = np.linspace(-0.5,0.5,1000)
  trn = Transit(per = 5., RpRs = 0.1, ecw = 0.5, esw = -0.5, b = 0., fullorbit = True, maxpts = 20000)
  
  params = ['unbinned', 'binned', 'E', 'f', 'r', 'b', 'x', 'y', 'z']
  truths = [0.998019215638, 0.998019220105, 1.51688474691, 2.27060195612, 
            11.8738101592, 2.60177797032, 0.0666039985104, 0.0, -11.4058661686]
  
  for param, truth in zip(params, truths):
    y = trn(time, param = param)
    i = np.trapz(y, time)
    np.testing.assert_array_almost_equal(i, truth)