#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
run.py
------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .plot import Plot
from .compute import Compute
import traceback
import sys
import pdb
import logging
log = logging.getLogger(__name__)

def Run(EPIC, debug = True, **kwargs):
  '''
  Compute and plot data for a given target.
  
  '''
  
  if debug:
    sys.excepthook = ExceptionHookPDB
  else:
    sys.excepthook = ExceptionHook
  
  data = Compute(EPIC, **kwargs)
  Plot(data)

def ExceptionHook(exctype, value, tb):
  '''
  A custom exception handler.
  
  '''
  
  for line in traceback.format_exception_only(exctype, value):
    log.error(line.replace('\n', ''))
  for line in traceback.format_tb(tb):
    log.error(line.replace('\n', ''))
  sys.__excepthook__(exctype, value, tb)

def ExceptionHookPDB(exctype, value, tb):
  '''
  A custom exception handler, with PDB post-mortem for debugging.
  
  '''
  
  for line in traceback.format_exception_only(exctype, value):
    log.error(line.replace('\n', ''))
  for line in traceback.format_tb(tb):
    log.error(line.replace('\n', ''))
  sys.__excepthook__(exctype, value, tb)
  pdb.pm()