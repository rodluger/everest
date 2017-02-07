#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import
from setuptools import setup, find_packages

# Hackishly inject a constant into builtins to enable importing of the
# module in "setup" mode. Stolen from `kplr`
import sys
if sys.version_info[0] < 3:
  import __builtin__ as builtins
else:
  import builtins
builtins.__EVEREST_SETUP__ = True
import everest

long_description = \
"""
EPIC Variability Extraction and Removal for Exoplanet Science Targets: 
A pipeline for de-trending `K2`, `Kepler`, and `TESS` light curves with 
pixel level decorrelation and Gaussian processes. The Python interface
allows easy access to the online EVEREST de-trended light curve catalog;
alternatively, users can de-trend their own light curves with customized
settings. Read the documentation at https://github.com/rodluger/everest
"""

# Setup!
setup(name = 'everest-pipeline',
      version = everest.__subversion__,
      description = 'EPIC Variability Extraction and Removal for Exoplanet Science Targets',
      long_description = long_description,
      classifiers = [
                      'Development Status :: 5 - Production/Stable',
                      'License :: OSI Approved :: MIT License',
                      'Programming Language :: Python',
                      'Programming Language :: Python :: 3',
                      'Topic :: Scientific/Engineering :: Astronomy',
                    ],
      url = 'http://github.com/rodluger/everest',
      author = 'Rodrigo Luger',
      author_email = 'rodluger@uw.edu',
      license = 'MIT',
      packages = ['everest', 'everest.missions', 'everest.missions.k2',
                  'everest.missions.kepler', 'everest.missions.tess'],
      install_requires = [
                          'numpy>=1.8',
                          'scipy',
                          'matplotlib',
                          'george==0.2.1',
                          'six',
                          'astropy',
                          'pysyzygy>=0.0.1',
                          'k2plr>=0.2.4',
                          'PyPDF2'
                         ],
      scripts=['bin/everest', 'bin/everest-stats', 'bin/everest-status'],
      include_package_data = True,
      zip_safe = False,
      test_suite='nose.collector',
      tests_require=['nose']
      )

# Set up the individual missions
from everest import missions
for mission in missions.Missions:
  getattr(missions, mission).Setup()
