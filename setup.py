#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import
from setuptools import setup, find_packages

# Hackishly inject a constant into builtins to enable importing of the
# module. Stolen from `kplr`
import sys
if sys.version_info[0] < 3:
  import __builtin__ as builtins
else:
  import builtins
builtins.__EVEREST_SETUP__ = True
import everest

# Check if the user has `pyfits` as part
# of the `astropy` distribution...
try:
  import astropy.io.fits
  pyfits = 'astropy'
except:
  pyfits = 'pyfits'

long_description = \
"""
EPIC Variability Extraction and Removal for Exoplanet Science Targets: 
A pipeline for de-trending `K2` light curves with pixel level decorrelation 
and Gaussian processes. Here you'll find the Python code used to generate 
the `everest` catalog, as well as tools for accessing and interacting 
with the de-trended light curves.
"""

# Setup!
setup(name = 'everest-pipeline',
      version = everest.__version__,
      description = 'EPIC Variability Extraction and Removal for Exoplanet Science Targets',
      long_description = long_description,
      classifiers = [
                      'Development Status :: 3 - Alpha',
                      'License :: OSI Approved :: MIT License',
                      'Programming Language :: Python',
                      'Programming Language :: Python :: 3',
                      'Topic :: Scientific/Engineering :: Astronomy',
                    ],
      url = 'http://github.com/rodluger/everest',
      author = 'Rodrigo Luger',
      author_email = 'rodluger@uw.edu',
      license = 'MIT',
      packages = ['everest'],
      install_requires = [
                          'numpy>=1.8',
                          'scipy',
                          'matplotlib',
                          'george==0.2.1',
                          'six',
                          pyfits,
                          'pysyzygy>=0.0.1',
                          'k2plr>=0.2.2',
                          'PyPDF2'
                         ],
      #dependency_links = [
      #                    'https://github.com/rodluger/pysyzygy/tarball/master#egg=pysyzygy-0.0.1',
      #                    'https://github.com/rodluger/k2plr/tarball/master#egg=k2plr-0.2.2'
      #                   ],
      scripts=['bin/everest', 'bin/everest-stats', 'bin/everest-status'],
      include_package_data = True,
      zip_safe = False,
      #test_suite='nose.collector',
      #tests_require=['nose']
      )