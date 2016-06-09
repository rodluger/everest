#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, unicode_literals
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

# Get the long description from the README
def readme():
  with open('README.md') as f:
    return f.read()

# Setup!
setup(name = 'everest',
      version = everest.__version__,
      description = 'EPIC Variability Extraction and Removal for Exoplanet Science Targets',
      long_description = readme(),
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
      packages = ['everest', 'everest.usertools'],
      install_requires = [
                          'numpy',
                          'scipy',
                          'matplotlib',
                          'george',
                          'sklearn',
                          'astroML',
                          'six',
                          pyfits,
                          'pysyzygy>=0.0.1',
                          'k2plr==0.2.1'
                         ],
      dependency_links = [
                          'https://github.com/rodluger/pysyzygy/tarball/master#egg=pysyzygy-0.0.1',
                          'https://github.com/rodluger/k2plr/tarball/master#egg=k2plr-0.2.1'
                         ],
      include_package_data = True,
      zip_safe = False)