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
builtins.__PYSYZYGY_SETUP__ = True
import pysyzygy

# Get the long description from the README
def readme():
  with open('README.md') as f:
    return f.read()

# Setup!
setup(name = 'pysyzygy',
      version = pysyzygy.__version__,
      description = 'Transit modeling in Python',
      long_description = readme(),
      classifiers = [
                      'Development Status :: 3 - Alpha',
                      'License :: OSI Approved :: MIT License',
                      'Programming Language :: Python',
                      'Programming Language :: Python :: 3',
                      'Topic :: Scientific/Engineering :: Astronomy',
                    ],
      url = 'http://github.com/rodluger/pysyzygy',
      author = 'Rodrigo Luger',
      author_email = 'rodluger@uw.edu',
      license = 'MIT',
      packages = [str('pysyzygy')],
      include_package_data = True,
      zip_safe = False,
      test_suite='nose.collector',
      tests_require=['nose']
      )