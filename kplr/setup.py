#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup
    setup

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

# Hackishly inject a constant into builtins to enable importing of the
# module.
if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins
builtins.__KPLR_SETUP__ = True
import kplr

setup(
    name="kplr",
    version=kplr.__version__,
    author="Daniel Foreman-Mackey",
    author_email="danfm@nyu.edu",
    packages=["kplr"],
    url="https://github.com/dfm/kplr",
    license="MIT",
    description="Tools for working with Kepler data in Python",
    long_description=open("README.rst").read(),
    install_requires=["six"],
    classifiers=[
        # "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
    ],
)
