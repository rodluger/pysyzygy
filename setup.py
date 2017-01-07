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

long_description = \
"""
A fast and general planet transit (syzygy) code written in C and in Python.
Pysyzygy computes fast lightcurves for the most general case of a massive, 
eccentric planet orbiting a limb-darkened star. The code is based
on the Mandel & Agol (2002) transit model.
"""

# Setup!
setup(name = 'pysyzygy',
      version = pysyzygy.__version__,
      description = 'Transit modeling in Python',
      long_description = long_description,
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
      keywords = 'exoplanets transits',
      packages = [str('pysyzygy')],
      include_package_data = True,
      install_requires = [
                          'numpy>=1.8',
                          'scipy',
                          'matplotlib',
                          'ctypes'
                         ],
      zip_safe = False,
      test_suite='nose.collector',
      tests_require=['nose']
      )