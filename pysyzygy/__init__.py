#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, unicode_literals
import os

# Version info
__version__ = "0.0.2"
__author__ = "Rodrigo Luger (rodluger@uw.edu)"
__copyright__ = "Copyright 2015 Rodrigo Luger"

# Check if the code needs to be compiled
if not os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'transitlib.so')):
  cwd = os.path.dirname(os.path.abspath(__file__))
  import subprocess
  subprocess.call(["make"], cwd = cwd)

# Was pysyzygy imported from setup.py?
try:
  __PYSYZYGY_SETUP__
except NameError:
  __PYSYZYGY_SETUP__ = False

# This is a regular pysyzygy run
if not __PYSYZYGY_SETUP__:
  from .transit import *
  from .plot import *
