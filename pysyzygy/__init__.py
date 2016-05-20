#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, unicode_literals
import os

# Check if the code needs to be compiled
if not os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'transitlib.so')):
  cwd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
  import subprocess
  subprocess.Popen(["make"], cwd = cwd)
  
from .transit import *
from .plot import *

__version__ = "0.0.1"
__author__ = "Rodrigo Luger (rodluger@uw.edu)"
__copyright__ = "Copyright 2015 Rodrigo Luger"