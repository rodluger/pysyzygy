#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
import os
PSZGPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))                # This is the path to the top-level pysyzygy folder
from .plot import *
from .transit import Transit

__version__ = "0.0.1"
__author__ = "Rodrigo Luger (rodluger@uw.edu)"
__copyright__ = "Copyright 2015 Rodrigo Luger"

