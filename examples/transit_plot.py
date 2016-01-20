#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
transit_plot.py
---------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import pysyzygy as ps
import numpy as np
import matplotlib.pyplot as pl

# Transit model kwargs
kwargs = dict(RpRs = 0.1, b = 0.5, 
              ecc = 0.5, w = 0.,
              MpMs = 0., rhos = 1.4,
              per = 5., t0 = 0.,
              q1 = 0.45, q2 = 0.3,
              maxpts = 20000,
              exptime = ps.KEPLONGEXP)

ps.PlotTransit() #**kwargs)