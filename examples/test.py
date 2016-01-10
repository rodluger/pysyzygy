#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test.py
-------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import pysyzygy as ps
import numpy as np
import matplotlib.pyplot as pl


trn = ps.Transit()
time = np.linspace(-1,1,1000)
pl.plot(time, trn(time))
pl.show()
