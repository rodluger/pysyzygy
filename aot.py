#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
aot.py
------

Produce plots for my AoT talk on July 22, 2015. Tailored to Hyak.

'''

import matplotlib as mpl; mpl.use('Agg')
import pysyzygy as ps
import matplotlib.pyplot as pl
import numpy as np
import sys
from emcee.utils import MPIPool

# Define some stuff
path = '/gscratch/vsm/rodluger/aot/'
kw = [{'per': 0.50, 'RpRs': 0.5,  'ecc': 0.1, 'w': 1.0,  
       'image_map': 'earth', 'starcolor': (1.0, 0.85, 0.1), 'b': 2.0},
      {'per': 0.50, 'RpRs': 0.25, 'ecc': 0.2, 'w': 0.,   
       'image_map': 'jupiter', 'starcolor': (1.0, 0.75, 0.1),  'b': 1.9},
      {'per': 0.50, 'RpRs': 0.5,  'ecc': 0.,  'w': 0.,   
       'image_map': 'neptune', 'starcolor': (1.0, 0.65, 0.2), 'b': 2.0},
      {'per': 0.50, 'RpRs': 0.35, 'ecc': 0.3, 'w': 4.0,  
       'image_map': 'earth', 'starcolor': (1.0, 0.35, 0.1), 'b': 1.5},
      {'per': 0.50, 'RpRs': 0.5,  'ecc': 0.3, 'w': 2.0,  
       'image_map': 'earth', 'starcolor': (1.0, 0.95, 0.0), 'b': 2.0},
      {'per': 0.50,'RpRs': 0.45, 'ecc': 0.1, 'w': 1.5,   
       'image_map': 'neptune', 'starcolor': (1.0, 0.85, 0.2), 'b': 1.8},
      {'per': 0.50, 'RpRs': 0.5,  'ecc': 0.2, 'w': 1.9,  
       'image_map': 'jupiter', 'starcolor': (1.0, 0.65, 0.1),  'b': 1.9},
      {'per': 0.50, 'RpRs': 0.45, 'ecc': 0.,  'w': 0,    
       'image_map': 'neptune', 'starcolor': (0.9, 0.85, 0.1), 'b': 0.2},
      {'per': 0.50, 'RpRs': 0.1,  'ecc': 0.,  'w': 0,    
       'image_map': 'jupiter', 'starcolor': (0.9, 0.95, 0.1),'b': 1.5}]
M0 = [1.3, 2.0, 0.0, 1.0, 4.0, 1.5, 5.8, 0.5, 3.3]

def plot(args):
  '''
  Plot a grid of planet orbits.
  
  '''
  f, M, long0 = args
  fig, ax = pl.subplots(3, 3)
  fig.set_size_inches(12,9)
  fig.subplots_adjust(hspace = 0.05, wspace = 0.05)
  ax = [x for y in ax for x in y]
  for axis, kwargs, m in zip(ax, kw, M0):
    ps.PlotImage(u1 = 1., u2 = 0., bkgimage = 'stars', ax = axis, long0 = long0,
                 rhos = 1.0, xlims = (-5,5), ylims = (-3,3), trail = True, 
                 M = (M + m) % (2*np.pi), **kwargs)
  fig.savefig(path + '%03d.png' % f, facecolor = 'black', bbox_inches = 'tight')
  pl.close()

# Set up parallel processing
try:
  pool = MPIPool()
  MAP = pool.map
  mpi = True
except:
  pool = None
  MAP = map
  mpi = False
if mpi:  
  if not pool.is_master():
    pool.wait()                                                                       
    sys.exit(0)

# PLOT 1: A grid of 9 planet orbits, only one of which is a transit
nsteps = 500
dpy = 4
frames = range(nsteps)
M = np.linspace(0,2*np.pi,nsteps,endpoint=False)
rotation = (np.linspace(0, -dpy, nsteps) % 1)*360 - 180
input = zip(frames, M, rotation)
MAP(plot, input)
if mpi:
  pool.close()


'''
# Main animation
ps.AnimateImage(per = 0.5, RpRs = 0.5, ecc = 0, rhos = 1.0,
               b = 0.8, u1 = 1., u2 = 0., delay = 20,
               bkgimage = 'stars', nsteps = 10,
               image_map = 'earth', size_inches = (12,9),
               lightcurve = True)
'''