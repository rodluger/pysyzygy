#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
grid.py
-------

Produce a grid of planet orbits, only one of which is a transit. If ``emcee``
is installed, attempts to use ``emcee.MPIPool()`` for multiprocessing.

Once all the PNGs have been generated, run

>>> ffmpeg -framerate 100 -i %04d.png -vcodec mpeg4 -b:v 10000k video.mp4

to produce an MPEG animation.

'''

import matplotlib as mpl; mpl.use('Agg')
import pysyzygy as ps
import matplotlib.pyplot as pl
import numpy as np
import sys

# Define some stuff
path = ''
kw = [{'RpRs': 0.5,  'ecc': 0.1, 'w': 1.0,  
       'image_map': 'earth', 'starcolor': (1.0, 0.85, 0.1), 'b': 2.0},
      {'RpRs': 0.25, 'ecc': 0.2, 'w': 0.,   
       'image_map': 'jupiter', 'starcolor': (1.0, 0.75, 0.1),  'b': 1.9},
      {'RpRs': 0.5,  'ecc': 0.,  'w': 0.,   
       'image_map': 'neptune', 'starcolor': (1.0, 0.65, 0.2), 'b': 2.0},
      {'RpRs': 0.35, 'ecc': 0.3, 'w': 4.0,  
       'image_map': 'earth', 'starcolor': (1.0, 0.35, 0.1), 'b': 1.5},
      {'RpRs': 0.5,  'ecc': 0.3, 'w': 2.0,  
       'image_map': 'earth', 'starcolor': (1.0, 0.95, 0.0), 'b': 2.0},
      {'RpRs': 0.45, 'ecc': 0.1, 'w': 1.5,   
       'image_map': 'neptune', 'starcolor': (1.0, 0.85, 0.2), 'b': 1.8},
      {'RpRs': 0.5,  'ecc': 0.2, 'w': 1.9,  
       'image_map': 'jupiter', 'starcolor': (1.0, 0.65, 0.1),  'b': 1.9},
      {'RpRs': 0.45, 'ecc': 0.,  'w': 0,    
       'image_map': 'neptune', 'starcolor': (0.9, 0.85, 0.1), 'b': 0.2},
      {'RpRs': 0.1,  'ecc': 0.,  'w': 0,    
       'image_map': 'neptune', 'starcolor': (1.0, 0.35, 0.1), 'b': 1.5}]
M0 = [1.3, 2.0, 0.0, 1.0, 4.0, 1.5, 5.8, 0.5, 3.3]
P = [1./4, 1./2, 1./4, 1./3, 1./4, 1./3, 1./4, 1., 1./2]

def plot(args):
  '''
  Plot a grid of planet orbits.
  
  '''
  f, time, long0 = args
  fig, ax = pl.subplots(3, 3)
  fig.set_size_inches(12,9)
  fig.subplots_adjust(hspace = 0.05, wspace = 0.05)
  ax = [x for y in ax for x in y]
  for axis, kwargs, m, p in zip(ax, kw, M0, P):
    M = (((time/p % 1.) * 2 * np.pi) + m) % (2 * np.pi)
    ps.PlotImage(per = 0.5, u1 = 1., u2 = 0., bkgimage = 'stars', ax = axis, 
                 long0 = long0, rhos = 1.0, xlims = (-5,5), ylims = (-3,3), 
                 trail = True, M = M, **kwargs)
  fig.savefig(path + '%04d.png' % f, facecolor = 'black', bbox_inches = 'tight')
  pl.close()

# Set up parallel processing
try:
  from emcee.utils import MPIPool
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

# Plot: A grid of 9 planet orbits, only one of which is a transit
nsteps = 2000
dpy = 4
frames = range(nsteps)
time = np.linspace(0,1,nsteps,endpoint=False)
rotation = (np.linspace(0, -dpy, nsteps) % 1)*360 - 180
input = zip(frames, time, rotation)
MAP(plot, input)

if mpi:
  pool.close()