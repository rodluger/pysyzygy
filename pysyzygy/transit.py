#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
transit.py
----------

A ``ctypes`` wrapper around a generalized C implementation of the 
Mandel & Agol (2002) transit model.

'''

import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer, as_ctypes
import platform
from pysyzygy import PSZGPATH

# Define errors
ERR_NONE             =   0                                                            # We're good!
ERR_NOT_IMPLEMENTED  =   1                                                            # Function/option not yet implemented
ERR_MAX_PTS          =   2                                                            # Maximum number of points exceeded in transit. Increase settings.maxpts.
ERR_KEPLER           =   3                                                            # Error in the Kepler solver; probably didn't converge
ERR_NO_TRANSIT       =   4                                                            # The planet doesn't transit the star
ERR_BAD_ECC          =   5                                                            # Bad value for eccentricity
ERR_RC               =   6                                                            # Error in rc() function
ERR_RJ               =   7                                                            # Error in rj() function
ERR_RF               =   8                                                            # Error in rf() function
ERR_RADIUS           =   9                                                            # Bad input radius
ERR_EXP_PTS          =   10                                                           # The number of exposure points cannot be odd
ERR_NOT_COMPUTED     =   11                                                           # User attempted to bin before computing
ERR_STAR_CROSS       =   12                                                           # Star-crossing orbit

# Define models
QUADRATIC  =              0
KIPPING    =              1
NONLINEAR  =              2
ECCENTRIC  =              3
CIRCULAR   =              4
RIEMANN    =              5
TRAPEZOID  =              6
SMARTINT   =              7
SLOWINT    =              8

# Cadences
KEPLONGEXP =              (1765.5/86400.)
KEPLONGCAD =              (1800./86400.)
KEPSHRTEXP =              (58.89/86400.)
KEPSHRTCAD =              (60./86400.)

# Array IDs
ARR_FLUX    =             0
ARR_BFLX    =             1
ARR_M       =             2
ARR_E       =             3
ARR_F       =             4
ARR_R       =             5
ARR_X       =             6
ARR_Y       =             7
ARR_Z       =             8
ARR_B       =             9

# Other
MAXTRANSITS =             500
TRANSITSARR =             ctypes.c_double * MAXTRANSITS

class TRANSIT(ctypes.Structure):
      _fields_ = [("model", ctypes.c_int),
                  ("bcirc", ctypes.c_double),
                  ("rhos", ctypes.c_double),
                  ("MpMs", ctypes.c_double),
                  ("esw", ctypes.c_double),
                  ("ecw", ctypes.c_double),
                  ("per", ctypes.c_double),
                  ("RpRs", ctypes.c_double),
                  ("t0", ctypes.c_double),
                  ("ntrans", ctypes.c_int),
                  ("_tN", TRANSITSARR)]
      
      def __init__(self, **kwargs):
        self.model = kwargs.pop('model', ECCENTRIC)
        self.bcirc = kwargs.pop('bcirc', 0.)
        self.rhos = kwargs.pop('rhos', 1.)
        self.MpMs = kwargs.pop('MpMs', 0.)
        self.esw = kwargs.pop('esw', 0.)
        self.ecw = kwargs.pop('ecw', 0.)
        self.per = kwargs.pop('per', 1.)
        self.RpRs = kwargs.pop('RpRs', 0.1)
        self.t0 = kwargs.pop('t0', 0.)
        self._tN_p = kwargs.pop('tN', [])                                             # The transit times. NOTE: Must be sorted!
        self._tN = TRANSITSARR(*self._tN_p)
        self.ntrans = len(self._tN_p)                                                 # Number of transits; only used if tN is set (i.e., for TTVs)
      
      @property
      def tN(self):
        return self._tN_p                                                             # The python-friendly list/array of transit times

      @tN.setter
      def tN(self, value):
        self._tN_p = value
        self.ntrans = len(self._tN_p)
        self._tN = TRANSITSARR(*self._tN_p)                                           # The pointer version that gets passed to C
      
      
class LIMBDARK(ctypes.Structure):
      _fields_ = [("ldmodel", ctypes.c_int),
                  ("u1", ctypes.c_double),
                  ("u2", ctypes.c_double),  
                  ("q1", ctypes.c_double),
                  ("q2", ctypes.c_double),  
                  ("c1", ctypes.c_double),
                  ("c2", ctypes.c_double),  
                  ("c3", ctypes.c_double),
                  ("c4", ctypes.c_double)]
                  
      def __init__(self, **kwargs):
        self.ldmodel = kwargs.pop('ldmodel', QUADRATIC)
        self.u1 = kwargs.pop('u1', 1.)
        self.u2 = kwargs.pop('u2', 0.)
        self.q1 = kwargs.pop('q1', 0.)
        self.q2 = kwargs.pop('q2', 0.)
        self.c1 = kwargs.pop('c1', 0.)
        self.c2 = kwargs.pop('c2', 0.)
        self.c3 = kwargs.pop('c3', 0.)
        self.c4 = kwargs.pop('c4', 0.)
                  
class ARRAYS(ctypes.Structure):
      _fields_ = [("npts", ctypes.c_int),
                  ("ipts", ctypes.c_int),
                  ("_time", ctypes.POINTER(ctypes.c_double)),
                  ("_flux", ctypes.POINTER(ctypes.c_double)),
                  ("_bflx", ctypes.POINTER(ctypes.c_double)),
                  ("_M", ctypes.POINTER(ctypes.c_double)),
                  ("_E", ctypes.POINTER(ctypes.c_double)),
                  ("_f", ctypes.POINTER(ctypes.c_double)),
                  ("_r", ctypes.POINTER(ctypes.c_double)),
                  ("_x", ctypes.POINTER(ctypes.c_double)),
                  ("_y", ctypes.POINTER(ctypes.c_double)),
                  ("_z", ctypes.POINTER(ctypes.c_double)),
                  ("_b", ctypes.POINTER(ctypes.c_double)),
                  ("_iarr", ctypes.POINTER(ctypes.c_double))]
                  
      def __init__(self, **kwargs):
        # self.M = as_ctypes(np.zeros(ndata))
        self.npts = 0
        self.ipts = 0
      
      @property
      def time(self):
        return np.array([self._time[i] for i in range(self.npts)])
      
      @property
      def flux(self):
        return np.array([self._flux[i] for i in range(self.npts)])
        
      @property
      def bflx(self):
        return np.array([self._bflx[i] for i in range(self.npts)])

      @property
      def M(self):
        return np.array([self._M[i] for i in range(self.npts)])
        
      @property
      def E(self):
        return np.array([self._E[i] for i in range(self.npts)])
        
      @property
      def f(self):
        return np.array([self._f[i] for i in range(self.npts)])
        
      @property
      def r(self):
        return np.array([self._r[i] for i in range(self.npts)])
      
      @property
      def x(self):
        return np.array([self._x[i] for i in range(self.npts)])
        
      @property
      def y(self):
        return np.array([self._y[i] for i in range(self.npts)])
      
      @property
      def z(self):
        return np.array([self._z[i] for i in range(self.npts)])
      
      @property
      def b(self):
        return np.array([self._b[i] for i in range(self.npts)])
      
      @property
      def iarr(self):
        return np.array([self._iarr[i] for i in range(self.ipts)])
             
class SETTINGS(ctypes.Structure):
      _fields_ = [("cadence", ctypes.c_double),
                  ("exptime", ctypes.c_double),
                  ("keptol", ctypes.c_double),
                  ("fullorbit", ctypes.c_int),
                  ("maxpts", ctypes.c_int),
                  ("exppts", ctypes.c_int),
                  ("binmethod", ctypes.c_int),
                  ("intmethod", ctypes.c_int),
                  ("maxkepiter", ctypes.c_int),
                  ("computed", ctypes.c_int),
                  ("binned", ctypes.c_int)]
      
      def __init__(self, **kwargs):
        self.cadence = kwargs.pop('cadence', KEPLONGCAD)                              # Long cadence dt
        self.exptime = kwargs.pop('exptime', KEPLONGEXP)                              # Long cadence integration time
        fullorbit = kwargs.pop('fullorbit', False)                                    # Compute full orbit or just the transits (default)
        if fullorbit: self.fullorbit = 1
        else: self.fullorbit = 0
        self.exppts = kwargs.pop('exppts', 50)                                        # Average flux over 10 points for binning
        self.maxpts = kwargs.pop('maxpts', 10000)                                     # Maximum length of arrays ( > exp_pts * transit duration / exptime )
        self.binmethod = kwargs.pop('binmethod', RIEMANN)                             # How to integrate when binning?
        self.intmethod = kwargs.pop('intmethod', SMARTINT)                            # Integration method
        self.keptol = kwargs.pop('keptol', 1.e-15)                                    # Kepler solver tolerance
        self.maxkepiter = kwargs.pop('maxkepiter', 100)                               # Maximum number of iterations in Kepler solver
        self.computed = 0
        self.binned = 0

# Check the OS
if platform.system() == "Darwin":
  lib = ctypes.CDLL(PSZGPATH + '/pysyzygy/transit_mac.so')
elif platform.system() == "Linux":
  lib = ctypes.CDLL(PSZGPATH + '/pysyzygy/transit_linux.so')
else:
  raise Exception("Unknown platform.")

# Declare the C functions
Compute = lib.Compute
Compute.restype = ctypes.c_int
Compute.argtypes = [ctypes.POINTER(TRANSIT), ctypes.POINTER(LIMBDARK), 
                    ctypes.POINTER(SETTINGS), ctypes.POINTER(ARRAYS)]

Bin = lib.Bin
Bin.restype = ctypes.c_int
Bin.argtypes = [ctypes.POINTER(TRANSIT), ctypes.POINTER(LIMBDARK), 
                ctypes.POINTER(SETTINGS), ctypes.POINTER(ARRAYS)]

_Interpolate = lib.Interpolate
_Interpolate.restype = ctypes.c_int
_Interpolate.argtypes = [ndpointer(dtype=ctypes.c_double),
                        ctypes.c_int,
                        ctypes.c_int,
                        ctypes.POINTER(TRANSIT), 
                        ctypes.POINTER(LIMBDARK), ctypes.POINTER(SETTINGS), 
                        ctypes.POINTER(ARRAYS)]

# Error handling
def RaiseError(err):
  if (err == ERR_NONE):
    return
  elif (err == ERR_NOT_IMPLEMENTED):
    raise Exception("Option not implemented.")
  elif (err == ERR_MAX_PTS):
    raise Exception("Maximum points in lightcurve exceeded.")  
  elif (err == ERR_NO_TRANSIT):
    raise Exception("Object does not transit the star.")  
  elif (err == ERR_BAD_ECC):
    raise Exception("Bad value for the eccentricity.")  
  elif (err == ERR_RC):
    raise Exception("Error in elliptic integral function RC().")  
  elif (err == ERR_RJ):
    raise Exception("Error in elliptic integral function RJ().") 
  elif (err == ERR_RF):
    raise Exception("Error in elliptic integral function RF().") 
  elif (err == ERR_RADIUS):
    raise Exception("Bad value for radius.") 
  elif (err == EXP_PTS):
    raise Exception("The number of exposure points must be even.") 
  elif (err == ERR_NOT_COMPUTED):
    raise Exception("Lightcurve must be computed before it can be binned.") 
  elif (err == ERR_STAR_CROSS):
    raise Exception("Star-crossing orbit.") 
  else:
    raise Excpetion("Error in transit computation.")

# User-friendly wrapper
class Transit():
  '''
  
  '''
  
  def __init__(self, **kwargs):
    self._arrays = ARRAYS(**kwargs)
    self._ldark = LIMBDARK(**kwargs)
    self._trans = TRANSIT(**kwargs)
    self._stngs = SETTINGS(**kwargs)
    
    '''
    if kwargs != {}:
      raise Exception("Unknown kwarg '%s'" % kwargs.keys()[0])
    '''
    
  def __call__(self, t, param):
    if param == 'flux':
      array = ARR_FLUX
    elif param == 'binned':
      array = ARR_BFLX
    elif param == 'M':
      array = ARR_M
    elif param == 'E':
      array = ARR_E
    elif param == 'f':
      array = ARR_F
    elif param == 'r':
      array = ARR_R
    elif param == 'x':
      array = ARR_X
    elif param == 'y':
      array = ARR_Y
    elif param == 'z':
      array = ARR_Z
    elif param == 'b':
      array = ARR_B
    else:
      RaiseError(ERR_NOT_IMPLEMENTED)
      
    err = _Interpolate(t, len(t), array, self._trans, self._ldark, self._stngs, self._arrays)
    if err != ERR_NONE: RaiseError(err)
    return self._arrays.iarr
    
if __name__ == '__main__':
  '''
  For debugging only
  
  '''
  import matplotlib.pyplot as pl
  from scipy import interpolate
  
  def PlotErrors():
    arr = ARRAYS()
    limbdark = LIMBDARK()
    transit = TRANSIT(ecw = 0.01, esw = 0.01, bcirc = 0.5, RpRs = 0.1, per = 5.0)
    settings = SETTINGS()
  
    # Compute low-res binned lightcurves with different settings
    t = [[]]*6
    b = [[]]*6
    for i, bm, ep in zip(range(6), [RIEMANN, TRAPEZOID]*3, [10, 10, 20, 20, 50, 50]):
      settings.binmethod=bm
      settings.exppts=ep
      Compute(transit, limbdark, settings, arr)
      Bin(transit, limbdark, settings, arr)
      t[i] = arr.time
      b[i] = arr.bflx
  
    # Compute a high-res binned lightcurve
    settings.binmethod=RIEMANN
    settings.maxpts=100000
    settings.exppts=10000
    Compute(transit, limbdark, settings, arr)
    Bin(transit, limbdark, settings, arr)
    B = interpolate.interp1d(arr.time, arr.bflx, bounds_error=False)
  
    for i, style, label in zip(range(6), ['r-', 'b-', 'r--', 'b--', 'r-.', 'b-.'], 
                                         ['R10', 'T10', 'R30', 'T30', 'R50', 'T50']):
      pl.plot(t[i], np.abs(b[i] - B(t[i])), style, label = label)
    pl.legend(loc='upper left')
    pl.yscale('log')
    pl.xlabel('Time (days)', fontsize=18)
    pl.ylabel('abs(error)', fontsize=18)
    pl.show()
  
  def PlotInterpolation():  
    t = np.arange(-2.5,2.5,KEPLONGCAD)
    t += 0.001 * np.random.randn(len(t))

    trn = Transit(ecw = 0.1, esw = 0.1, bcirc = 0.3, per = 1.3, fullorbit = True)
    pl.plot(t, trn(t, 'x'), 'r.')
    pl.show()
  
  PlotInterpolation()