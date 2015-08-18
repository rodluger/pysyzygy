# TODO: Make npts a user option; make arrays dynamically allocated
# TODO: Fix zero eccentricity issues!

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
ERR_PER              =   13                                                           # Bad period
ERR_RHOS_ARS         =   14                                                           # Must specify either rhos or aRs!
ERR_RHOS             =   15                                                           # Bad rhos
ERR_ECC_W            =   16                                                           # Bad eccentricity/omega
ERR_LD               =   17                                                           # Bad limb darkening coeffs
ERR_T0               =   18                                                           # Bad t0

# Define models
QUADRATIC  =              0
KIPPING    =              1
NONLINEAR  =              2
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
MAXPTS      =             10000
TRANSITSARR =             ctypes.c_double * MAXTRANSITS
G           =             6.672e-8
DAYSEC      =             86400.

class TRANSIT(ctypes.Structure):
      _fields_ = [("bcirc", ctypes.c_double),
                  ("rhos", ctypes.c_double),
                  ("MpMs", ctypes.c_double),
                  ("esw", ctypes.c_double),
                  ("ecw", ctypes.c_double),
                  ("per", ctypes.c_double),
                  ("RpRs", ctypes.c_double),
                  ("t0", ctypes.c_double),
                  ("ecc", ctypes.c_double),
                  ("w", ctypes.c_double),
                  ("aRs", ctypes.c_double),
                  ("ntrans", ctypes.c_int),
                  ("_tN", TRANSITSARR)]
      
      def __init__(self, **kwargs):
        self.bcirc = np.nan
        self.rhos = np.nan
        self.MpMs = np.nan
        self.esw = np.nan
        self.ecw = np.nan
        self.per = np.nan
        self.RpRs = np.nan
        self.t0 = np.nan
        self.ecc = np.nan
        self.w = np.nan
        self.aRs = np.nan
        self._tN_p = []
        self._tN = TRANSITSARR(*self._tN_p)
        self.ntrans = len(self._tN_p)     
        self.update(**kwargs)
        
      def update(self, **kwargs):
        '''
        
        '''
        
        self.bcirc = kwargs.pop('bcirc', self.bcirc)
        b = kwargs.pop('b', None)
        if b is not None: 
          self.bcirc = b
        self.rhos = kwargs.pop('rhos', self.rhos)
        self.MpMs = kwargs.pop('MpMs', self.MpMs)
        self.esw = kwargs.pop('esw', self.esw)
        self.ecw = kwargs.pop('ecw', self.ecw)
        if self.esw == 0. and self.ecw == 0.:                                         # DEBUG: Fix this
          self.esw = 1.e-10
          self.ecw = 1.e-10
        self.per = kwargs.pop('per', self.per)
        self.RpRs = kwargs.pop('RpRs', self.RpRs)
        self.t0 = kwargs.pop('t0', self.t0)
        self.ecc = kwargs.pop('ecc', self.ecc)
        if self.ecc == 0.:
          self.ecc = 1.e-10
        self.w = kwargs.pop('w', self.w)
        self.aRs = kwargs.pop('aRs', self.aRs)
        tN = kwargs.pop('tN', None)
        if tN is not None:
          self._tN_p = tN                                                             # The transit times. NOTE: Must be sorted!
          self._tN = TRANSITSARR(*self._tN_p)
          self.ntrans = len(self._tN_p)                                               # Number of transits; only used if tN is set (i.e., for TTVs)
      
      @property
      def tN(self):
        return self._tN_p                                                             # The python-friendly list/array of transit times

      @tN.setter
      def tN(self, value):
        self._tN_p = value
        self.ntrans = len(self._tN_p)
        self._tN = TRANSITSARR(*self._tN_p)                                           # The pointer version that gets passed to C
      
      @property
      def tdur(self):
        '''
        The approximate transit duration for the general case of an eccentric orbit
        
        '''
        ecc = self.ecc if not np.isnan(self.ecc) else np.sqrt(self.ecw**2 + self.esw**2)
        aRs = ((G * self.rhos * (1. + self.MpMs) * 
              (self.per * DAYSEC)**2.) / (3. * np.pi))**(1./3.)
        inc = np.arccos(self.bcirc/aRs)
        becc = self.bcirc * (1 - ecc**2)/(1 - self.esw)
        tdur = self.per / 2. / np.pi * np.arcsin(((1. + self.RpRs)**2 -
               becc**2)**0.5 / (np.sin(inc) * aRs))
        tdur *= np.sqrt(1. - ecc**2.)/(1. - self.esw)
        return tdur
           
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
        self.ldmodel = QUADRATIC
        self.u1 = np.nan
        self.u2 = np.nan
        self.q1 = np.nan
        self.q2 = np.nan
        self.c1 = np.nan
        self.c2 = np.nan
        self.c3 = np.nan
        self.c4 = np.nan
        self.update(**kwargs)
        
      def update(self, **kwargs):
        '''
        
        '''
        
        self.ldmodel = kwargs.pop('ldmodel', self.ldmodel)
        self.u1 = kwargs.pop('u1', self.u1)
        self.u2 = kwargs.pop('u2', self.u1)
        self.q1 = kwargs.pop('q1', self.q1)
        self.q2 = kwargs.pop('q2', self.q2)
        self.c1 = kwargs.pop('c1', self.c1)
        self.c2 = kwargs.pop('c2', self.c2)
        self.c3 = kwargs.pop('c3', self.c3)
        self.c4 = kwargs.pop('c4', self.c4)
                  
class ARRAYS(ctypes.Structure):
      _fields_ = [("nstart", ctypes.c_int),
                  ("nend", ctypes.c_int),
                  ("ipts", ctypes.c_int),
                  ("_time", ctypes.c_double * MAXPTS),
                  ("_flux", ctypes.c_double * MAXPTS),
                  ("_bflx", ctypes.c_double * MAXPTS),
                  ("_M", ctypes.c_double * MAXPTS),
                  ("_E", ctypes.c_double * MAXPTS),
                  ("_f", ctypes.c_double * MAXPTS),
                  ("_r", ctypes.c_double * MAXPTS),
                  ("_x", ctypes.c_double * MAXPTS),
                  ("_y", ctypes.c_double * MAXPTS),
                  ("_z", ctypes.c_double * MAXPTS),
                  ("_b", ctypes.c_double * MAXPTS),
                  ("_iarr", ctypes.POINTER(ctypes.c_double))]
                  
      def __init__(self, **kwargs):        
        self._time = as_ctypes(np.zeros(MAXPTS))
        self._flux = as_ctypes(np.zeros(MAXPTS))
        self._bflx = as_ctypes(np.zeros(MAXPTS))
        self._M = as_ctypes(np.zeros(MAXPTS))
        self._E = as_ctypes(np.zeros(MAXPTS))
        self._f = as_ctypes(np.zeros(MAXPTS))
        self._r = as_ctypes(np.zeros(MAXPTS))
        self._x = as_ctypes(np.zeros(MAXPTS))
        self._y = as_ctypes(np.zeros(MAXPTS))
        self._z = as_ctypes(np.zeros(MAXPTS))
        self._b = as_ctypes(np.zeros(MAXPTS))     
        
        self.nstart = 0
        self.nend = MAXPTS
        self.ipts = 0
      
      @property
      def time(self):
        return np.array([self._time[i] for i in range(self.nstart, self.nend)])
      
      @property
      def flux(self):
        return np.array([self._flux[i] for i in range(self.nstart, self.nend)])
        
      @property
      def bflx(self):
        return np.array([self._bflx[i] for i in range(self.nstart, self.nend)])

      @property
      def M(self):
        return np.array([self._M[i] for i in range(self.nstart, self.nend)])
        
      @property
      def E(self):
        return np.array([self._E[i] for i in range(self.nstart, self.nend)])
        
      @property
      def f(self):
        return np.array([self._f[i] for i in range(self.nstart, self.nend)])
        
      @property
      def r(self):
        return np.array([self._r[i] for i in range(self.nstart, self.nend)])
      
      @property
      def x(self):
        return np.array([self._x[i] for i in range(self.nstart, self.nend)])
        
      @property
      def y(self):
        return np.array([self._y[i] for i in range(self.nstart, self.nend)])
      
      @property
      def z(self):
        return np.array([self._z[i] for i in range(self.nstart, self.nend)])
      
      @property
      def b(self):
        return np.array([self._b[i] for i in range(self.nstart, self.nend)])
      
      @property
      def iarr(self):
        return np.array([self._iarr[i] for i in range(self.ipts)])
             
class SETTINGS(ctypes.Structure):
      _fields_ = [("cadence", ctypes.c_double),
                  ("exptime", ctypes.c_double),
                  ("keptol", ctypes.c_double),
                  ("fullorbit", ctypes.c_int),
                  ("exppts", ctypes.c_int),
                  ("binmethod", ctypes.c_int),
                  ("intmethod", ctypes.c_int),
                  ("maxkepiter", ctypes.c_int),
                  ("computed", ctypes.c_int),
                  ("binned", ctypes.c_int)]
      
      def __init__(self, **kwargs):
        self.cadence = KEPLONGCAD
        self.exptime = KEPLONGEXP
        self.fullorbit = 0
        self.exppts = 50
        self.binmethod = RIEMANN
        self.intmethod = SMARTINT
        self.keptol = 1.e-15
        self.maxkepiter = 100
        self.update(**kwargs)
      
      def update(self, **kwargs):
        '''
        
        '''
        
        self.cadence = kwargs.pop('cadence', self.cadence)                            # Long cadence dt
        self.exptime = kwargs.pop('exptime', self.exptime)                            # Long cadence integration time
        self.fullorbit = 1 if kwargs.pop('fullorbit', self.fullorbit) else 0          # Compute full orbit or just the transits (default)
        self.exppts = kwargs.pop('exppts', self.exppts)                               # Average flux over this many points for binning
        self.binmethod = kwargs.pop('binmethod', self.binmethod)                      # How to integrate when binning?
        self.intmethod = kwargs.pop('intmethod', self.intmethod)                      # Integration method
        self.keptol = kwargs.pop('keptol', self.keptol)                               # Kepler solver tolerance
        self.maxkepiter = kwargs.pop('maxkepiter', self.maxkepiter)                   # Maximum number of iterations in Kepler solver
        self.computed = 0
        self.binned = 0

# Check the OS
if platform.system() == "Darwin":
  try:
    lib = ctypes.CDLL(PSZGPATH + '/pysyzygy/transit_mac.so')
  except:
    raise Exception("Can't find .so file; please type ``make`` to compile the code.")
elif platform.system() == "Linux":
  try:
    lib = ctypes.CDLL(PSZGPATH + '/pysyzygy/transit_linux.so')
  except:
    raise Exception("Can't find .so file; please type ``make`` to compile the code.")
else:
  raise Exception("Unknown platform.")

# Declare the C functions; user should access these through the Transit() class below
_Compute = lib.Compute
_Compute.restype = ctypes.c_int
_Compute.argtypes = [ctypes.POINTER(TRANSIT), ctypes.POINTER(LIMBDARK), 
                    ctypes.POINTER(SETTINGS), ctypes.POINTER(ARRAYS)]

_Bin = lib.Bin
_Bin.restype = ctypes.c_int
_Bin.argtypes = [ctypes.POINTER(TRANSIT), ctypes.POINTER(LIMBDARK), 
                ctypes.POINTER(SETTINGS), ctypes.POINTER(ARRAYS)]

_Interpolate = lib.Interpolate
_Interpolate.restype = ctypes.c_int
_Interpolate.argtypes = [ndpointer(dtype=ctypes.c_double),
                        ctypes.c_int,
                        ctypes.c_int,
                        ctypes.POINTER(TRANSIT), 
                        ctypes.POINTER(LIMBDARK), ctypes.POINTER(SETTINGS), 
                        ctypes.POINTER(ARRAYS)]

_dbl_free = lib.dbl_free
_dbl_free.argtypes = [ctypes.POINTER(ctypes.c_double)]

# Error handling
def RaiseError(err):
  if (err == ERR_NONE):
    return
  elif (err == ERR_NOT_IMPLEMENTED):
    raise Exception("Option not implemented.")
  elif (err == ERR_MAX_PTS):
    raise Exception("Maximum points in lightcurve exceeded. " + 
                    "Try decreasing `exppts`, increasing `exptime`, or recompiling " +
                    "the code with a larger value for `MAXPTS`.")  
  elif (err == ERR_NO_TRANSIT):
    raise Exception("Object does not transit the star.")  
  elif (err == ERR_BAD_ECC):
    raise Exception("Bad value for ``ecc``.")  
  elif (err == ERR_RC):
    raise Exception("Error in elliptic integral function RC().")  
  elif (err == ERR_RJ):
    raise Exception("Error in elliptic integral function RJ().") 
  elif (err == ERR_RF):
    raise Exception("Error in elliptic integral function RF().") 
  elif (err == ERR_RADIUS):
    raise Exception("Bad value for ``RpRs``.") 
  elif (err == ERR_EXP_PTS):
    raise Exception("The number of exposure points must be even.") 
  elif (err == ERR_NOT_COMPUTED):
    raise Exception("Lightcurve must be computed before it can be binned.") 
  elif (err == ERR_STAR_CROSS):
    raise Exception("Star-crossing orbit.") 
  elif (err == ERR_RHOS_ARS):
    raise Exception("Must specify one of ``rhos`` or ``aRs``.") 
  elif (err == ERR_RHOS):
    raise Exception("Bad value for ``per``.") 
  elif (err == ERR_ECC_W):
    raise Exception("Bad value for ``esw`` or ``ecw``.") 
  elif (err == ERR_LD):
    raise Exception("Bad value for the limb darkening coefficients.") 
  elif (err == ERR_T0):
    raise Exception("Bad value for ``t0``.")
  else:
    raise Excpetion("Error in transit computation.")

def AddTTVs(tN, ttv_amp = 0.1, ttv_noise = 0.02, ttv_phi = 0., ttv_per = 4.321):
  '''
  
  '''
  ttv_amp = np.atleast_1d(ttv_amp)
  ttv_per = np.atleast_1d(ttv_per)
  if len(ttv_amp) > 1 and ttv_phi == 0.:
    ttv_phi = np.zeros_like(ttv_amp)
  ttv_phi = np.atleast_1d(ttv_phi)
  
  ntrans = len(tN)
  ttv = ttv_noise * np.random.randn(ntrans - 1)                                       # Add white noise      

  for amp, per, phi in zip(ttv_amp, ttv_per, ttv_phi):
    ttv += amp * np.sin(phi + 2 * np.pi / per * np.array(range(ntrans - 1)))          # Perfectly periodic
  
  ttv = np.concatenate(([0], ttv))                                                    # First transit is not shifted
  return tN + ttv

class Transit():
  '''
  A user-friendly wrapper around the ``ctypes`` routines.
  
  '''
  
  def __init__(self, **kwargs):
          
    self.arrays = ARRAYS()
    self.limbdark = LIMBDARK()
    self.transit = TRANSIT()
    self.settings = SETTINGS()
    self.update(**kwargs)
  
  def update(self, **kwargs):
    '''
    
    '''
    
    valid = [y[0] for x in [TRANSIT, LIMBDARK, SETTINGS] for y in x._fields_]         # List of valid kwargs
    valid += ['b', 'tN']                                                              # These are special!
    for k in kwargs.keys():
      if k not in valid:
        raise Exception("Invalid kwarg '%s'." % k)  
  
    if ('q1' in kwargs.keys()) and ('q2' in kwargs.keys()):
      kwargs.update({'ldmodel': KIPPING})
    elif ('c1' in kwargs.keys()) and ('c2' in kwargs.keys()) and \
         ('c3' in kwargs.keys()) and ('c4' in kwargs.keys()):
      kwargs.update({'ldmodel': NONLINEAR})
    
    self.limbdark.update(**kwargs)
    self.transit.update(**kwargs)
    self.settings.update(**kwargs)
    
  
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
      
    err = _Interpolate(t, len(t), array, self.transit, self.limbdark, self.settings, self.arrays)
    if err != ERR_NONE: RaiseError(err)
    res = self.arrays.iarr
    self.Clean()
    return res
  
  def Compute(self):
    err = _Compute(self.transit, self.limbdark, self.settings, self.arrays)
    if err != ERR_NONE: RaiseError(err)    

  def Bin(self):
    err = _Bin(self.transit, self.limbdark, self.settings, self.arrays)
    if err != ERR_NONE: RaiseError(err)
  
  def Clean(self):
    '''
    Frees the memory used by the interpolated array.
    
    '''
    
    _dbl_free(self.arrays._iarr)
  
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
      _Compute(transit, limbdark, settings, arr)
      _Bin(transit, limbdark, settings, arr)
      t[i] = arr.time
      b[i] = arr.bflx
  
    # Compute a high-res binned lightcurve
    settings.binmethod=RIEMANN
    settings.maxpts=100000
    settings.exppts=10000
    _Compute(transit, limbdark, settings, arr)
    _Bin(transit, limbdark, settings, arr)
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
  
  def DianaTest():
    import timeit
    import __builtin__
    
    t = np.load("/Users/rod/Desktop/kep35_rod.npz")['t'][:3000]
    trn = Transit(rhos = 0.25, ecw = 0.0086125, esw = 0.1399, per = 20.7337445, t0 = 132.846716, MpMs = 0.79/0.8877, bcirc = 0.04999, q1 = 0.127, q2 = 0.999, RpRs = 0.818, exppts = 30)
    
    def compute():
      trn = Transit(rhos = 0.25, ecw = 0.0086125, esw = 0.1399, per = 20.7337445, t0 = 132.846716, MpMs = 0.79/0.8877, bcirc = 0.04999, q1 = 0.127, q2 = 0.999, RpRs = 0.818, exppts = 30)
      trn.Compute()
      trn.Bin()
    
    __builtin__.__dict__.update(locals())
    print timeit.timeit('compute()', number=100)/100
    
    pl.plot(t, trn(t, 'binned'), 'r.')

    pl.show()
  
  DianaTest()