# MAC:
# >>> gcc -fPIC -c transit.c
# >>> gcc -shared -Wl,-install_name,transit_mac.so -o transit_mac.so transit.o -lc
# >>> rm transit.o
#
# LINUX:
# >>> gcc -fPIC -Wl,-Bsymbolic-functions -c -O3 transit.c
# >>> gcc -shared -O3 -Wl,-Bsymbolic-functions,-soname,transit_linux.so -o transit_linux.so transit.o -lc
# >>> rm transit.o

import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer, as_ctypes
import platform

# Define models
QUADRATIC =              0
KIPPING   =              1
NONLINEAR =              2
ECCENTRIC =              3
CIRCULAR  =              4
RIEMANN   =              5
TRAPEZOID =              6

class TRANSIT(ctypes.Structure):
      _fields_ = [("model", ctypes.c_int),
                  ("bcirc", ctypes.c_double),
                  ("rhos", ctypes.c_double),
                  ("MpMs", ctypes.c_double),
                  ("esw", ctypes.c_double),
                  ("ecw", ctypes.c_double),
                  ("per", ctypes.c_double),
                  ("RpRs", ctypes.c_double)]
      
      def __init__(self, **kwargs):
        self.model = kwargs.get('model', ECCENTRIC)
        self.bcirc = kwargs.get('bcirc', 0.)
        self.rhos = kwargs.get('rhos', 1.)
        self.MpMs = kwargs.get('MpMs', 0.)
        self.esw = kwargs.get('esw', 0.)
        self.ecw = kwargs.get('ecw', 0.)
        self.per = kwargs.get('per', 1.)
        self.RpRs = kwargs.get('RpRs', 0.1)
      
class LIMBDARK(ctypes.Structure):
      _fields_ = [("model", ctypes.c_int),
                  ("u1", ctypes.c_double),
                  ("u2", ctypes.c_double),  
                  ("q1", ctypes.c_double),
                  ("q2", ctypes.c_double),  
                  ("c1", ctypes.c_double),
                  ("c2", ctypes.c_double),  
                  ("c3", ctypes.c_double),
                  ("c4", ctypes.c_double)]
                  
      def __init__(self, **kwargs):
        self.model = kwargs.get('model', QUADRATIC)
        self.u1 = kwargs.get('u1', 1.)
        self.u2 = kwargs.get('u2', 0.)
        self.q1 = kwargs.get('q1', 0.)
        self.q2 = kwargs.get('q2', 0.)
        self.c1 = kwargs.get('c1', 0.)
        self.c2 = kwargs.get('c2', 0.)
        self.c3 = kwargs.get('c3', 0.)
        self.c4 = kwargs.get('c4', 0.)
                  
class ARRAYS(ctypes.Structure):
      _fields_ = [("npts", ctypes.c_int),
                  ("ipts", ctypes.c_int),
                  ("_time", ctypes.POINTER(ctypes.c_double)),
                  ("_flux", ctypes.POINTER(ctypes.c_double)),
                  ("_bflx", ctypes.POINTER(ctypes.c_double)),
                  ("_iflx", ctypes.POINTER(ctypes.c_double)),
                  ("_M", ctypes.POINTER(ctypes.c_double)),
                  ("_E", ctypes.POINTER(ctypes.c_double)),
                  ("_f", ctypes.POINTER(ctypes.c_double)),
                  ("_r", ctypes.POINTER(ctypes.c_double)),
                  ("_x", ctypes.POINTER(ctypes.c_double)),
                  ("_y", ctypes.POINTER(ctypes.c_double)),
                  ("_z", ctypes.POINTER(ctypes.c_double)),
                  ("_b", ctypes.POINTER(ctypes.c_double))]
                  
      def __init__(self):
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
      def iflx(self):
        return np.array([self._iflx[i] for i in range(self.ipts)])
        
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
       
class SETTINGS(ctypes.Structure):
      _fields_ = [("exptime", ctypes.c_double),
                  ("keptol", ctypes.c_double),
                  ("maxpts", ctypes.c_int),
                  ("exppts", ctypes.c_int),
                  ("binmethod", ctypes.c_int),
                  ("maxkepiter", ctypes.c_int),
                  ("computed", ctypes.c_int),
                  ("binned", ctypes.c_int)]
      
      def __init__(self, **kwargs):
        self.exptime = kwargs.get('exptime', 1765.5/86400)                            # Long cadence
        self.maxpts = kwargs.get('maxpts', 10000)                                     # Maximum length of arrays
        self.exppts = kwargs.get('exppts', 10)                                        # Average flux over 10 points for binning
        self.binmethod = kwargs.get('binmethod', RIEMANN)                             # How to integrate when binning?
        self.keptol = kwargs.get('keptol', 1.e-15)                                    # Kepler solver tolerance
        self.maxkepiter = kwargs.get('maxkepiter', 100)                               # Maximum number of iterations in Kepler solver
        self.computed = 0
        self.binned = 0
        
if platform.system() == "Darwin":
  lib = ctypes.CDLL('transit_mac.so')
elif platform.system() == "Linux":
  lib = ctypes.CDLL('transit_linux.so')
else:
  raise Exception("Unknown platform.")

Compute = lib.Compute
Compute.restype = ctypes.c_int
Compute.argtypes = [ctypes.POINTER(TRANSIT), ctypes.POINTER(LIMBDARK), 
                    ctypes.POINTER(SETTINGS), ctypes.POINTER(ARRAYS)]

Bin = lib.Bin
Bin.restype = ctypes.c_int
Bin.argtypes = [ctypes.POINTER(TRANSIT), ctypes.POINTER(LIMBDARK), 
                ctypes.POINTER(SETTINGS), ctypes.POINTER(ARRAYS)]

Interpolate = lib.Interpolate
Interpolate.restype = ctypes.c_int
Interpolate.argtypes = [ndpointer(dtype=ctypes.c_double),
                        ctypes.c_int,
                        ctypes.POINTER(TRANSIT), 
                        ctypes.POINTER(LIMBDARK), ctypes.POINTER(SETTINGS), 
                        ctypes.POINTER(ARRAYS)]

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
  
    for i, style, label in zip(range(6), ['r-', 'b-', 'r--', 'b--', 'r-.', 'b-.'], ['R10', 'T10', 'R30', 'T30', 'R50', 'T50']):
      pl.plot(t[i], np.abs(b[i] - B(t[i])), style, label = label)
    pl.legend(loc='upper left')
    pl.yscale('log')
    pl.xlabel('Time (days)', fontsize=18)
    pl.ylabel('abs(error)', fontsize=18)
    pl.show()
  
  def PlotInterpolation():  
    arr = ARRAYS()
    limbdark = LIMBDARK()
    transit = TRANSIT(ecw = 0.01, esw = 0.01, bcirc = 0.5, RpRs = 0.1, per = 5.0)
    settings = SETTINGS()
  
    ipts = 100
    t = np.linspace(-0.1,0.1,ipts)
    Interpolate(t, ipts, transit, limbdark, settings, arr)
    pl.plot(t, arr.iflx, 'b.')
  
    pl.show()
  
  PlotInterpolation()