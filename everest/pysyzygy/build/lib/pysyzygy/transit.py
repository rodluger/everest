#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`transit.py` - Transit models
-------------------------------------

A :py:mod:`ctypes` wrapper around a generalized C implementation of the 
`Mandel and Agol (2002) <http://adsabs.harvard.edu/abs/2002ApJ...580L.171M>`_
transit model

.. warning:: The longitude of pericenter `w` may be defined differently here than in other transit codes; \
             watch out for a possible offset of :math:`\pi` from what you're used to.

.. todo::
   - Add nonlinear limb darkening
   - Add secondary eclipses
   
'''

from __future__ import division, print_function, absolute_import, unicode_literals
import ctypes
import numpy as np
import os
from numpy.ctypeslib import ndpointer, as_ctypes

# Define errors
_ERR_NONE             =   0                                                           # We're good!
_ERR_NOT_IMPLEMENTED  =   1                                                           # Function/option not yet implemented
_ERR_MAX_PTS          =   2                                                           # Maximum number of points exceeded in transit. Increase settings.maxpts.
_ERR_KEPLER           =   3                                                           # Error in the Kepler solver; probably didn't converge
_ERR_NO_TRANSIT       =   4                                                           # The planet doesn't transit the star
_ERR_BAD_ECC          =   5                                                           # Bad value for eccentricity
_ERR_RC               =   6                                                           # Error in rc() function
_ERR_RJ               =   7                                                           # Error in rj() function
_ERR_RF               =   8                                                           # Error in rf() function
_ERR_RADIUS           =   9                                                           # Bad input radius
_ERR_EXP_PTS          =   10                                                          # The number of exposure points cannot be odd
_ERR_NOT_COMPUTED     =   11                                                          # User attempted to bin before computing
_ERR_STAR_CROSS       =   12                                                          # Star-crossing orbit
_ERR_PER              =   13                                                          # Bad period
_ERR_RHOS_ARS         =   14                                                          # Must specify either rhos or aRs!
_ERR_RHOS             =   15                                                          # Bad rhos
_ERR_ECC_W            =   16                                                          # Bad eccentricity/omega
_ERR_LD               =   17                                                          # Bad limb darkening coeffs
_ERR_T0               =   18                                                          # Bad t0

# Define models
QUADRATIC  =              0
KIPPING    =              1
NONLINEAR  =              2
RIEMANN    =              5
TRAPEZOID  =              6
SMARTINT   =              7
SLOWINT    =              8
MDFAST     =              9
NEWTON     =              10

# Cadences
KEPLONGEXP =              (1765.5/86400.)
KEPLONGCAD =              (1800./86400.)
KEPSHRTEXP =              (58.89/86400.)
KEPSHRTCAD =              (60./86400.)

# Array IDs
_ARR_FLUX    =             0
_ARR_BFLX    =             1
_ARR_M       =             2
_ARR_E       =             3
_ARR_F       =             4
_ARR_R       =             5
_ARR_X       =             6
_ARR_Y       =             7
_ARR_Z       =             8
_ARR_B       =             9

# Other
MAXTRANSITS =             500
TRANSITSARR =             ctypes.c_double * MAXTRANSITS
G           =             6.672e-8
DAYSEC      =             86400.

class TRANSIT(ctypes.Structure):
      '''
      The class containing all the transiting planet parameters
      
      '''
      
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
        self._tN_p = []
        self._tN = TRANSITSARR(*self._tN_p)
        self.ntrans = len(self._tN_p)     
        self.update(**kwargs)
        
      def update(self, **kwargs):
        '''
        
        '''
        
        self.MpMs = kwargs.pop('MpMs', 0.)
        self.per = kwargs.pop('per', 10.)
        self.RpRs = kwargs.pop('RpRs', 0.1)
        
        self.bcirc = kwargs.pop('bcirc', 0.)
        b = kwargs.pop('b', None)                                                     # User may specify either ``b`` or ``bcirc``
        if b is not None: 
          self.bcirc = b
        
        self.rhos = kwargs.pop('rhos', 1.4)                                           # User may specify either ``rhos`` or ``aRs``
        self.aRs = kwargs.pop('aRs', np.nan)
        if not np.isnan(self.aRs):
          self.rhos = np.nan
               
        self.ecc = kwargs.pop('ecc', 0.)                                              # User may specify (``esw`` and ``ecw``) or (``ecc`` and ``w``)
        self.w = kwargs.pop('w', 0.)
        self.esw = kwargs.pop('esw', np.nan)
        self.ecw = kwargs.pop('ecw', np.nan)
        if (not np.isnan(self.esw)) and (not np.isnan(self.ecw)):
          self.ecc = np.nan
          self.w = np.nan
                
        self.t0 = kwargs.pop('t0', 0.)                                                # User may specify either ``t0`` or ``times``
        tN = kwargs.pop('times', None)
        if tN is not None:
          self.t0 = np.nan
          self._tN_p = tN                                                             # The transit times. NOTE: Must be sorted!
          self._tN = TRANSITSARR(*self._tN_p)
          self.ntrans = len(self._tN_p)                                               # Number of transits; only used if tN is set (i.e., for TTVs)
      
      @property
      def times(self):
        return self._tN_p                                                             # The python-friendly list/array of transit times

      @times.setter
      def times(self, value):
        self._tN_p = value
        self.ntrans = len(self._tN_p)
        self._tN = TRANSITSARR(*self._tN_p)                                           # The pointer version that gets passed to C
      
      @property
      def duration(self):
        '''
        The approximate transit duration for the general case of an eccentric orbit
        
        '''
        ecc = self.ecc if not np.isnan(self.ecc) else np.sqrt(self.ecw**2 + self.esw**2)
        esw = self.esw if not np.isnan(self.esw) else ecc * np.sin(self.w)
        aRs = ((G * self.rhos * (1. + self.MpMs) * 
              (self.per * DAYSEC)**2.) / (3. * np.pi))**(1./3.)
        inc = np.arccos(self.bcirc/aRs)
        becc = self.bcirc * (1 - ecc**2)/(1 - esw)
        tdur = self.per / 2. / np.pi * np.arcsin(((1. + self.RpRs)**2 -
               becc**2)**0.5 / (np.sin(inc) * aRs))
        tdur *= np.sqrt(1. - ecc**2.)/(1. - esw)
        return tdur
           
class LIMBDARK(ctypes.Structure):
      '''
      The class containing the limb darkening parameters
      
      '''
      
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
        self.update(**kwargs)
        
      def update(self, **kwargs):
        '''
        
        '''
        
        self.ldmodel = kwargs.pop('ldmodel', QUADRATIC)
        self.u1 = kwargs.pop('u1', 0.40)                                              # ~ The sun seen by Kepler (http://arxiv.org/pdf/0912.2274.pdf)
        self.u2 = kwargs.pop('u2', 0.26)
        self.q1 = kwargs.pop('q1', np.nan)
        self.q2 = kwargs.pop('q2', np.nan)
        self.c1 = kwargs.pop('c1', np.nan)
        self.c2 = kwargs.pop('c2', np.nan)
        self.c3 = kwargs.pop('c3', np.nan)
        self.c4 = kwargs.pop('c4', np.nan)
        
        # If other coeffs are set, clear the defaults
        if (not np.isnan(self.q1)) or (not np.isnan(self.c1)):
          self.u1 = np.nan
          self.u2 = np.nan 
                  
class ARRAYS(ctypes.Structure):
      '''
      The class that stores the input and output arrays
      
      '''
      
      _fields_ = [("nstart", ctypes.c_int),
                  ("nend", ctypes.c_int),
                  ("ipts", ctypes.c_int),
                  ("_calloc", ctypes.c_int),
                  ("_balloc", ctypes.c_int),
                  ("_ialloc", ctypes.c_int),
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
        self.nstart = 0
        self.nend = 0
        self.ipts = 0
        self._calloc = 0
        self._balloc = 0
        self._ialloc = 0
      
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
      '''
      The class that contains the light curve settings
      
      '''
      _fields_ = [("exptime", ctypes.c_double),
                  ("keptol", ctypes.c_double),
                  ("fullorbit", ctypes.c_int),
                  ("maxpts", ctypes.c_int),
                  ("exppts", ctypes.c_int),
                  ("binmethod", ctypes.c_int),
                  ("intmethod", ctypes.c_int),
                  ("maxkepiter", ctypes.c_int),
                  ("computed", ctypes.c_int),
                  ("binned", ctypes.c_int),
                  ("kepsolver", ctypes.c_int)]
      
      def __init__(self, **kwargs):
        self.exptime = KEPLONGEXP
        self.fullorbit = 0
        self.maxpts = 10000
        self.exppts = 50
        self.binmethod = RIEMANN
        self.intmethod = SMARTINT
        self.keptol = 1.e-15
        self.maxkepiter = 100
        self.kepsolver = NEWTON
        self.update(**kwargs)
      
      def update(self, **kwargs):
        '''
        
        '''
        
        self.exptime = kwargs.pop('exptime', self.exptime)                         # Long cadence integration time
        self.fullorbit = 1 if kwargs.pop('fullorbit', self.fullorbit) else 0          # Compute full orbit or just the transits (default)
        self.maxpts = kwargs.pop('maxpts', self.maxpts)                               # Maximum number of points in arrays (for mem. allocation)
        self.exppts = kwargs.pop('exppts', self.exppts)                               # Average flux over this many points for binning
        self.binmethod = kwargs.pop('binmethod', self.binmethod)                      # How to integrate when binning?
        self.intmethod = kwargs.pop('intmethod', self.intmethod)                      # Integration method
        self.keptol = kwargs.pop('keptol', self.keptol)                               # Kepler solver tolerance
        self.maxkepiter = kwargs.pop('maxkepiter', self.maxkepiter)                   # Maximum number of iterations in Kepler solver
        self.kepsolver = kwargs.pop('kepsolver', self.kepsolver)                      # Newton solver or fast M&D solver?
        self.computed = 0
        self.binned = 0

# Load the C library
try:
  lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'transitlib.so'))
except:
  raise Exception("Can't find `transitlib.so`; please run `make` to compile `pysyzygy`.")

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
  if (err == _ERR_NONE):
    return
  elif (err == _ERR_NOT_IMPLEMENTED):
    raise Exception("Option not implemented.")
  elif (err == _ERR_MAX_PTS):
    raise Exception("Maximum points in lightcurve exceeded. " + 
                    "Try decreasing `exppts`, increasing `exptime`, or increasing "+
                    "`maxpts`.")  
  elif (err == _ERR_NO_TRANSIT):
    raise Exception("Object does not transit the star.")  
  elif (err == _ERR_BAD_ECC):
    raise Exception("Bad value for ``ecc``.")  
  elif (err == _ERR_RC):
    raise Exception("Error in elliptic integral function RC().")  
  elif (err == _ERR_RJ):
    raise Exception("Error in elliptic integral function RJ().") 
  elif (err == _ERR_RF):
    raise Exception("Error in elliptic integral function RF().") 
  elif (err == _ERR_RADIUS):
    raise Exception("Bad value for ``RpRs``.") 
  elif (err == _ERR_EXP_PTS):
    raise Exception("The number of exposure points must be even.") 
  elif (err == _ERR_NOT_COMPUTED):
    raise Exception("Lightcurve must be computed before it can be binned.") 
  elif (err == _ERR_STAR_CROSS):
    raise Exception("Star-crossing orbit.") 
  elif (err == _ERR_RHOS_ARS):
    raise Exception("Must specify one of ``rhos`` or ``aRs``.") 
  elif (err == _ERR_RHOS):
    raise Exception("Bad value for ``per``.") 
  elif (err == _ERR_ECC_W):
    raise Exception("Bad value for ``esw`` or ``ecw``.") 
  elif (err == _ERR_LD):
    raise Exception("Bad value for the limb darkening coefficients.") 
  elif (err == _ERR_T0):
    raise Exception("Bad value for ``t0``.")
  elif (err == _ERR_KEPLER):
    raise Exception("Error in Kepler solver.")
  else:
    raise Exception("Error in transit computation (%d)." % err)

class Transit():
  '''
  A user-friendly wrapper around the :py:class:`ctypes` routines.
  
  :param kwargs: The :py:mod:`pysyzygy` keyword arguments. These are:
  
    - **b** or **bcirc** - The (circular) impact parameter. Default `0.`
    - **MpMs** - The planet-star mass ratio. Default `0.`
    - **per** - The planet orbital period in days. Default `10.`
    - **RpRs** - The planet-star radius ratio. Default `0.1`
    - **rhos** or **aRs** - The stellar density in `g/cm^3` or the semi-major axis-stellar radius ratio. \
                            Default is `rhos = 1.4`, the density of the Sun
    - **ecc** and **w** or **esw** and **ecw** - The eccentricity and the longitude of pericenter in radians, \
                                                 or the two eccentricity vectors. Default is `ecc = 0.` and `w = 0.`
    - **t0** or **times** - The time of first transit, or the time of each of the transits (in case \
                            they are not periodic) in days. Default is `t0 = 0.`
  
    - **ldmodel** - The limb darkening model. Default `ps.QUADRATIC` (only option for now!)
    - **u1** and **u2** or **q1** and **q2** - The quadratic limb darkening parameters (u1, u2) or the \
                                               modified quadratic limb darkening parameters (q1, q2) \
                                               from `Kipping (2013) <http://dx.doi.org/10.1093/mnras/stt1435>`_. \
                                               Default is `u1 = 0.40` and `u2 = 0.26`
    
    - **exptime** - The exposure time in days for binning the model. Default `ps.KEPLONGEXP`
    - **fullorbit** - Compute the orbital parameters for the entire orbit? Only useful if \
                      you're interested in the full arrays of orbital parameters. Default `False`
    - **maxpts** - Maximum number of points in the model. Increase this if you're getting errors. Default `10,000`
    - **exppts** - The number of exposure points per cadence when binning the model. Default `50`
    - **binmethod** - The binning method. Default `ps.RIEMANN` (recommended)
    - **intmethod** - The integration method. Default `ps.SMARTINT` (recommended)
    - **keptol** - The tolerance of the Kepler solver. Default `1.e-15`
    - **maxkepiter** - Maximum number of iterations in the Kepler solver. Default `100`
    - **kepsolver** - The Kepler solver to use. Default `ps.NEWTON` (recommended)

  Once a :py:class:`Transit` model is instantiated, it may be called as follows:
  
  .. code-block::python
      
      import pysyzygy as ps
      import numpy as np
      
      # Instantiate the model
      trn = ps.Transit()
      
      # The observation times
      time = np.linspace(0, 10, 100)
      
      # Obtain the actual model, evaluated on the time array
      model = trn(time)
  
  '''
  
  def __init__(self, **kwargs):
    self.arrays = ARRAYS()
    self.limbdark = LIMBDARK()
    self.transit = TRANSIT()
    self.settings = SETTINGS()
    self.update(**kwargs)
  
  def update(self, **kwargs):
    '''
    Update the transit keyword arguments
    
    '''
    
    if kwargs.get('verify_kwargs', True):
      valid = [y[0] for x in [TRANSIT, LIMBDARK, SETTINGS] for y in x._fields_]       # List of valid kwargs
      valid += ['b', 'times']                                                         # These are special!
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
  
  def __call__(self, t, param = 'binned'):
    if param == 'binned':
      array = _ARR_BFLX
    elif param == 'unbinned':
      array = _ARR_FLUX
    elif param == 'M':
      array = _ARR_M
    elif param == 'E':
      array = _ARR_E
    elif param == 'f':
      array = _ARR_F
    elif param == 'r':
      array = _ARR_R
    elif param == 'x':
      array = _ARR_X
    elif param == 'y':
      array = _ARR_Y
    elif param == 'z':
      array = _ARR_Z
    elif param == 'b':
      array = _ARR_B
    else:
      RaiseError(_ERR_NOT_IMPLEMENTED)
    
    # Ensure the time is a float array
    if not (type(t) is np.ndarray):
      t = np.array(t, dtype = 'float64')
    elif t.dtype != 'float64':
      t = np.array(t, dtype = 'float64')
    
    err = _Interpolate(t, len(t), array, self.transit, self.limbdark, self.settings, 
                       self.arrays)
    if err != _ERR_NONE: RaiseError(err)
    res = self.arrays.iarr
    _dbl_free(self.arrays._iarr)
    self.arrays._ialloc = 0
    return res
  
  def Compute(self):
    '''
    Computes the light curve model
    
    '''
    
    err = _Compute(self.transit, self.limbdark, self.settings, self.arrays)
    if err != _ERR_NONE: RaiseError(err)    

  def Bin(self):
    '''
    Bins the light curve model to the provided time array
    
    '''
    
    err = _Bin(self.transit, self.limbdark, self.settings, self.arrays)
    if err != _ERR_NONE: RaiseError(err)
  
  def Free(self):
    '''
    Frees the memory used by all of the dynamically allocated C arrays.
    
    '''

    if self.arrays._calloc:
      _dbl_free(self.arrays._time)
      _dbl_free(self.arrays._flux)
      _dbl_free(self.arrays._bflx)
      _dbl_free(self.arrays._M)
      _dbl_free(self.arrays._E)
      _dbl_free(self.arrays._f)
      _dbl_free(self.arrays._r)
      _dbl_free(self.arrays._x)
      _dbl_free(self.arrays._y)
      _dbl_free(self.arrays._z)
      self.arrays._calloc = 0
    if self.arrays._balloc:  
      _dbl_free(self.arrays._b)
      self.arrays._balloc = 0
    if self.arrays._ialloc:
      _dbl_free(self.arrays._iarr)
      self.arrays._ialloc = 0
  
  def __del__(self):
    '''
    Free the C arrays when the last reference to the class goes out of scope!
    
    '''
    self.Free()