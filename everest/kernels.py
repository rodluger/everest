#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
kernels.py
----------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .utils import LatexExp, LatexExpSq
import numpy as np
import george
import logging
log = logging.getLogger(__name__)

class Kernel(object):
  '''
  A generic kernel class. This allows us to use a kernel
  for both fitting the autocorrelation function and instantiating
  a ``george.GP`` object.
  
  '''
  
  def __init__(self):
    self.kernels = []
    self.params = np.array([], dtype = float)
    self.parnames = np.array([], dtype = str)
  
  def __mul__(self, other):
    return Prod(self, other)
  
  def __rmul__(self, other):
    return Prod(other, self)
  
  def __add__(self, other):
    return Sum(self, other)
  
  def __radd__(self, other):
    return Sum(other, self)
  
  def __getitem__(self, i):
    return self.params[i]
  
  def __setitem__(self, i, p):
    params = self.params
    params[i] = p
    self.params = params
    for k in self.kernels:
      k.params = p[:len(k.params)]
      p = p[len(k.params):]
  
  def __call__(self, t):
    return self._function(t)
  
  def george_kernel(self):
    # Return the george equivalent
    return self._george()
  
class Sum(Kernel):
  '''
  A kernel sum.
  
  '''
  
  def __init__(self, *kernels):
    super(Sum, self).__init__()
    kernels = list(kernels)
    for i, k in enumerate(kernels):
      if type(k) is float or type(k) is int:
        kernels[i] = Const(k)
    self.kernels = kernels
    self._params = np.concatenate([k.params for k in kernels])
    self.parnames = np.concatenate([k.parnames for k in kernels])
  
  @property
  def params(self):
    return self._params
    
  @params.setter
  def params(self, params):
    self._params = params
    for k in self.kernels:
      k.params = params[:len(k.params)]
      params = params[len(k.params):]
      
  def __repr__(self):
    return r' $+$ '.join([k.__repr__() for k in self.kernels])
  
  def _function(self, t):
    res = np.zeros_like(t)
    for k in self.kernels:
      res += k._function(t)
    return res
  
  def _george(self):
    res = None
    for k in self.kernels:
      if res is None:
        res = k._george()
      else:
        res += k._george()
    return res 
  
class Prod(Kernel):
  '''
  A kernel product.
  
  '''
    
  def __init__(self, *kernels):
    super(Prod, self).__init__()
    kernels = list(kernels)
    for i, k in enumerate(kernels):
      if type(k) is float or type(k) is int:
        kernels[i] = Const(k)
    self.kernels = kernels
    self._params = np.concatenate([k.params for k in kernels])
    self.parnames = np.concatenate([k.parnames for k in kernels])
  
  @property
  def params(self):
    return self._params
    
  @params.setter
  def params(self, params):
    self._params = params
    for k in self.kernels:
      k.params = params[:len(k.params)]
      params = params[len(k.params):]
  
  def __repr__(self):
    return r' $\times$ '.join([k.__repr__() for k in self.kernels])
  
  def _function(self, t):
    res = np.ones_like(t)
    for k in self.kernels:
      res *= k._function(t)
    return res
  
  def _george(self):
    res = None
    for k in self.kernels:
      if res is None:
        res = k._george()
      else:
        res *= k._george()
    return res 

class Const(Kernel): 
  '''
  A constant kernel.
  
  '''
    
  def __init__(self, amp = 1.):
    super(Const, self).__init__()
    self.params = np.array([amp])
    self.parnames = np.array(['amp'])
    
  def __repr__(self):
    return LatexExpSq(self.params[0])  
  
  def _function(self, t):
    return self.params[0] ** 2
  
  def _george(self):
    return george.kernels.ConstantKernel(self.params[0] ** 2)

class White(Kernel): 
  '''
  A white noise kernel.
  
  '''
    
  def __init__(self, amp = 1.):
    super(White, self).__init__()
    self.params = np.array([amp])
    self.parnames = np.array(['amp'])
    
  def __repr__(self):
    return r'$\mathbf{White}$(%s)' % LatexExpSq(self.params[0]) 
  
  def _function(self, t):
    return (np.array(t) == 0) * self.params[0] ** 2
  
  def _george(self):
    return george.kernels.WhiteKernel(self.params[0] ** 2)
    
class Exp(Kernel): 
  '''
  An exponential kernel.
  
  '''
    
  def __init__(self, tau = 1.):
    super(Exp, self).__init__()
    self.params = np.array([tau])
    self.parnames = np.array(['tau'])
    
  def __repr__(self):
    return r'$\mathbf{Exp}$(%s)' % LatexExp(self.params[0])
  
  def _function(self, t):
    return np.exp(-np.abs(t) / self.params[0])
  
  def _george(self):
    # As far as I understand, ExpKernel expects the **square** of the timescale,
    # even though it's never squared in the actual expression. TODO: Confirm this.
    return george.kernels.ExpKernel(self.params[0] ** 2)

class Cos(Kernel):  
  '''
  A cosine kernel.
  
  '''
    
  def __init__(self, per = 1.):
    super(Cos, self).__init__()
    self.params = np.array([per])
    self.parnames = np.array(['per'])
    
  def __repr__(self):
    return r'$\mathbf{Cos}$(%s)' % LatexExp(self.params[0])  
  
  def _function(self, t):
    return np.cos(2 * np.pi / self.params[0] * t)
  
  def _george(self):
    return george.kernels.CosineKernel(self.params[0])

class Mat32(Kernel):
  '''
  A Matern-3/2 kernel.
  
  '''
    
  def __init__(self, tau = 1.):
    super(Mat32, self).__init__()
    self.params = np.array([tau])
    self.parnames = np.array(['tau'])
  
  def __repr__(self):
    return r'$\mathbf{Mat}$(%s)' % LatexExpSq(self.params[0])
  
  def _function(self, t):
    return (1 + np.sqrt(3) * t / self.params[0]) * np.exp(-np.sqrt(3) * t / self.params[0])
  
  def _george(self):
    return george.kernels.Matern32Kernel(self.params[0] ** 2)

# Our basis set of kernel models used for GP optimization.
KernelModels = [

  White() + Const() * Exp(),
  White() + Const() * Exp() * Cos(),
  White() + Const() * Exp() * Cos() * Cos(),
  White() + Const() * Exp() * Cos() * Cos() * Cos(),
  
  White() + Const() * Exp() + Const() * Cos(),
  White() + Const() * Exp() + Const() * Exp() * Cos(),
  White() + Const() * Exp() * Cos() + Const() * Exp() * Cos(),
  
  White() + Const() * Mat32(),
  White() + Const() * Mat32() * Cos(),
  White() + Const() * Mat32() * Cos() * Cos(),
  White() + Const() * Mat32() * Cos() * Cos() * Cos(),
  
  White() + Const() * Mat32() + Const() * Cos(),
  White() + Const() * Mat32() + Const() * Mat32() * Cos(),
  White() + Const() * Mat32() * Cos() + Const() * Mat32() * Cos(),

  White() + Const() * Exp() + Const() * Mat32() * Cos(),
  White() + Const() * Mat32() + Const() * Exp() * Cos(),
  White() + Const() * Exp() * Cos() + Const() * Mat32() * Cos()
  
]