from __future__ import division, print_function, absolute_import, unicode_literals
from .utils import Chunks, RMS
from .data import GetData, GetGoodNeighbors
from .gp import GetGP
from .model import Model
import numpy as np
import matplotlib.pyplot as pl
from scipy.signal import savgol_filter
from itertools import combinations_with_replacement as multichoose

class cPLD(Model):
  '''
  EXPERIMENTAL
  
  '''
  
  def __init__(self, EPIC, neighbors = 3, pld_order = 3):
    '''
    
    '''
    
    # Initialize
    self.EPIC = EPIC
    
    # Load target data
    data = GetData(self.EPIC, interp_nans = True, interp_outliers = False)
    self.time = data.time
    self.model = np.zeros_like(self.time)
    self.fpix = data.fpix
    self.fraw = np.sum(self.fpix, axis = 1)
    self.fpix_err = data.fpix_err
    self.fraw_err = np.sqrt(np.sum(self.fpix_err ** 2, axis = 1))
    self.mask = data.mask
    self.M = lambda x: np.delete(x, self.mask, axis = 0)
    self.aperture = data.aperture
        
    # Get target GP (TODO)
    self.gp = GetGP(self.EPIC)
    self.K = np.diag(self.M(self.fraw_err) ** 2)
    self.K += self.gp.get_matrix(self.M(self.time))
        
    # Get neighbors?
    self.X = np.empty((self.time.shape[0], 0), dtype = 'float64')
    for i, star in enumerate(GetGoodNeighbors(self.EPIC)):
      X1 = star.fpix / np.sum(star.fpix, axis = 1).reshape(-1, 1)
      for n in range(1, pld_order + 1):  
        x = np.product(list(multichoose(X1.T, n)), axis = 1).T
        self.X = np.hstack([self.X, x])
      if i >= neighbors - 1:
        break
    
    # debug. Beat 69.83 ppm --> 77.30 (63.04) --> 50ish
    niter = 5
    
    # First order PLD
    lam = 1e9
    self.X0 = np.array(self.X)
    for j in range(niter):
      print("Iteration %d/%d..." % (j + 1, niter))
      # Compute normalized pixels, ``X1``
      A, B, _ = self._precompute()
      X1 = np.empty((self.time.shape[0], 0), dtype = 'float64')
      for i in range(self.fpix.shape[1]):
        f = self.M(self.fpix[:,i]) - np.median(self.fpix[:,i])
        W = np.linalg.solve(self.K + lam * A, f)
        model = np.dot(lam * B, W)
        x = (self.fpix[:,i] / (self.fpix[:,i] - model))        
        mask = np.where(np.isinf(x))[0]
        if len(mask):
          t_ = np.delete(self.time, mask)
          x_ = np.delete(x, mask, axis = 0)
          x[mask] = np.interp(self.time[mask], t_, x_)
        X1 = np.hstack([X1, x.reshape(-1, 1)])
      # Add in the neighbors
      self.X = np.hstack([self.X0, X1])
      
    # Second order PLD
    lam = 1e0
    X2 = np.product(list(multichoose(X1.T, 2)), axis = 1).T
    self.X = np.hstack([self.X0, X1, X2])
    for j in range(niter):
      print("Iteration %d/%d..." % (j + 1, niter))
      # Compute normalized pixels, ``X1``
      A, B, _ = self._precompute()
      X1 = np.empty((self.time.shape[0], 0), dtype = 'float64')
      for i in range(self.fpix.shape[1]):
        f = self.M(self.fpix[:,i]) - np.median(self.fpix[:,i])
        W = np.linalg.solve(self.K + lam * A, f)
        model = np.dot(lam * B, W)
        x = (self.fpix[:,i] / (self.fpix[:,i] - model))        
        mask = np.where(np.isinf(x))[0]
        if len(mask):
          t_ = np.delete(self.time, mask)
          x_ = np.delete(x, mask, axis = 0)
          x[mask] = np.interp(self.time[mask], t_, x_)
        X1 = np.hstack([X1, x.reshape(-1, 1)])
      X2 = np.product(list(multichoose(X1.T, 2)), axis = 1).T
      self.X = np.hstack([self.X0, X1, X2])
    
    self.cross_validate(niter = 3, lambda_arr = np.logspace(-10, 20, 30))
    import pdb; pdb.set_trace()

class Outliers(Model):
  '''
  
  '''
  
  def __init__(self, EPIC, pld_order = 3, gp = None, sigma = 5, maxiter = 25, **kwargs):
    '''

    '''
    
    # Initialize
    self.EPIC = EPIC
    self.lambda_opt = 1.e9
    
    # Load target data
    data = GetData(self.EPIC)
    self.time = data.time
    self.model = np.zeros_like(self.time)
    self.fpix = data.fpix
    self.fraw = np.sum(self.fpix, axis = 1)
    self.fpix_err = data.fpix_err
    self.fraw_err = np.sqrt(np.sum(self.fpix_err ** 2, axis = 1))
    self.nanmask = data.nanmask
    self.outmask = data.outmask
    self.aperture = data.aperture
    
    # Get target GP
    if gp is None:
      amp = np.nanmedian([np.nanstd(y) for y in Chunks(self.fraw, int(2. / np.nanmedian(self.time[1:] - self.time [:-1])))])
      self.gp = george.GP(amp ** 2 * george.kernels.Matern32Kernel(2. ** 2))
    self.K = np.diag(self.fraw_err ** 2)
    self.K += self.gp.get_matrix(self.time)
            
    # Get target pixels
    log.info("Normalizing target basis vectors...")
    self.X = np.empty((self.time.shape[0], 0), dtype = 'float64')
    X1 = self.fpix / self.flux.reshape(-1, 1)
    for n in range(1, pld_order + 1):  
      x = np.product(list(multichoose(X1.T, n)), axis = 1).T
      self.X = np.hstack([self.X, x])

    # Cross-validate to get lambda
    self.cross_validate(**kwargs)

    # Now do iterative sigma clipping
    log.info("Clipping outliers...")
    log.info('Iteration #00: %d outliers.' % len(self.outmask))
    count = 0
    outmask = np.array([-1])
    while not np.array_equal(self.outmask, outmask):
    
      # Reset; the loop ends when we get the same outliers twice in a row
      outmask = np.array(self.outmask)
      count += 1
      if count > maxiter:
        # We will continue, even though there may be issues
        log.error('Maximum number of iterations in ``Outliers()`` exceeded.')
        break
      
      # Compute the model to get the flux
      self.compute()
      
      # Clip!
      f = self.savgol(self.flux)
      med = np.nanmedian(f)
      MAD = 1.4826 * np.nanmedian(np.abs(f - med))
      self.outmask = np.where((f > med + sigma * MAD) | (f < med - sigma * MAD))[0]

      # Log
      log.info('Iteration #%02d: %d outliers.' % (count, len(self.outmask)))
      
      # DEBUG
      fig = pl.figure()
      pl.plot(self.M(self.time), self.M(self.flux), 'b.', alpha = 0.3)
      pl.plot(self.time[self.outmask], self.flux[self.outmask], 'r.', alpha = 0.3)
    pl.show()
    
    @property
    def outliers(self):
      '''
      
      '''
      
      return self.outmask


  class nPLD(Model):
  '''
  Neighbors PLD. **TODO**
  
  '''
  
  def __init__(self, EPIC, neighbors = 0, pld_order = 3, riter = 3, **kwargs):
    '''
    
    '''
    
    # Initialize
    self.EPIC = EPIC
    self.lambda_opt = 1.e9
    self.pld_order = pld_order
    self.riter = riter
    self.neighbors = neighbors
    self.fig = [None for i in range(self.riter + 1)]
    self.ax = [None for i in range(self.riter + 1)]
    
    # Load target data
    log.info("Loading target data...")
    data = GetData(self.EPIC)
    self.time = data.time
    self.model = np.zeros_like(self.time)
    self.fpix = data.fpix
    self.fraw = np.sum(self.fpix, axis = 1)
    self.fpix_err = data.fpix_err
    self.fraw_err = np.sqrt(np.sum(self.fpix_err ** 2, axis = 1))
    self.nanmask = data.nanmask
    self.badmask = data.badmask
    self.outmask = np.array([], dtype = int)
    self.aperture = data.aperture
    self.cdpp6 = np.nan
    self.cdpp6v = np.nan
    self.kernel_params = None
    log.info("Raw: CDPP = %.2f ppm" % (self.get_cdpp6()))
        
    # Get target GP
    self.compute_linear()
    log.info("Linear: CDPP = %.2f ppm" % (self.get_cdpp6()))
    self.K, self.kernel_params = GetCovariance(self.time, self.flux, self.fraw_err, 
                                               mask = self.mask, 
                                               guess = self.kernel_params, **kwargs)
    
    # Get neighbor pixels
    self.XNeigh = np.empty((self.time.shape[0], 0), dtype = 'float64')
    if self.neighbors > 0:
      for i, star in enumerate(GetGoodNeighbors(self.EPIC)):
        
        # TODO!!!
        raise Exception('Code something up to obtain the neighbor outliers!')
        star.fpix = Interpolate(star.time, mask, star.fpix)

        # Get PLD vectors
        log.info("Computing neighbor basis vectors (%d/%d)..." % (i + 1, neighbors))
        X1 = star.fpix / np.sum(star.fpix, axis = 1).reshape(-1, 1)
        for n in range(1, self.pld_order + 1):  
          x = np.product(list(multichoose(X1.T, n)), axis = 1).T
          self.XNeigh = np.hstack([self.XNeigh, x])
        if i == self.neighbors - 1:
          break
    
    # The main loop
    for r in range(self.riter):
    
      if r == 0:
        # Get outliers using a linear PLD model (so we don't have to cross-validate)
        self.get_outliers(linear = True, **kwargs)
        # Reset the model
        self.model = np.zeros_like(self.time)
      else:
        # Get outliers
        self.get_outliers(linear = False, **kwargs)
        # Update the GP
        self.K, self.kernel_params = GetCovariance(self.time, self.flux, self.fraw_err, 
                                                   mask = self.mask, 
                                                   guess = self.kernel_params, **kwargs)
    
      # Get target pixels
      log.info("Normalizing target basis vectors...")
      self.XSelf = np.empty((self.time.shape[0], 0), dtype = 'float64')
      X1 = self.fpix / self.flux.reshape(-1, 1)
      for n in range(1, self.pld_order + 1):  
        x = np.product(list(multichoose(X1.T, n)), axis = 1).T
        self.XSelf = np.hstack([self.XSelf, x])
      
      # Combine them and cross-validate to get optimal lambda
      self.X = np.hstack([self.XSelf, self.XNeigh])
      self.fig[r], self.ax[r] = self.cross_validate(**kwargs)
      log.info("Step %d/%d: CDPP = %.2f ppm (%.2f ppm)" % (r + 1, self.riter, self.cdpp6, self.cdpp6v))
      
    # Plot the final light curve
    self.fig[self.riter], self.ax[self.riter] = self.plot()