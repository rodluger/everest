# ----- These need work! -------

class s3nPLD(nPLD):
  '''
  Separately-regularized nPLD.
  
  '''
        
  def __init__(self, *args, **kwargs):
    '''
    
    '''
    
    # Initialize
    kwargs.update({'pld_order': kwargs.get('pld_order', 3) * 2})
    super(s3nPLD, self).__init__(*args, **kwargs)
        
  def get_X(self):
    '''
  
    '''
      
    if not self.is_parent:
      log.info("Computing the design matrix...")
    if self.recursive:
      X1 = self.fpix / self.flux.reshape(-1, 1)
    else:
      X1 = self.fpix / self.fraw.reshape(-1, 1)
    
    for i in range(self.pld_order // 2): 
      n = 2 * i
      if (self._X[n] is None) and ((n == self.lam_idx) or (self.lam[0][n] is not None)):
        self._X[n] = np.product(list(multichoose(X1.T, i + 1)), axis = 1).T
        self.optimize_gp = True
      if (self._X[n + 1] is None) and ((n + 1 == self.lam_idx) or (self.lam[0][n + 1] is not None)):  
        self._X[n + 1] = np.array(self._XNeighbors[i])
        self.optimize_gp = False
  
  def get_sc_X(self):
    '''
    
    '''
    
    # TODO!
    return None
  
  def plot_weights(self):
    '''
    
    '''
    
    # TODO!
    pass

class MotionVectorPLD(Model):
  '''
  
  '''
        
  def __init__(self, *args, **kwargs):
    '''
    
    '''
    
    # Initialize
    super(MotionVectorPLD, self).__init__(*args, **kwargs)
    
    # Check for saved model
    if self.load_model():
      return
    
    # Read model-specific kwargs
    self.motion_order = kwargs.get('motion_order', 10)
    assert (self.motion_order > 0), "Invalid value for the de-trending order."
    pld_order = kwargs.get('pld_order', 3)
    assert (pld_order > 0), "Invalid value for the de-trending order."
    self.nblock = pld_order + 1
    
    # Run
    self.run()
      
  def get_X(self, order = 1):
    '''
  
    '''
    
    log.info("Computing the design matrix...")
    self._X = []
    pld_order = order - 1
    assert (self.Xpos is not None) and (self.Ypos is not None), \
           'Position information not available for this target.'
    
    # Motion vectors
    XM = np.hstack([self.Xpos.reshape(-1, 1), self.Ypos.reshape(-1, 1)])
    for n in range(1, self.motion_order + 1):  
      x = np.product(list(multichoose(XM.T, n)), axis = 1).T
      self._X.append(x)
    
    # Pixels
    XP = self.fpix / self.flux.reshape(-1, 1)
    for n in range(1, pld_order + 1):  
      x = np.product(list(multichoose(XP.T, n)), axis = 1).T
      self._X.append(x)

class SingleLambdaPLD(Model):
  '''
  
  '''
        
  def __init__(self, *args, **kwargs):
    '''
    
    '''
    
    # Initialize
    super(SingleLambdaPLD, self).__init__(*args, **kwargs)
    
    # Check for saved model
    if self.load_model():
      return
    
    # Read model-specific kwargs
    self.pld_order = kwargs.get('pld_order', 3)
    assert (self.pld_order > 0), "Invalid value for the de-trending order."
    self.nblock = kwargs.get('liter', 2)
    self._get_X = True
    
    # Run
    self.run()
      
  def get_X(self, order = None):
    '''
  
    '''
    
    if self._get_X:
      log.info("Computing the design matrix...")
      X1 = self.fpix / self.flux.reshape(-1, 1)
      self._X = np.empty(shape = (self.fpix.shape[0], 0), dtype = 'float64')
      for n in range(1, self.pld_order + 1):  
        x = np.product(list(multichoose(X1.T, n)), axis = 1).T
        self._X = np.hstack([self._X, x])
      self._get_X = False

  def precompute(self):
    '''
    Pre-compute some expensive matrices.
    
    '''
    
    m = self.M(np.arange(len(self.time)))
    self._A = np.dot(self.X[m], self.X[m].T)
    self._B = np.dot(self.X, self.X[m].T)
    self._mK = self.M(self.M(self.K, axis = 0), axis = 1)
    self._f = self.M(self.fraw) - np.nanmedian(self.fraw)
    
  def compute(self, precompute = True):
    '''
    Compute the model for the current value of lambda.
    
    '''
    
    if precompute:
      self.precompute()
    W = np.linalg.solve(self._mK + self.lam[-1] * self._A, self._f)
    self.model = np.dot(self.lam[-1] * self._B, W) 
    self.model -= np.nanmedian(self.model)
            
  def _precompute(self, mask):
    '''
    Pre-compute the matrices A and B (cross-validation step only).
    
    '''
    
    # Mask transits and outliers
    flux = self.M(self.fraw)
    m1 = self.M(np.arange(len(self.time)))
    K = self.M(self.M(self.K, axis = 0), axis = 1)
    med = np.nanmedian(flux)
    
    # Now mask the validation set
    M = lambda x, axis = 0: np.delete(x, mask, axis = axis)
    m2 = M(m1)
    mK = M(M(K, axis = 0), axis = 1)
    f = M(flux) - med
    
    # Pre-compute the matrices
    A = np.dot(self.X[m2], self.X[m2].T)
    B = np.dot(self.X[m1], self.X[m2].T)
    
    return A, B, mK, f
    
  def _compute(self, A, B, mK, f):
    '''
    Compute the model (cross-validation step only).
    
    '''

    W = np.linalg.solve(mK + self.lam[-1] * A, f)
    model = np.dot(self.lam[-1] * B, W)
    model -= np.nanmedian(model)
    
    return model
    
class RefinedPLD(Model):
  '''
  Re-runs PLD to refine the choice of `lambda` for each PLD order.
  Doesn't typically improve things by that much, meaning the optimization is
  in fact separable in the different PLD orders.
  
  '''
  
  def __init__(self, *args, **kwargs):
    '''
    
    '''
    
    # Initialize
    super(RefinedPLD, self).__init__(*args, **kwargs)
    
    # Check for saved model
    if self.load_model():
      return
    
    # Get parent
    self.parent = kwargs.get('parent', 'StandardPLD')
    
    # Load the parent model: note that most kwargs get overwritten
    self.clobber = False
    if self.load_model(self.parent):
      self.lam_idx = -1
      self.run()
      return
    else:
      log.error('Unable to locate `%s` model for target.' % self.parent)
  
  def load_tpf(self):
    '''
    
    '''
    
    # We don't want to re-load the raw data!
    pass
      
  def get_X(self):
    '''
  
    '''
      
    if not self.is_child:
      log.info("Computing the design matrix...")
    X1 = self.fpix / self.fraw.reshape(-1, 1)
    for n in range(self.maxblocks): 
      if (self._X[n] is None) and ((n == self.lam_idx) or (self.lam[0][n] is not None)):
        self._X[n] = np.product(list(multichoose(X1.T, n + 1)), axis = 1).T