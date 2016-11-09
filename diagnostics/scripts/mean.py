import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as pl
from everest2.config import EVEREST_DAT
from everest2.math import Chunks
from everest2.missions.k2.aux import GetK2Campaign, GetNeighboringChannels
from wpca import PCA, WPCA
import george
from scipy.signal import savgol_filter, medfilt, periodogram
import george

# Get the light curves
npc = 5
channels = range(10)
campaign = 6.0
model = 'SegmentedPLD'
all = GetK2Campaign(campaign)
stars = np.array([s[0] for s in all if s[2] in channels and 
        os.path.exists(os.path.join(EVEREST_DAT, 'k2', 'c%02d' % int(campaign),
        ('%09d' % s[0])[:4] + '00000', 
        ('%09d' % s[0])[4:], model + '.npz'))], dtype = int)
N = len(stars)
Y = []
W = []
M = []
for n in range(N):
  nf = os.path.join(EVEREST_DAT, 'k2', 'c%02d' % int(campaign),
                 ('%09d' % stars[n])[:4] + '00000', 
                 ('%09d' % stars[n])[4:], model + '.npz')
  data = np.load(nf)
  m = np.array(list(set(np.concatenate([data['outmask'], data['badmask'], data['nanmask']]))), dtype = int)
  t = data['time']
  y = data['fraw'] - data['model']  
  y = np.interp(t, np.delete(t, m), np.delete(y, m))
  
  # Skip really bad light curves
  my = np.nanmedian(y)
  if len(np.where((y < 0.8 * my) | (y > 1.2 * my))[0]):
    continue
  
  # Normalize
  y /= np.nanmedian(y)
  
  # Tag outliers
  w = np.ones_like(t); w[m] = 0

  # Cumulative lists
  M.append(np.nanmedian(data['fraw']))
  Y.append(y)
  W.append(w)
  
Y = np.array(Y)
W = np.array(W)
bad = np.where(np.sum(W, axis = 0) == 0)
mu = np.average(Y, axis = 0, weights = np.sqrt(M)).reshape(-1,1)
X = np.hstack([mu, np.ones(mu.shape[0]).reshape(-1, 1)])

for y in Y:
  
  gp = george.GP(george.kernels.WhiteKernel(0.))
  gp.compute(t)
  A = np.dot(X.T, gp.solver.apply_inverse(X))
  B = np.dot(X.T, gp.solver.apply_inverse(y))
  C = np.linalg.solve(A, B)
  pcamodel = np.dot(X, C)
  
  pl.plot(t, y, 'k.', ms = 3)
  pl.plot(t, pcamodel, 'r-')
  pl.show()
