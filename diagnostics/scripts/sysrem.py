import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as pl
from everest2.config import EVEREST_DAT
from everest2.missions.k2.aux import GetK2Campaign, GetNeighboringChannels
from wpca import PCA, WPCA
import george
from scipy.signal import savgol_filter, medfilt, periodogram
import george

# NOTE:
# Channel 30 is interesting. Stars #9 and #34 single-handedly screw up PCA.

# User
npc = 2
campaign = 6.0
model = 'SegmentedPLD'
ni = 7
nj = 6
wout = False

# Setup
X = None
W = None
fig1, ax1 = pl.subplots(ni, nj, sharex = True)
fig1.subplots_adjust(wspace = 0, hspace = 0)
ax1 = ax1.flatten()
[axis.axis('off') for axis in ax1]
all = GetK2Campaign(campaign)
stars = np.array([s[0] for s in all], dtype = int)
channels = np.array([s[2] for s in all], dtype = int)
N = 3000

# Loop
j = 0
for i in range(N):
  sys.stdout.write('\rProcessing target %d/%d...' % (j + 1, N))
  sys.stdout.flush()
  try:
    nf = os.path.join(EVEREST_DAT, 'k2', 'c%02d' % int(campaign),
                   ('%09d' % stars[i])[:4] + '00000', 
                   ('%09d' % stars[i])[4:], model + '.npz')
  
    data = np.load(nf)
  except:
    continue
    
  # Get arrays
  m = np.array(list(set(np.concatenate([data['outmask'], data['badmask'], data['nanmask']]))), dtype = int)
  t = data['time']
  y = data['fraw'] - data['model']  
    
  # debug
  if channels[i] not in [30, 31, 32, 33]:
    continue  
    
  # Throw out light curves with more than 20% variation above/below the median,
  # as they could mess up PCA
  z = np.delete(y, m)
  mz = np.nanmedian(z)
  if len(np.where((z < 0.8 * mz) | (z > 1.2 * mz))[0]):
    continue
    
  # Weighting
  w = np.ones_like(t) * np.sqrt(mz)
  y = np.interp(t, np.delete(t, m), np.delete(y, m))
  y /= mz
    
  if X is None:
    X = y.reshape(-1, 1)
    W = w.reshape(-1, 1)
  else:
    X = np.hstack([X, y.reshape(-1, 1)])
    W = np.hstack([W, w.reshape(-1, 1)])
  
  # Plot
  if j < ni * nj:
    ax1[j].plot(t, y, 'k.', alpha = 0.3, markersize = 2)    
  j += 1

# Remove global outliers (TODO: add them back in below)
bad = np.where(~W.any(axis=1))[0]
W = np.delete(W, bad, axis = 0)
X = np.delete(X, bad, axis = 0)
t = np.delete(t, bad)

# PCA
print("\nRunning PCA...")
pca = WPCA(n_components = npc)
xpca = pca.fit_transform(X, weights = W)
XPCA = np.hstack([xpca, np.ones(xpca.shape[0]).reshape(-1, 1)])

# Plot the principal components
fig2, ax2 = pl.subplots(npc, sharex = True)
for n in range(npc):
  ax2[n].plot(t, XPCA[:,n], 'k-')
 
# Show 
print("Explained variance: ", pca.explained_variance_ratio_)
pl.show()
quit()

# Fit the principal components to the first few targets
for i in range(20):
  nf = os.path.join(EVEREST_DAT, 'k2', 'c%02d' % int(campaign),
                   ('%09d' % stars[i])[:4] + '00000', 
                   ('%09d' % stars[i])[4:], model + '.npz')
  data = np.load(nf)
  m = np.array(list(set(np.concatenate([data['outmask'], data['badmask'], data['nanmask']]))), dtype = int)
  t = data['time']
  y = data['fraw'] - data['model'] 
  e = data['fraw_err']
  _, amp, tau = data['kernel_params']
  #gp = george.GP(amp ** 2 * george.kernels.Matern32Kernel(tau ** 2))
  gp = george.GP(george.kernels.WhiteKernel(0.))
  gp.compute(np.delete(t, m), np.delete(e, m))
  mX = np.delete(XPCA, m, axis = 0)
  A = np.dot(mX.T, gp.solver.apply_inverse(mX))
  B = np.dot(mX.T, gp.solver.apply_inverse(np.delete(y, m)))
  C = np.linalg.solve(A, B)
  pcamodel = np.dot(XPCA, C)
  pcamodel -= np.nanmedian(pcamodel)
  fig3, ax3 = pl.subplots(2)
  ax3[0].plot(np.delete(t, m), np.delete(y, m), 'k.', alpha = 0.3)
  ax3[0].plot(np.delete(t, m), np.delete(pcamodel, m) + np.nanmedian(y), 'r-', alpha = 0.3)
  ax3[1].plot(np.delete(t, m), np.delete(y, m) - np.delete(pcamodel, m), 'k.', alpha = 0.3)
  pl.show()
  pl.close()

