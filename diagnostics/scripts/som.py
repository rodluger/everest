import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
from everest2.config import EVEREST_DAT
from everest2.math import Chunks
from everest2.missions.k2.aux import GetK2Campaign, GetNeighboringChannels
from wpca import PCA, WPCA
import george
from scipy.signal import savgol_filter, medfilt, periodogram
import george

def som_like(residuals, **kwargs):
  '''
  SOM likelihood function. Nothing fancy.
  
  '''
  
  return -0.5 * np.sum(residuals ** 2)

# User
npc = 3
ki = 3
kj = 1
T = 100
alpha_0 = 0.1
sig_0 = max(ki, kj)
periodic = False
channels = range(99)
campaign = 6.0
model = 'SegmentedPLD'
red = False
calc = False
smooth = True

# Compute the common design matrix
if calc:

  # Get the light curves
  all = GetK2Campaign(campaign)
  stars = np.array([s[0] for s in all if s[2] in channels and 
          os.path.exists(os.path.join(EVEREST_DAT, 'k2', 'c%02d' % int(campaign),
          ('%09d' % s[0])[:4] + '00000', 
          ('%09d' % s[0])[4:], model + '.npz'))], dtype = int)
  N = len(stars)
  nflx = []
  flux = []
  ferr = []
  time = None
  for n in range(N):
    nf = os.path.join(EVEREST_DAT, 'k2', 'c%02d' % int(campaign),
                   ('%09d' % stars[n])[:4] + '00000', 
                   ('%09d' % stars[n])[4:], model + '.npz')
    # Get the data
    data = np.load(nf)
    t = data['time']
    if time is None:
      time = t
    fraw_err = data['fraw_err']
    
    # Get de-trended light curve
    y = data['fraw'] - data['model'] 
    
    # Interpolate over outliers 
    m = np.array(list(set(np.concatenate([data['outmask'], data['badmask'], data['nanmask']]))), dtype = int)
    y = np.interp(t, np.delete(t, m), np.delete(y, m))
    
    # Throw out light curves with more than 20% variation 
    # above/below the median, as they could mess up PCA
    my = np.nanmedian(y)
    if len(np.where((y < 0.8 * my) | (y > 1.2 * my))[0]):
      continue
    
    # Append to our lists
    flux.append(y)
    ferr.append(fraw_err)
    nflx.append(((y - np.min(y)) / (np.max(y) - np.min(y))))

  # Initialize the SOM
  nflx = np.array(nflx)
  flux = np.array(flux)
  ferr = np.array(ferr)
  N = len(nflx)
  t = 0
  alpha = alpha_0
  sig = sig_0
  SOM = np.random.random((ki, kj, len(time))) 

  # Do T iterations of the SOM
  for t in range(T):
  
    print("Iteration %d/%d..." % (t + 1, T))
  
    # Loop over all light curves
    for n in range(N):
  
      # Find best matching pixel for this light curve
      bi, bj = np.unravel_index(np.argmax([som_like(nflx[n] - SOM[i][j]) for i in range(ki) for j in range(kj)]), (ki, kj))
    
      # Flip light curve upside-down and see if we get a better match
      bi_, bj_ = np.unravel_index(np.argmax([som_like(1 - nflx[n] - SOM[i][j]) for i in range(ki) for j in range(kj)]), (ki, kj))
      if som_like(1 - nflx[n] - SOM[bi_][bj_]) > som_like(nflx[n] - SOM[bi][bj]):
        nflx[n] = 1 - nflx[n]
        bi, bj = bi_, bj_
    
      # Update the Kohonen layer
      for i in range(ki):
        for j in range(kj):
        
          # Compute distance
          if periodic:
            dx2 = min((bi - i) ** 2, (ki - (bi - i)) ** 2)
            dy2 = min((bj - j) ** 2, (kj - (bj - j)) ** 2)
          else:
            dx2 = (bi - i) ** 2
            dy2 = (bj - j) ** 2
          d2 = dx2 + dy2
        
          # Update the pixel
          SOM[i][j] += alpha * np.exp(-d2 / (2 * sig ** 2)) * (nflx[n] - SOM[i][j])
  
    # Update the params
    sig = sig_0 * np.exp(-t / T * np.log(max(ki, kj)))
    alpha = alpha_0 * (1 - t / T)

  # Find the best matching pixel for each light curve
  pixel_lcs = [[[] for j in range(kj)] for i in range(ki)]
  evratio = [[0 for j in range(kj)] for i in range(ki)]
  for n in range(N):
    bi, bj = np.unravel_index(np.argmin([np.sum((nflx[n] - SOM[i][j]) ** 2) for i in range(ki) for j in range(kj)]), (ki, kj))
    pixel_lcs[bi][bj].append(n)

  # Get the top principal components in each pixel
  pc = [[[] for j in range(kj)] for i in range(ki)]
  for i in range(ki):
    for j in range(kj):
      n = pixel_lcs[i][j]
      pca = WPCA(n_components = npc)
      pc[i][j] = pca.fit_transform((flux[n] / np.nanmedian(flux[n], axis = 1).reshape(-1, 1)).T, weights = np.array([np.nanmedian(np.sqrt(f)) * np.ones_like(f) for f in flux[n]]).T).T
      evratio[i][j] = pca.explained_variance_ratio_[0]
      
      # Apply a smoothing filter?
      if smooth:
        for n in range(npc):
          pc[i][j][n] = savgol_filter(pc[i][j][n], 99, 2)
      
  # Plot
  fig, ax = pl.subplots(1, figsize = (ki * 5, kj * 5))
  for i in range(ki):
    for j in range(kj):
      for n in pixel_lcs[i][j]:
        ax.plot(np.linspace(i + 0.5 - 0.4, i + 0.5 + 0.4, len(time)), j + 0.5 - 0.4 + 0.8 * nflx[n], 'k-', alpha = min(1., 100. / N))
      for p, c in zip(pc[i][j], ['r', 'b', 'g', 'y', 'purple']):
        ax.plot(np.linspace(i + 0.5 - 0.4, i + 0.5 + 0.4, len(time)), j + 0.5 - 0.4 + 0.8 * (p - np.min(p)) / (np.max(p) - np.min(p)), ls = '-', color = c)
      ax.annotate('%d (%d)' % (len(pixel_lcs[i][j]), len(pixel_lcs[i][j]) * evratio[i][j]), xy = (i + 0.5, j + 0.9), xycoords = 'data', va = 'center', ha = 'center', zorder = 99, fontsize = 18)
  ax.get_xaxis().set_major_locator(MaxNLocator(integer=True))
  ax.get_yaxis().set_major_locator(MaxNLocator(integer=True))
  ax.grid(True, which = 'major', axis = 'both', linestyle = '-')
  ax.set_xlim(0,ki)
  ax.set_ylim(0,kj)
  pl.show()

  # Get the principal components in the best SOM pixel
  # and create a simple design matrix for fitting
  i, j = np.unravel_index(np.argmax([len(pixel_lcs[i][j]) * evratio[i][j] for i in range(ki) for j in range(kj)]), (ki, kj))
  X = np.hstack([pc[i][j].T, np.ones_like(time).reshape(-1, 1)])

  # Save
  np.savez('som.npz', X = X, flux = flux, ferr = ferr, time = time)

else:

  # Load
  data = np.load('som.npz')
  X, flux, ferr, time = data['X'], data['flux'], data['ferr'], data['time']

# Now apply the fit to each light curve
for n in range(len(flux)):

  # GP regression
  if red:
    gp = george.GP((0.005 * np.nanmedian(flux[n])) ** 2 * george.kernels.Matern32Kernel(1000. ** 2))
  else:
    gp = george.GP(george.kernels.WhiteKernel(0.))
  gp.compute(time, ferr[n])
  A = np.dot(X.T, gp.solver.apply_inverse(X))
  B = np.dot(X.T, gp.solver.apply_inverse(flux[n]))
  C = np.linalg.solve(A, B)
  pcamodel = np.dot(X, C)
  if red:
    # Add the red component to the model prediction
    y, _ = gp.predict(flux[n] - pcamodel, time)
    pcamodel += y
  
  fig, ax = pl.subplots(2, figsize = (12, 6))
  ax[0].plot(time, flux[n], 'k.', markersize = 3)
  ax[0].plot(time, pcamodel, 'r-')
  ax[1].plot(time, flux[n] - pcamodel + np.nanmedian(flux[n]), 'k.', markersize = 3)
  pl.show()