import numpy as np
import everest
from everest.gp import GetCovariance
from everest.transit import TransitShape, TransitModel, Get_rhos
import matplotlib.pyplot as pl
from scipy.linalg import cholesky, cho_solve
from tqdm import tqdm
import os

# The injected depth
depth = 0.001

# Load the light curve
star = everest.Everest(201367065, quiet = True)
overfit = star.overfit(plot = False)

# Compute
if not os.path.exists('overfit.npz'):

    # Pre-compute some stuff
    med = np.nanmedian(star.flux)
    mtime = star.apply_mask(star.time)
    merr = star.apply_mask(star.fraw_err)
    npts = len(mtime)
    inj = TransitShape(depth = 1)
    K = GetCovariance(star.kernel, star.kernel_params, mtime, merr)
    CK = cholesky(K)
    fpix = np.array(star.fpix)

    # Unmasked Overfitting Metric (Brute force)
    UOMbf = np.zeros(npts)

    # Masked Overfitting Metric (Brute force)
    MOMbf = np.zeros(npts)

    # Joint Overfitting Metric (Brute force)
    JOMbf = np.zeros(npts)

    # Loop over the light curve
    for n, t0 in tqdm(enumerate(star.apply_mask(star.time)), total = npts):
  
      # Inject the synthetic transit
      tau = inj(star.time, t0)
      star.fpix = fpix * (1. + depth * tau.reshape(-1,1))
      star.fraw = np.sum(star.fpix, axis = 1)
      star.get_norm()
      
      # Compute the unmasked overfitting metric
      star.compute()
      mflux = star.apply_mask(star.flux)
      mtau = star.apply_mask(tau)
      d = np.dot(mtau, cho_solve((CK, False), mflux)) \
        / np.dot(mtau, cho_solve((CK, False), mtau)) \
        / med
      UOMbf[n] = 1 - d / depth

      # Compute the masked overfitting metric
      star.transitmask = np.where(tau < 0)[0]
      star.compute()
      star.transitmask = []
      mflux = star.apply_mask(star.flux)
      d = np.dot(mtau, cho_solve((CK, False), mflux)) \
        / np.dot(mtau, cho_solve((CK, False), mtau)) \
        / med
      MOMbf[n] = 1 - d / depth
     
      # Compute the joint overfitting metric
      star.transit_model = TransitModel('injected', sig_RpRs = depth / 100.,
                                        t0 = t0, per = 3.56789,
                                        rhos = Get_rhos(0.1, per = 3.56789),
                                        single = True)
      star.compute_joint()
      d = star.transit_depth
      star.transit_model = None 
      JOMbf[n] = 1 - d / depth

    # Save!
    np.savez('overfit.npz', UOMbf = UOMbf, MOMbf = MOMbf, JOMbf = JOMbf)

else:
    
    data = np.load('overfit.npz')
    UOMbf = data['UOMbf']
    MOMbf = data['MOMbf']
    JOMbf = data['JOMbf']

# Plot
fig, ax = pl.subplots(3, figsize = (10, 6))

# Unmasked
ax[0].plot(star.apply_mask(star.time), star.apply_mask(1 - (overfit._O2 + overfit._O3 / depth) / overfit._O1), 'k-', lw = 1)
ax[0].plot(star.apply_mask(star.time), UOMbf, 'r-', lw = 1, alpha = 0.5)
ax[0].set_xlim(star.apply_mask(star.time)[30],star.apply_mask(star.time)[-1])
ax[0].set_ylim(-1.5, 1.5)

# Masked
ax[1].plot(star.apply_mask(star.time), star.apply_mask(overfit._O5), 'k-', lw = 1)
ax[1].plot(star.apply_mask(star.time), MOMbf, 'r-', lw = 1, alpha = 0.5)
ax[1].set_xlim(star.apply_mask(star.time)[30],star.apply_mask(star.time)[-1])
#ax[1].set_ylim(-1.5, 1.5)

# Joint
ax[2].plot(star.apply_mask(star.time), JOMbf, 'r-', lw = 1, alpha = 0.5)
ax[2].set_xlim(star.apply_mask(star.time)[30],star.apply_mask(star.time)[-1])
#ax[2].set_ylim(-1.5, 1.5)

pl.show()