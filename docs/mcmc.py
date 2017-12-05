"""MCMC example for transit fitting."""

import matplotlib.pyplot as pl
from everest import Everest, TransitModel
import numpy as np
import emcee
from tqdm import tqdm
from corner import corner


def lnprior(x):
    """Return the log prior given parameter vector `x`."""
    per, t0, b = x
    if b < -1 or b > 1:
        return -np.inf
    elif per < 7 or per > 10:
        return -np.inf
    elif t0 < 1978 or t0 > 1979:
        return -np.inf
    else:
        return 0.


def lnlike(x, star):
    """Return the log likelihood given parameter vector `x`."""
    ll = lnprior(x)
    if np.isinf(ll):
        return ll, (np.nan, np.nan)
    per, t0, b = x
    model = TransitModel('b', per=per, t0=t0, b=b, rhos=10.)(star.time)
    like, d, vard = star.lnlike(model, full_output=True)
    ll += like
    return ll, (d,)


# Initialize the everest model
star = Everest(201635569)

# Set up the MCMC sampler
params = ['Period (days)', 't0 (BJD - 2456811)', 'Impact parameter']
blobs = ['Depth']
nsteps = 100
nburn = 30
nwalk = 10
ndim = len(params)
nblobs = len(blobs)
sampler = emcee.EnsembleSampler(nwalk, ndim, lnlike, args=[star])
x0 = [[8.368 + 0.01 * np.random.randn(),
       1978.4513 + 0.01 * np.random.randn(),
       0. + 0.1 * np.random.randn()] for k in range(nwalk)]
blobs0 = [[0.] for k in range(nwalk)]

# Run!
for i in tqdm(sampler.sample(x0, iterations=nsteps, blobs0=blobs0),
              total=nsteps):
        pass

# Re-scale the transit time for prettier axes labels
sampler.chain[:, :, 1] -= 1978.

# Take the absolute value of the impact parameter for plotting
sampler.chain[:, :, 2] = np.abs(sampler.chain[:, :, 2])

# Add the blobs to the chain for plotting
chain = np.concatenate((sampler.chain,
                        np.array(sampler.blobs).swapaxes(0, 1)), axis=2)

# Plot the chains
fig1, ax = pl.subplots(ndim + nblobs, figsize=(6, 6))
ax[-1].set_xlabel("Iteration", fontsize=14)
for n in range(ndim + nblobs):
    for k in range(nwalk):
        ax[n].plot(chain[k, :, n], alpha=0.3, lw=1)
    ax[n].set_ylabel((params + blobs)[n], fontsize=10)
    ax[n].margins(0, None)
    ax[n].axvline(nburn, color='b', alpha=0.5, lw=1, ls='--')

# Plot the posterior distributions
samples = chain[:, nburn:, :].reshape(-1, ndim + nblobs)
fig2 = corner(samples, labels=params + blobs)
fig2.set_size_inches(6, 6)
pl.show()
