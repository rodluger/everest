import everest
import numpy as np
import matplotlib.pyplot as pl
from itertools import combinations_with_replacement as multichoose
import george
import os

EPIC = 201862715

if not os.path.exists('%s.npz' % EPIC):
  data = everest.missions.k2.GetData(EPIC, cadence = 'sc')
  time = np.delete(data.time, np.append(data.nanmask, data.badmask))
  fpix = np.delete(data.fpix, np.append(data.nanmask, data.badmask), axis = 0)
  fpix_err = np.delete(data.fpix_err, np.append(data.nanmask, data.badmask), axis = 0)
  fraw = np.sum(fpix, axis = 1)
  ferr = np.sqrt(np.sum(fpix_err ** 2, axis = 1))
  inds = np.where((np.abs(time - 2003.84) < 1.5))
  time = time[inds]
  fpix = fpix[inds]
  fraw = fraw[inds]
  ferr = ferr[inds]
  np.savez('%s.npz' % EPIC, time = time, fpix = fpix, fraw = fraw, ferr = ferr)

data = np.load('%s.npz' % EPIC)
time = data['time']
fpix = data['fpix']
fraw = data['fraw']
ferr = data['ferr']
mask = np.where((time > 2003.78) & (time < 2003.92))[0]

#time = np.delete(time, mask)
#fpix = np.delete(fpix, mask, axis = 0)
#fraw = np.delete(fraw, mask)
#ferr = np.delete(ferr, mask)

mask = np.where(time < 2003.)[0]

def TV(x):
  '''
  
  '''
  
  return np.sum(np.abs(np.diff(x)))

def X(i, j = slice(None, None, None)):
  '''
  
  '''

  X1 = fpix[j] / fraw[j].reshape(-1, 1)
  return np.product(list(multichoose(X1.T, i + 1)), axis = 1).T

def compute(lam):
  '''
  
  '''

  # Covariance
  K = np.diag(np.delete(ferr, mask) ** 2)
  gp = george.GP(1000. ** 2 * george.kernels.Matern32Kernel(3. ** 2))
  K += gp.get_matrix(np.delete(time, mask))

  # The X^2 matrices
  A = np.zeros((len(fraw) - len(mask), len(fraw) - len(mask)))
  B = np.zeros((len(fraw), len(fraw) - len(mask)))
  
  # Loop over all orders
  for n in range(len(lam)):
    XA = X(n, np.delete(np.arange(len(time)), mask))
    XB = X(n)
    A += lam[n] * np.dot(XA, XA.T)
    B += lam[n] * np.dot(XB, XA.T)
  
  W = np.linalg.solve(K + A, np.delete(fraw, mask) - np.nanmedian(np.delete(fraw, mask)))
  model = np.dot(B, W)
  model -= np.nanmedian(model)
  
  return model


lambda_arr = np.array([1,5,7,8,9,10,12,14,15,16], dtype = float)
tv = np.zeros_like(lambda_arr)
for i, l in enumerate(lambda_arr):
  print(l)
  model = compute([10**l])
  tv[i] = TV((fraw - model)[mask])
  
pl.plot(lambda_arr, tv, 'ro')
pl.plot(lambda_arr, tv, 'r-', alpha = 0.5)
pl.axhline(TV(fraw[mask]), color = 'k')
pl.show()
quit()


l = 14
model = compute([10**l, 10**l, 10**l])
fig, ax = pl.subplots(2, figsize = (13, 7), sharex = True, sharey = True)
ax[0].plot(time, fraw, 'k.', alpha = 0.3)
ax[1].plot(time, fraw - model, 'k.', alpha = 0.3)
pl.show()