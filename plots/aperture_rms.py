import os
import matplotlib.pyplot as pl
import numpy as np

path = '/Users/nks1994/Documents/Research/everest/output/C01/201367065/'
rms00 = np.load(path + 'ap00/data.npz')['rms'][3]
rms01 = np.load(path + 'ap01/data.npz')['rms'][3]
rms02 = np.load(path + 'ap02/data.npz')['rms'][3]
rms03 = np.load(path + 'ap03/data.npz')['rms'][3]
rms04 = np.load(path + 'ap04/data.npz')['rms'][3]
rms05 = np.load(path + 'ap05/data.npz')['rms'][3]
rms10 = np.load(path + 'ap10/data.npz')['rms'][3]
rms11 = np.load(path + 'ap11/data.npz')['rms'][3]
rms12 = np.load(path + 'ap12/data.npz')['rms'][3]
rms13 = np.load(path + 'ap13/data.npz')['rms'][3]
rms14 = np.load(path + 'ap14/data.npz')['rms'][3]
rms15 = np.load(path + 'ap15/data.npz')['rms'][3]
rms16 = np.load(path + 'ap16/data.npz')['rms'][3]
rms17 = np.load(path + 'ap17/data.npz')['rms'][3]
rms18 = np.load(path + 'ap18/data.npz')['rms'][3]
rms19 = np.load(path + 'ap19/data.npz')['rms'][3]

rms = [rms00, rms01, rms02, rms03, rms04, rms05, None, None, None, None]

rms_fit = [rms10, rms11, rms12, rms13, rms14, rms15, rms16, rms17, rms18, rms19]

ap = [10,11,12,13,14,15,16,17,18,19]

pl.xlabel('Aperture Number')
pl.ylabel('RMS')

print(rms,rms_fit)
pl.plot(ap, rms_fit, 'k',marker='o'); pl.plot(ap, rms, 'r',marker='o'); pl.show()

# import everest
# everest.RunSingle(<K2 ID number>)