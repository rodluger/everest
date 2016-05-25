import os
import matplotlib.pyplot as pl
import numpy as np

path = '/Users/nks1994/Documents/Research/everest/output/C01/201384232/'
# path = '/Users/nks1994/Documents/Research/everest/output/C01/201367065/'
# path = '/Users/nks1994/Documents/Research/everest/output/C01/201208431/'


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

rms = [rms01, rms02, rms03, rms04, rms05, np.nan, np.nan, np.nan, np.nan]

rms_fit = [rms11, rms12, rms13, rms14, rms15, rms16, rms17, rms18, rms19]

ap = [11,12,13,14,15,16,17,18,19]

pl.title('EPIC 201384232')
pl.xlabel('Aperture Number')
pl.ylabel('RMS')

print(rms,rms_fit)
pl.plot(ap, rms_fit, 'k',marker='o');
pl.plot(ap, rms, 'r',marker='o');
# pl.show()

for n,i,x in zip(ap,rms,rms_fit):
	print("{:>5d} {:>10.5f} {:>10.5f}".format(n,i,x))

# import everest
# everest.RunSingle(<K2 ID number>)