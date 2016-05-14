import os
import matplotlib.pyplot as pl
import numpy as np

data = np.load('/Users/nks1994/Documents/Research/everest/output/C01/201367065/default/data.npz')

r4 = [0,0,0,0,0,data['rms'][3],0,0,0,0]
ap = [10,11,12,13,14,15,16,17,18,19]

pl.xlabel('Aperture Number')
pl.ylabel('RMS')

pl.plot(ap, r4, 'k.'); pl.show()

# import everest
# everest.RunSingle(<K2 ID number>)