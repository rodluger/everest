import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
from everest2.config import EVEREST_DAT
from everest2.math import Chunks
from everest2.missions.k2.aux import GetK2Campaign, GetNeighboringChannels
from scipy.signal import savgol_filter

# User
channels = range(99)
campaign = 6.0
model = 'SegmentedPLD'

# Get the light curves
all = GetK2Campaign(campaign)
stars = np.array([s[0] for s in all if s[2] in channels and 
        os.path.exists(os.path.join(EVEREST_DAT, 'k2', 'c%02d' % int(campaign),
        ('%09d' % s[0])[:4] + '00000', 
        ('%09d' % s[0])[4:], model + '.npz'))], dtype = int)
cdpp6 = []
for n in range(len(stars)):
  nf = os.path.join(EVEREST_DAT, 'k2', 'c%02d' % int(campaign),
                 ('%09d' % stars[n])[:4] + '00000', 
                 ('%09d' % stars[n])[4:], model + '.npz')
  data = np.load(nf)
  
  if (data['cdpp6'] >= 10) and (data['cdpp6'] <= 30):
    
    t = data['time']
    y = data['fraw'] - data['model'] 
    m = np.array(list(set(np.concatenate([data['outmask'], data['badmask'], data['nanmask']]))), dtype = int)
    y = np.interp(t, np.delete(t, m), np.delete(y, m))
    y /= np.nanmedian(y)
    pl.plot(t, y, 'k-', alpha = 0.05)
    
    

pl.show()