import matplotlib.pyplot as pl
from matplotlib.ticker import FuncFormatter
import numpy as np


    

fig, ax = pl.subplots(1)

ax.plot(np.linspace(0,100.,10), np.linspace(1,20,10))

ax.get_yaxis().set_major_formatter(FluxFormatter)
pl.show()