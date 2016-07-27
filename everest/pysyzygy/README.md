py·sy·zy·gy
-----------
/ˈpīsizijē/
<a href="https://raw.githubusercontent.com/rodluger/pysyzygy/master/img/pysyzygy.mp3"><img style="float: right;" src="img/speaker.png?raw=True"/></a>
<a href="https://raw.githubusercontent.com/rodluger/pysyzygy/master/LICENSE"><img align="right" src="https://img.shields.io/badge/license-MIT-blue.svg"/></a>
<a href="https://coveralls.io/github/rodluger/pysyzygy?branch=master"><img align="right" src="https://coveralls.io/repos/github/rodluger/pysyzygy/badge.svg?branch=master"/></a>
<a href="https://travis-ci.org/rodluger/pysyzygy"><img align="right" src="https://travis-ci.org/rodluger/pysyzygy.svg?branch=master"/></a>

**1.** A fast and general planet transit ([syzygy](http://en.wikipedia.org/wiki/Syzygy_%28astronomy%29)) code written in C and in Python.

**2.** ``pysyzygy`` computes **fast** lightcurves for the most general case of a *massive*, *eccentric* planet orbiting a limb-darkened star. Here's a sample output image of an assymetric transit:

![transit](img/transit.png?raw=True)

Installation
============
Clone the repository and run

```bash
python setup.py install
```

Calling pysyzygy...
===================

... is super easy.

```python
import pysyzygy as ps
import numpy as np
import matplotlib.pyplot as pl

# Instantiate a transit object
trn = ps.Transit(t0 = 0.5, RpRs = 0.1, per = 1.234) 

# Now evaluate the light curve on a grid of observation times
t = np.arange(0., 3.5, ps.transit.KEPSHRTCAD)
flux = trn(t)

# Plot the light curve
fig, ax = pl.subplots(figsize = (12, 5))
fig.subplots_adjust(bottom = 0.15)
ax.plot(t, flux, 'b.')
ax.set_xlabel('Time (days)', fontsize = 18)
ax.set_ylabel('Relative flux', fontsize = 18)
ax.margins(None, 0.5)
pl.show()
```     

![transit](img/hotjup.png?raw=True)

Notes
=====

More detailed documentation coming soon. For now, check out the [examples](examples) directory for
some cool things you can do with ``pysyzygy``.

Feel free to change, adapt, or incorporate this code into your project, but please make sure to cite this repository, as well as [Mandel and Agol (2002)](http://adsabs.harvard.edu/abs/2002ApJ...580L.171M>), the transit model on which ``pysyzygy`` is based.
