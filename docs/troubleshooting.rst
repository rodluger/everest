Troubleshooting
===============

.. contents::
   :local:

The :py:mod:`george` package keeps raising errors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:mod:`everest` was coded using :py:mod:`george 0.2.1`. If you have the development
version of :py:mod:`george` installed, you'll get tons of errors. Consider installing
:py:mod:`george 0.2.1` in a separate directory and adding it to your :py:obj:`$PATH`
before calling :py:mod:`everest`.
    
I find that :py:mod:`everest` decreases transit depths!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
This can happen if the transit isn't properly masked during the de-trending step, since
PLD is trying as hard as it can to remove short-timescale features from the data. The way
around this is to explicitly mask any known transits during de-trending. This can be
done by instantiating a :py:obj:`everest.Everest` object and specifying a transit mask:

.. code-block:: python
  
  from everest import Everest
  import matplotlib.pyplot as pl

  # Instantiate to download the target's FITS file
  star = Everest(205071984)

  # Specify a custom mask with the period, time
  # of first transit, and transit duration for each of
  # the planets transiting the target star
  star.set_mask(transits = [( 8.992, 2067.93, 0.25),
                            (20.661, 2066.42, 0.25),
                            (31.716, 2070.79, 0.25)])

  # Plot the de-trended light curve
  star.plot()

  # Plot the whitened, folded light curve
  star.plot_folded()

  # Show the plots
  pl.show()

This should result in the following plots:
  
.. figure:: usertools.jpg
   :width: 700px
   :align: center
   :height: 100px
   :figclass: align-center

I get errors when :py:mod:`everest` tries to plot things
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Consider using a different :py:mod:`matplotlib`
`backend <http://matplotlib.org/api/index_backend_api.html>`_. You can specify
this by calling

.. code-block:: python
  
  import matplotlib
  matplotlib.use('backend_name')

at the top of your script or by editing your ``.matplotlibrc`` file. I've had the 
most success with the ``QtAgg`` backend.
If you're running :py:mod:`everest` non-interactively on a cluster, you might want
to choose a non-interactive backend such as ``agg``.
