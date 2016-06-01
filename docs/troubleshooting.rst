Troubleshooting
===============

.. contents::
   :local:
    
I find that :py:mod:`everest` decreases transit depths!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
   This can happen if the transit isn't properly masked during the de-trending step, since
   PLD is trying as hard as it can to remove short-timescale features from the data. The way
   around this is to explicitly mask any known transits during de-trending. This can be
   done by calling :py:func:`everest.usertools.Detrend` and specifying a mask:
   
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
     :width: 600px
     :align: center
     :height: 100px
     :figclass: align-center
  
   The :py:obj:`res` object returned by :py:func:`everest.usertools.Detrend` is a
   tuple containing the time array, the de-trended flux array, and the indices
   of the light curve that were masked.