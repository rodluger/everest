Troubleshooting
===============

.. contents::
   :local:
    
I find that :py:mod:`everest` decreases transit depths!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
   This can happen if the transit isn't properly masked during the de-trending step, since
   PLD is trying as hard as it can to remove short-timescale features from the data. The way
   around this is to explicitly mask any known transits during de-trending. This can be
   done by calling :py:func:`everest.tools.Detrend` and specifying a mask:
   
   .. code-block:: python
      
      from everest import usertools as eu
      
      # Create a custom mask
      mask = eu.Mask(transits = [( 8.992, 2067.93, 0.25),
                                 (20.661, 2066.42, 0.25),
                                 (31.716, 2070.79, 0.25)])
      
      # Detrend with this mask
      time, flux, mask_inds = eu.Detrend(205071984, mask = mask)
      
      
      