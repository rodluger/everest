Troubleshooting
===============

**The function** :py:func:`everest.plot.Plot` **spits out an error when trying to save an image to disk**
   
   If you're on a Mac, try changing your :py:mod:`matplotlib` backend to something other
   than `MacOSX`, or try changing the `figext` keyword in `kwargs.py` to something other
   than `jpg`. I've had lots of trouble with the `MaxOSX` backend lately...

**I find that** :py:mod:`everest` **decreases transit depths!**
   
   This can happen if the transit ins't properly masked during the de-trending step, since
   PLD is trying as hard as it can to remove short-timescale features from the data. The way
   around this is to explicitly mask any known transits during de-trending. *more soon*