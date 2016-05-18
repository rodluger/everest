Troubleshooting
===============

**I find that** :py:mod:`everest` **decreases transit depths!**
   
   This can happen if the transit ins't properly masked during the de-trending step, since
   PLD is trying as hard as it can to remove short-timescale features from the data. The way
   around this is to explicitly mask any known transits during de-trending. *more soon*