Code Engine
===========

The modules below house all of the mission-independent routines for 
de-trending light curves. The core of the code is in :doc:`basecamp <basecamp>` and
:doc:`detrender <detrender>`. The classes implemented in these two modules call the
rest of the code below. For mission-specific routines, which include downloading
and pre-processing the raw data and formatting the final output, see :doc:`missions <missions>`.

.. toctree::
   :maxdepth: 2
   
   basecamp
   config
   detrender
   dvs
   fits
   gp
   inject
   math
   pool
   transit
   utils