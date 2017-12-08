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
   masksolve
   math
   pool
   standalone
   transit
   utils

.. raw:: html

  <script>
    (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
    (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
    m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
    })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

    ga('create', 'UA-47070068-3', 'auto');
    ga('send', 'pageview');

  </script>
