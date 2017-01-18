Overview
========

**EVEREST** (EPIC Variability Extraction and Removal for Exoplanet Science Targets) 
is an open-source pipeline for removing instrumental noise from light curves. 
It was originally developed for the :py:obj:`K2` mission, but will soon support
additional missions such as :py:obj:`Kepler` and :py:obj:`TESS`.
**EVEREST** exploits correlations across the pixels on the CCD to remove 
systematics introduced by the spacecraft's pointing error. For :py:obj:`K2`, it
yields light curves with precision comparable to that of the original Kepler mission. 
Here we provide 
detailed documentation of the **EVEREST** pipeline, which was coded in Python and is 
freely available on `github <https://github.com/rodluger/everest>`_.

.. toctree::
   :maxdepth: 1
   
   installation
   using_everest
   pipeline

.. raw:: html

  <script>
    (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
    (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
    m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
    })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

    ga('create', 'UA-47070068-2', 'auto');
    ga('send', 'pageview');

  </script>