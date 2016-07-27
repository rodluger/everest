:py:mod:`kwargs.py` - User options
----------------------------------

This file is located in ``~/.everest/`` and contains the kwargs passed to 
:py:func:`everest.compute.Compute()` when
running a cluster job with :py:func:`everest.run.RunCampaign()`.
See the :py:func:`everest.compute.Compute()` documentation for info on each
of these values.

.. code-block:: python
   
  from __future__ import division, print_function, absolute_import, unicode_literals
  import logging
  import numpy as np

  kwargs = dict(
                  run_name = 'default',
                  apnum = 15, 
                  outlier_sigma = 5,
                  mask_times = [], 
                  pld_order = 3,
                  optimize_npc = True,
                  ps_iter = 30, 
                  ps_masks = 10, 
                  npc_arr = np.arange(25, 260, 10),
                  scatter_alpha = 0.,
                  gp_iter = 2,
                  inject = {}, 
                  fig_ext = 'jpg',
                  jpeg_quality = 30,
                  log_level = logging.DEBUG, 
                  screen_level = logging.CRITICAL
               )

.. raw:: html

  <script>
    (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
    (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
    m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
    })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

    ga('create', 'UA-47070068-2', 'auto');
    ga('send', 'pageview');
  </script>