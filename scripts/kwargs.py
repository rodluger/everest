#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
kwargs.py
---------

This file contains the kwargs passed to :py:func:`everest.compute.Compute()` when
running a cluster job with :py:func:`everest.run.RunCampaign()`.
See the :py:func:`everest.compute.Compute()` documentation for info on each
of these values.

.. literalinclude:: ../scripts/kwargs.py
   :lines: 17-38
   
'''

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