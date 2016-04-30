#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
kwargs.py
---------

This file contains the kwargs passed to ``everest.Compute()`` when
running a cluster job with ``everest.RunCampaign()``.

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
                ps_iter = 50, 
                ps_masks = 10, 
                npc_arr = np.arange(25, 200, 10),
                scatter_alpha = 0.,
                gp_iter = 2,
                inject = {}, 
                log_level = logging.DEBUG, 
                screen_level = logging.CRITICAL
             )