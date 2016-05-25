# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["KBJD_ZERO", "EXPOSURE_TIMES"]

# The zero point of the Kepler time basis (KBJD).
KBJD_ZERO = 2454833.0

# The (approximate) length of a Kepler exposure in seconds for short and long
# cadence measurements.
EXPOSURE_TIMES = [54.2, 1626.0]
