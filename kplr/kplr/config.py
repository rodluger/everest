# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["KPLR_ROOT"]

import os

KPLR_ROOT = os.path.expanduser(os.environ.get("KPLR_DATA_DIR",
                                              os.path.join("~", ".kplr")))
