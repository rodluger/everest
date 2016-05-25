# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import unittest

try:
    from unittest import mock
except ImportError:
    import mock

from kplr.mast import Adapter


class AdapterTestCase(unittest.TestCase):
    def test_default_conversion_success(self):
        adapter = Adapter({})
        row = {"Ang Sep (')": "0.55"}
        output_row = adapter(row)
        self.assertAlmostEqual(output_row["angular_separation"], 0.55)

    def test_default_conversion_failure(self):
        adapter = Adapter({})
        row = {"Ang Sep (')": "notanumber"}
        output_row = adapter(row)
        self.assertIsNone(output_row["angular_separation"])

    def test_custom_conversion(self):
        adapter = Adapter({
            "Kepler ID": ("kepid", int),
        })
        row = {"Kepler ID": "666"}
        output_row = adapter(row)
        self.assertEqual(output_row["kepid"], 666)

    @mock.patch('kplr.mast.logging')
    def test_unknown_parameter_logs_warning(self, mock_logging):
        adapter = Adapter({
            "Kepler ID": ("kepid", int),
        })
        row = {"Some other parameter": "42"}
        adapter(row)
        mock_logging.warn.assert_called_once()
