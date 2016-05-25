# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import unittest

try:
    from unittest import mock
except ImportError:
    import mock

from kplr.api import API, Model, KOI, Planet, Star, LightCurve
from kplr.config import KPLR_ROOT


class ApiTestCase(unittest.TestCase):
    def test_default_data_root(self):
        api = API()
        self.assertEqual(api.data_root, KPLR_ROOT)

    def test_custom_data_root(self):
        api = API("/home/data/")
        self.assertEqual(api.data_root, "/home/data/")

    def test_data_root_in_str_and_repr(self):
        api = API()
        self.assertIn(api.data_root, str(api))
        self.assertIn(api.data_root, repr(api))

    def test_munge_dict_int_value(self):
        api = API()
        row = {"key": "666"}
        new_row = api._munge_dict(row)
        self.assertEqual(new_row["key"], 666)

    def test_munge_dict_float_value(self):
        api = API()
        row = {"key": "66.6"}
        new_row = api._munge_dict(row)
        self.assertAlmostEqual(new_row["key"], 66.6)

    def test_munge_dict_text_value(self):
        api = API()
        row = {"key": "value"}
        new_row = api._munge_dict(row)
        self.assertEqual(new_row["key"], "value")

    def test_munge_dict_empty_value(self):
        api = API()
        row = {"key": ""}
        new_row = api._munge_dict(row)
        self.assertIsNone(new_row["key"])


class TestModel(Model):
    _id = "\"{kepler_name}\""


class ModelTestCase(unittest.TestCase):
    def setUp(self):
        self.mock_api = mock.MagicMock(spec=API)
        self.params = {
            "kepler_name": "Kepler-32 f",
            "kepid": 9787239,
        }
        self.model = TestModel(self.mock_api, self.params)

    def test_setting_params_in_init(self):
        self.assertEqual(self.model.kepler_name, self.params["kepler_name"])

    def test_class_in_str_and_repr(self):
        self.assertIn("TestModel", str(self.model))
        self.assertIn("TestModel", repr(self.model))

    def test_name_in_str_and_repr(self):
        self.assertIn("Kepler-32 f", str(self.model))
        self.assertIn("Kepler-32 f", repr(self.model))

    def test_get_light_curves(self):
        self.model.get_light_curves(short_cadence=False, fetch=True,
                                    clobber=True)
        self.mock_api.light_curves.assert_called_once_with(
            self.model.kepid,
            short_cadence=False,
            fetch=True,
            clobber=True,
        )

    def test_get_target_pixel_files(self):
        self.model.get_target_pixel_files(short_cadence=False, fetch=True,
                                          clobber=True)
        self.mock_api.target_pixel_files.assert_called_once_with(
            self.model.kepid,
            short_cadence=False,
            fetch=True,
            clobber=True,
        )


class KOIModelTestCase(unittest.TestCase):
    def setUp(self):
        self.mock_api = mock.MagicMock(spec=API)
        self.mock_api.star.return_value = mock.MagicMock(spec=Star)
        self.params = {
            "kepid": "9787239",
            "kepoi_name": "K00952.01",
        }
        self.koi = KOI(self.mock_api, self.params)

    def test_kepoi_name_in_str_and_repr(self):
        self.assertIn(self.params["kepoi_name"], str(self.koi))

    def test_star_property_is_cached(self):
        self.assertIsNotNone(self.koi.star)
        self.assertTrue(self.mock_api.star.called)
        self.mock_api.star.reset_mock()
        self.assertIsNotNone(self.koi.star)
        self.assertFalse(self.mock_api.star.called)


class PlanetModelTestCase(unittest.TestCase):
    def setUp(self):
        self.mock_api = mock.MagicMock(spec=API)
        self.mock_api.koi.return_value = mock.MagicMock(spec=KOI)
        self.mock_api.star.return_value = mock.MagicMock(spec=Star)
        self.params = {
            "kepid": "9787239",
            "kepler_name": "Kepler-32 f",
            "koi_number": "952.05",
        }
        self.planet = Planet(self.mock_api, self.params)

    def test_kepler_name_in_str_and_repr(self):
        self.assertIn(self.params["kepler_name"], str(self.planet))

    def test_koi_property_is_cached(self):
        self.assertIsNotNone(self.planet.koi)
        self.assertTrue(self.mock_api.koi.called)
        self.mock_api.koi.reset_mock()
        self.assertIsNotNone(self.planet.koi)
        self.assertFalse(self.mock_api.koi.called)

    def test_star_property_is_cached(self):
        self.assertIsNotNone(self.planet.star)
        self.assertTrue(self.mock_api.star.called)
        self.mock_api.star.reset_mock()
        self.assertIsNotNone(self.planet.star)
        self.assertFalse(self.mock_api.star.called)


class StarModelTestCase(unittest.TestCase):
    def setUp(self):
        self.mock_api = mock.MagicMock(spec=API)
        self.mock_api.kois.return_value = [mock.MagicMock(spec=KOI)]
        self.params = {
            "kic_kepler_id": "9787239",
            "kepid": "9787239",
        }
        self.star = Star(self.mock_api, self.params)

    def test_kic_kepler_id_in_str_and_repr(self):
        self.assertIn(self.params["kic_kepler_id"], str(self.star))

    def test_kois_property_is_cached(self):
        self.assertIsNotNone(self.star.kois)
        self.assertTrue(self.mock_api.kois.called)
        self.mock_api.kois.reset_mock()
        self.assertIsNotNone(self.star.kois)
        self.assertFalse(self.mock_api.kois.called)


class LightCurveTestCase(unittest.TestCase):
    def setUp(self):
        self.mock_api = mock.MagicMock(spec=API)
        self.mock_api.data_root = "/home/data/"
        self.params = {
            "sci_data_set_name": "KPLR009787239-2009166043257",
            "ktc_kepler_id": "9787239",
            "ktc_target_type": "LC",
        }
        self.lightcurve = LightCurve(self.mock_api, self.params)

    def test_filename(self):
        filename = self.lightcurve.filename
        self.assertTrue(filename.endswith(".fits"))
        self.assertIn(self.mock_api.data_root, filename)
        self.assertIn(self.lightcurve.sci_data_set_name.lower(), filename)

    def test_url(self):
        url = self.lightcurve.url
        self.assertTrue(url.endswith(".fits"))
        self.assertIn("lightcurves", url)
        self.assertIn(self.lightcurve.kepid, url)
