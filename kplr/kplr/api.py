# -*- coding: utf-8 -*-
'''
:py:mod:`api.py` - The kplr API
-------------------------------

This :py:mod:`kplr` package was forked from `<https://github.com/dfm/kplr>`_ and
modified to allow easy access to different `K2` databases.
Read the original documentation `here <http://dan.iel.fm/kplr/>`_.

'''

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["API", "KOI", "Planet", "Star", "LightCurve", "TargetPixelFile"]

import os
import re
import json
import shutil
import logging
from itertools import product
from functools import partial
from tempfile import NamedTemporaryFile
import random

import six
from six.moves import urllib

# Optional dependencies.
try:
    import pyfits
except ImportError:
    try:
        import astropy.io.fits as pyfits
    except ImportError:
        pyfits = None

try:
    import fitsio
except ImportError:
    fitsio = None

try:
    import numpy as np
except ImportError:
    np = None

try:
    import matplotlib.pyplot as pl
except ImportError:
    pl = None

try:
    from tornado import gen
    import tornado.ioloop
    from tornado.httpclient import AsyncHTTPClient, HTTPRequest
except ImportError:
    AsyncHTTPClient = None

    def async_download(*args, **kwargs):
        raise ImportError("You need to install tornado to do an async "
                          "download")
else:
    def handle_response(obj, response):
        if response.error:
            response.rethrow()
        obj._save_fetched_file(response.body)

    @gen.coroutine
    def _async_download(lcs, clobber=False):
        client = AsyncHTTPClient(max_clients=10)
        to_fetch = [l for l in lcs if clobber or not l.cache_exists]
        responses = yield [client.fetch(HTTPRequest(l.url))
                           for l in to_fetch]
        [handle_response(to_fetch[i], r) for i, r in enumerate(responses)]
        raise gen.Return()

    def async_download(lcs, clobber=False):
        f = partial(_async_download, lcs, clobber=clobber)
        ioloop = tornado.ioloop.IOLoop.instance()
        ioloop.run_sync(f)

from .config import KPLR_ROOT
from . import mast

class API(object):
    """
    Interface with MAST and Exoplanet Archive APIs.

    :param data_root: (optional)
        The local base directory where any data should be downloaded to. This
        can also be set using the ``KPLR_ROOT`` environment variable. The
        default value is ``~/.kplr``.

    """

    ea_url = ("http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI"
              "/nph-nstedAPI")
    mast_url = "http://archive.stsci.edu/{0}/{1}/search.php"

    def __init__(self, data_root=None):
        self.data_root = data_root
        if data_root is None:
            self.data_root = KPLR_ROOT

    def __str__(self):
        return "<API(data_root=\"{0}\")>".format(self.data_root)

    def __unicode__(self):
        return self.__str__()

    def __repr__(self):
        return self.__str__()

    def ea_request(self, table, sort=None, **params):
        """
        Submit a request to the Exoplanet Archive API and return a dictionary.

        :param table:
            The table that you want to search.

        :param params:
            Any other search parameters.

        """
        params["table"] = table

        # Deal with sort order.
        if sort is not None:
            if isinstance(sort, six.string_types):
                params["order"] = sort
            else:
                params["order"] = sort[0]
                if sort[1] == -1:
                    params["order"] += "+desc"

        # Format the URL in the *horrifying* way that EA needs it to be...
        # they don't un-escape the HTTP parameters!!!
        payload = ["{0}={1}".format(k, urllib.parse.quote_plus(v, "\"'+"))
                   for k, v in params.items()]

        # Send the request.
        encoded_data = "&".join(payload).encode("ascii")
        r = urllib.request.Request(self.ea_url, data=encoded_data)
        handler = urllib.request.urlopen(r)
        code = handler.getcode()
        txt = handler.read().decode("ascii")

        # Hack because Exoplanet Archive doesn't return HTTP errors.
        if int(code) != 200 or "ERROR" in txt:
            full_url = handler.geturl() + "?" + "&".join(payload)
            raise APIError(code, full_url, txt)

        # Parse the CSV output.
        csv = txt.splitlines()
        columns = csv[0].split(",")
        result = []
        for line in csv[1:]:
            result.append(dict(zip(columns, line.split(","))))

        return [self._munge_dict(row) for row in result]

    def mast_request(self, category, adapter=None, sort=None, mission="kepler",
                     **params):
        """
        Submit a request to the MAST API and return a dictionary of parameters.

        :param category:
            The table that you want to search.

        :param params:
            Any other search parameters.

        """
        params["action"] = params.get("action", "Search")
        params["outputformat"] = "JSON"
        params["coordformat"] = "dec"
        if params.get("selectedColumnsCsv", None) is not None:
          params["verb"] = 3
        
        # Deal with sort order.
        if sort is not None:
            if isinstance(sort, six.string_types):
                params["ordercolumn1"] = sort
            else:
                params["ordercolumn1"] = sort[0]
                if sort[1] == -1:
                    params["descending1"] = "on"

        # Send the request.
        encoded_data = urllib.parse.urlencode(params).encode("ascii")
        r = urllib.request.Request(self.mast_url.format(mission, category),
                                   data=encoded_data)
        handler = urllib.request.urlopen(r)
        code = handler.getcode()
        txt = handler.read().decode("ascii")
        if int(code) != 200:
            full_url = handler.geturl() + "?" + urllib.parse.urlencode(params)
            raise APIError(code, full_url, txt)

        # Check for no rows found.
        if "no rows" in txt:
            return []

        # Parse the JSON.
        try:
            result = json.loads(txt)
        except ValueError:
            full_url = handler.geturl() + "?" + urllib.parse.urlencode(params)
            raise APIError(code, full_url,
                           "No JSON object could be decoded.\n" + txt)

        # Fake munge the types if no adapter was provided.
        if adapter is None:
            return [self._munge_dict(row) for row in result]

        return [adapter(row) for row in result]

    def _munge_dict(self, row):
        """
        Iterate through a dictionary and try to interpret the data types in a
        sane way.

        :param row:
            A dictionary of (probably) strings.

        :returns new_row:
            A dictionary with the same keys as ``row`` but reasonably typed
            values.

        """
        tmp = {}
        for k, v in row.items():
            # Good-god-what-type-is-this-parameter?!?
            try:
                tmp[k] = int(v)
            except ValueError:
                try:
                    tmp[k] = float(v)
                except ValueError:
                    tmp[k] = v

            # Empty entries are mapped to None.
            if v == "":
                tmp[k] = None

        return tmp

    def kois(self, **params):
        """
        Get a list of KOIs from The Exoplanet Archive.

        :param params:
            The search parameters for the Exoplanet Archive API.

        """
        params["select"] = params.get("select", "*")
        return [KOI(self, k) for k in self.ea_request("cumulative", **params)]

    def koi(self, koi_number, **params):
        """
        Find a single KOI given a KOI number (e.g. 145.01).

        :param koi_number:
            The number identifying the KOI. This should be a float with the
            ``.0N`` for some value of ``N``.

        """
        kois = self.kois(where="kepoi_name+like+'K{0:08.2f}'"
                         .format(float(koi_number)), **params)
        if not len(kois):
            raise ValueError("No KOI found with the number: '{0}'"
                             .format(koi_number))
        return kois[0]

    def planets(self, **params):
        """
        Get a list of confirmed (Kepler) planets from MAST.

        :param params:
            The search parameters for the MAST API.

        """
        planets = self.mast_request("confirmed_planets",
                                    adapter=mast.planet_adapter, **params)
        return [Planet(self, p) for p in planets]

    def planet(self, name):
        """
        Get a planet by the Kepler name (e.g. "6b" or "Kepler-62b").

        :param name:
            The name of the planet.

        """
        # Parse the planet name.
        matches = re.findall("([0-9]+)[-\s]*([a-zA-Z])", name)
        if len(matches) != 1:
            raise ValueError("Invalid planet name '{0}'".format(name))
        kepler_name = "Kepler-{0} {1}".format(*(matches[0]))

        # Query the API.
        planets = self.planets(kepler_name=kepler_name, max_records=1)
        if not len(planets):
            raise ValueError("No planet found with the name: '{0}'"
                             .format(kepler_name))
        return planets[0]

    def stars(self, **params):
        """
        Get a list of KIC targets from MAST. Only return up to 100 results by
        default.

        :param params:
            The query parameters for the MAST API.

        """
        params["max_records"] = params.pop("max_records", 100)
        stars = self.mast_request("kic10", adapter=mast.star_adapter,
                                  **params)
        return [Star(self, s) for s in stars]

    def star(self, kepid):
        """
        Get a KIC target by id from MAST.

        :param kepid:
            The integer ID of the star in the KIC.

        """
        stars = self.stars(kic_kepler_id=kepid, max_records=1)
        if not len(stars):
            raise ValueError("No KIC target found with id: '{0}'"
                             .format(kepid))
        return stars[0]

    def kepler_star_list(self, quarter = 8, **params):
        """
        Return a list of all the KIC target IDs for which data was collected
        in a given quarter.
        
        """

        params["selectedColumnsCsv"] = "ktc_kepler_id"
        params["ordercolumn1"] = "ktc_kepler_id"
        params["ktc_target_type"] = "LC"
        params["sci_data_quarter"] = str(quarter)
        params["max_records"] = params.pop("max_records", 999999)
        stars = self.mast_request("data_search", adapter=mast.mini_kic_adapter,
                                  mission="kepler", **params)
        return [star["ktc_kepler_id"] for star in stars]
        
    def k2_stars(self, **params):
        """
        Get a list of EPIC targets from MAST. Only return up to 100 results by
        default.

        :param params:
            The query parameters for the MAST API.

        """
        params["max_records"] = params.pop("max_records", 100)
        stars = self.mast_request("epic", adapter=mast.epic_adapter,
                                  mission="k2", **params)
        return [K2Star(self, s) for s in stars]

    def k2_star_list(self, **params):
        """
        Return a dict of all the EPIC target IDs
        for which LC data was collected, indexed by campaign,
        and randomly shuffled.
        
        Returns objects that are listed as either "star" or "\\null".
        
        """

        params["selectedColumnsCsv"] = "ktc_k2_id,sci_campaign"
        params["ktc_target_type"] = "LC"
        params["max_records"] = params.pop("max_records", 999999)
        params["objtype"] = "\\null"
        nulls = self.mast_request("data_search", adapter=mast.mini_dataset_adapter,
                                  mission="k2", **params)
        params["objtype"] = "star"
        stars = self.mast_request("data_search", adapter=mast.mini_dataset_adapter,
                                  mission="k2", **params)
        stars = stars + nulls
        random.shuffle(stars)
        
        # Sort them into campaigns
        c = [[] for i in range(99)]
        for star in stars:
          c[star["sci_campaign"]].append(star["ktc_k2_id"])
        
        # Create a dict
        res = {}
        for campaign in range(99):
          if len(c[campaign]):
            res.update({campaign: c[campaign]})

        return res

    def k2_star_mags(self, stars_per_mag = 50, mags = range(8,18), **params):
        """
        Returns the EPIC numbers of ``stars_per_mag`` random stars in each 
        Kepler magnitude bin in the range ``mags``.
        
        """
        params["selectedColumnsCsv"] = "id"
        params["ordercolumn1"] = "id"
        params["max_records"] = params.pop("max_records", 999999)
        params["k2_avail_flag"] = "1"
        targets = []
        for kp in mags:
          params["kp"] = "%d..%d" % (kp, kp + 1)
          stars = self.mast_request("epic", adapter=mast.mini_epic_adapter,
                                  mission="k2", **params)
          random.shuffle(stars)
          targets.append([star["id"] for star in stars[:stars_per_mag]])
        return targets

    def k2_star(self, id_):
        """
        Get a EPIC target by id from MAST.

        :param id:
            The integer ID of the star in the EPIC.

        """
        stars = self.k2_stars(id=id_, max_records=1)
        if not len(stars):
            raise ValueError("No EPIC target found with id: '{0}'"
                             .format(id_))
        return stars[0]

    def _data_search(self, kepler_id=None, short_cadence=True,
                     adapter=mast.dataset_adapter, **params):
        """
        Run a generic data search on MAST to return a list of dictionaries
        describing the data products.

        :param kepler_id:
            The KIC ID of the target star.

        :param short_cadence:
            A boolean flag that determines whether or not the short cadence
            data should be included.

        :param params:
            Any other search parameters.

        """
        if kepler_id is not None:
            params["ktc_kepler_id"] = kepler_id
        if not short_cadence:
            params["ktc_target_type"] = "LC"

        data_list = self.mast_request("data_search", adapter=adapter,
                                      **params)
        return data_list

    def light_curves(self, kepler_id=None, short_cadence=True, fetch=False,
                     clobber=False, async=False, **params):
        """
        Find the set of light curves associated with a KIC target.

        :param kepler_id:
            The KIC ID of the target star.

        :param short_cadence:
            A boolean flag that determines whether or not the short cadence
            data should be included. (default: True)

        :param fetch:
            A boolean flag that determines whether or not the data file should
            be downloaded.

        :param clobber:
            A boolean flag that determines whether or not the data file should
            be overwritten even if it already exists.

        :param async:
            If ``True``, download the files asynchronously using Tornado.

        :param params:
            Other search parameters to be passed to the MAST data search.

        """
        lcs = [LightCurve(self, d) for d in self._data_search(kepler_id,
               short_cadence=short_cadence, **params)]
        if fetch:
            if async:
                async_download(lcs, clobber=clobber)
            else:
                [l.fetch(clobber=clobber) for l in lcs]
        return lcs

    def target_pixel_files(self, kepler_id=None, short_cadence=True,
                           fetch=False, clobber=False, async=False,
                           cls=None, **params):
        """
        Find the set of target pixel files associated with a KIC target.

        :param kepler_id:
            The KIC ID of the target star.

        :param short_cadence:
            A boolean flag that determines whether or not the short cadence
            data should be included. (default: True)

        :param fetch:
            A boolean flag that determines whether or not the data file should
            be downloaded.

        :param clobber:
            A boolean flag that determines whether or not the data file should
            be overwritten even if it already exists.

        :param async:
            If ``True``, download the files asynchronously using Tornado.

        :param params:
            Other search parameters to be passed to the MAST data search.

        """
        if cls is None:
            cls = TargetPixelFile
        tpfs = [cls(self, d) for d in self._data_search(kepler_id,
                short_cadence=short_cadence, **params)]
        if fetch:
            if async:
                async_download(tpfs, clobber=clobber)
            else:
                [l.fetch(clobber=clobber) for l in tpfs]
        return tpfs


class APIError(Exception):
    """
    Exception raised when an API request fails.

    :param code:
        The HTTP status code that caused the failure.

    :param url:
        The endpoint (with parameters) of the request.

    :param txt:
        A human readable description of the error.

    """

    def __init__(self, code, url, txt):
        super(APIError, self).__init__(("The API returned code {0} for URL: "
                                        "'{1}' with message:\n{2}")
                                       .format(code, url, txt))
        self.code = code
        self.txt = txt
        self.url = url


class Model(object):
    """
    An abstract base class that provides functions for converting the JSON
    dictionaries returned by the API into Python objects with the correct
    properties.

    :param api:
        A reference to the :class:`API` object that made the request.

    :param params:
        The dictionary of values returned by the API.

    """

    # A format string that generates a unique identifier for the model that
    # should be overloaded by subclasses.
    _id = "..."

    def __init__(self, api, params):
        self.api = api
        self.params = params
        for k, v in params.items():
            setattr(self, k, v)

        self._name = self._id.format(**params)

    def __str__(self):
        return "<{0}({1})>".format(self.__class__.__name__, self._name)

    def __unicode__(self):
        return self.__str__()

    def __repr__(self):
        return self.__str__()

    def get_light_curves(self, **kwargs):
        """
        Get a list of light curve datasets for the model and optionally
        download the FITS files. See :func:`API.light_curves` for parameter
        documentation.

        """
        return self.api.light_curves(self.kepid, **kwargs)

    def get_target_pixel_files(self, **kwargs):
        """
        Get a list of target pixel datasets for the model and optionally
        download the FITS files. See :func:`API.target_pixel_files` for
        parameter documentation.

        """
        return self.api.target_pixel_files(self.kepid, **kwargs)


class KOI(Model):
    """
    A model specifying a Kepler Object of Interest (KOI).

    """

    _id = "\"{kepoi_name}\""

    def __init__(self, *args, **params):
        super(KOI, self).__init__(*args, **params)
        self._star = None

    @property
    def star(self):
        """
        The :class:`Star` entry from the Kepler Input Catalog associated with
        this object.

        """
        if self._star is None:
            self._star = self.api.star(self.kepid)
        return self._star


class Planet(Model):
    """
    A confirmed planet from the `MAST confirmed_planets table
    <http://archive.stsci.edu/search_fields.php?mission=kepler_cp>`_. This
    table has far less—and far less accurate—information than the KOI table
    so it's generally a good idea to use the ``koi`` property to access the
    catalog values.

    """

    _id = "\"{kepler_name}\""

    def __init__(self, *args, **params):
        super(Planet, self).__init__(*args, **params)
        self._koi = None
        self._star = None

    @property
    def koi(self):
        """
        The :class:`KOI` entry that led to this planet. The KOI table is much
        more complete so the use of this object tends to be preferred over the
        built in :class:`Planet` property values.

        """
        if self._koi is None:
            self._koi = self.api.koi(self.koi_number)
        return self._koi

    @property
    def star(self):
        """
        The :class:`Star` entry from the Kepler Input Catalog associated with
        this object.

        """
        if self._star is None:
            self._star = self.api.star(self.kepid)
        return self._star


class Star(Model):
    """
    A star from the `Kepler Input Catalog (KIC)
    <http://archive.stsci.edu/search_fields.php?mission=kic10>`_.

    """

    _id = "{kic_kepler_id}"

    def __init__(self, *args, **params):
        super(Star, self).__init__(*args, **params)
        self.kepid = self.kic_kepler_id
        self._kois = None

    @property
    def kois(self):
        """
        The list of :class:`KOI` entries found in this star's light curve.

        """
        if self._kois is None:
            self._kois = self.api.kois(where="kepid like '{0}'"
                                       .format(self.kepid))
        return self._kois


class K2Star(Model):
    """
    A star from the `K2 EPIC Catalog (EPIC)
    <http://archive.stsci.edu/search_fields.php?mission=epic>`_.

    """

    _id = "{id}"

    def get_light_curves(self, **kwargs):
        """
        Get a list of light curve datasets for the model and optionally
        download the FITS files. See :func:`API.light_curves` for parameter
        documentation.

        """
        raise NotImplementedError("There aren't any light curves for K2 stars")

    def get_target_pixel_files(self, **kwargs):
        """
        Get a list of target pixel datasets for the model and optionally
        download the FITS files. See :func:`API.target_pixel_files` for
        parameter documentation.

        """
        return self.api.target_pixel_files(ktc_k2_id=self.id, mission="k2",
                                           adapter=mast.k2_dataset_adapter,
                                           cls=K2TargetPixelFile,
                                           **kwargs)


class _datafile(Model):

    _id = "\"{sci_data_set_name}_{ktc_target_type}\""

    kepid_template = "{ktc_kepler_id:09d}"

    product = None
    suffixes = None
    filetype = None

    def __init__(self, *args, **params):
        super(_datafile, self).__init__(*args, **params)
        self.kepid = self.kepid_template.format(**(self.params))
        suffix = self.suffixes[int(self.ktc_target_type != "LC")]
        self._filename = "{0}_{1}{2}".format(self.sci_data_set_name,
                                             suffix, self.filetype).lower()

    @classmethod
    def local(cls, fn):
        self = cls.__new__(cls)
        self._name = "'local:{0}'".format(fn)
        self.base_dir = ""
        self._filename = fn
        return self

    @property
    def base_dir(self):
        """
        The local base directory for these types of products.

        """
        return os.path.join(self.api.data_root, "data", self.product,
                            self.kepid)

    @property
    def filename(self):
        """
        The local filename of the data file. This file is only guaranteed to
        exist after ``fetch()`` has been called.

        """
        return os.path.join(self.base_dir, self._filename)

    @property
    def url(self):
        """
        The remote URL for the data file on the MAST servers.

        """
        base_url = "http://archive.stsci.edu/pub/kepler/{0}/{1}/{2}/{3}"
        return base_url.format(self.product, self.kepid[:4], self.kepid,
                               self._filename)

    def open(self, clobber=False, **kwargs):
        """
        Open the FITS data file and return the ``pyfits.HDUList``. This will
        download the file if it isn't already saved locally.

        :param clobber:
            Overwrite the local file even if it exists? This can be helpful if
            the file gets corrupted somehow.

        :param kwargs:
            Any keyword arguments that you would like to pass to the
            :func:`pyfits.open` function.

        """
        if pyfits is None:
            raise ImportError("The pyfits module is required to read data "
                              "files.")
        fn = self.filename
        self.fetch(clobber=clobber)
        return pyfits.open(fn, **kwargs)

    def read(self, clobber=False, **kwargs):
        """
        Read a FITS file using `fitsio <https://github.com/esheldon/fitsio>`_.

        :param clobber:
            Overwrite the local file even if it exists? This can be helpful if
            the file gets corrupted somehow.

        :param kwargs:
            Any keyword arguments that you would like to pass to the
            :func:`fitsio.read` function.

        """
        if fitsio is None:
            raise ImportError("The fitsio module is required to read data "
                              "files.")
        fn = self.filename
        self.fetch(clobber=clobber)
        return fitsio.read(fn, **kwargs)

    @property
    def cache_exists(self):
        return os.path.exists(self.filename)

    def fetch(self, clobber=False):
        """
        Download the data file from the server and save it locally. The local
        file will be saved in the directory specified by the ``data_root``
        property of the API.

        :param clobber:
            Should an existing local file be overwritten? (default: False)

        """
        # Check if the file already exists.
        filename = self.filename
        if self.cache_exists and not clobber:
            logging.info("Found local file: '{0}'".format(filename))
            return self

        # Fetch the remote file.
        url = self.url
        logging.info("Downloading file from: '{0}'".format(url))
        r = urllib.request.Request(url)
        handler = urllib.request.urlopen(r)
        code = handler.getcode()
        if int(code) != 200:
            raise APIError(code, url, "")
        return self._save_fetched_file(handler.read())

    def _save_fetched_file(self, data):
        # Make sure that the root directory exists.
        try:
            os.makedirs(self.base_dir)
        except os.error:
            pass

        # Save the contents of the file.
        filename = self.filename
        logging.info("Saving file to: '{0}'".format(filename))

        # Atomically write to disk.
        # http://stackoverflow.com/questions/2333872/ \
        #        atomic-writing-to-file-with-python
        f = NamedTemporaryFile("wb", delete=False)
        f.write(data)
        f.flush()
        os.fsync(f.fileno())
        f.close()
        shutil.move(f.name, filename)

        return self

    def plot(self):
        if np is None:
            raise ImportError("numpy is required for plotting.")
        if pl is None:
            raise ImportError("matplotlib is required for plotting.")


class LightCurve(_datafile):
    """
    A reference to a light curve dataset on the MAST severs. This object
    handles local caching of the file in a strict directory structure.

    """

    product = "lightcurves"
    suffixes = ["llc", "slc"]
    filetype = ".fits"

    def plot(self):
        """
        Make a simple diagnostic plot of the light curve.

        :return fig:
            A :class:`matplotlib.Figure` object containing the plot.

        """
        super(LightCurve, self).plot()

        # Load the data.
        with self.open() as f:
            data = f[1].data
        time, sapflux, pdcflux, qual = (data["time"], data["sap_flux"],
                                        data["pdcsap_flux"],
                                        data["sap_quality"])

        # Set up the figure.
        fig, axes = pl.subplots(2, 1, figsize=(6, 6))
        fig.subplots_adjust(left=0.17, top=0.99, right=0.99,
                            wspace=0.0, hspace=0.0)

        # Plot the data.
        m = np.isfinite(time)
        xlim = [np.min(time[m]), np.max(time[m])]
        for i, (f, nm) in enumerate(zip([sapflux, pdcflux],
                                        ["sap flux", "pdc flux"])):
            ax = axes[i]
            m = np.isfinite(time) * np.isfinite(f)
            m1 = m * (qual == 0)
            m2 = m * (qual != 0)
            mu = np.median(f[m1])
            f1 = (f[m1] / mu - 1) * 1e6
            ax.plot(time[m1], f1, ".k", ms=3, alpha=0.6)
            ax.plot(time[m2], (f[m2] / np.median(f[m2]) - 1) * 1e6, ".r", ms=3,
                    alpha=0.6)
            ax.set_xlim(xlim)
            ax.set_ylim(np.min(f1), np.max(f1))
            ax.annotate("relative " + nm + " [ppm]",
                        xy=(1, 1), xycoords="axes fraction",
                        xytext=(-3, -3), textcoords="offset points",
                        horizontalalignment="right", verticalalignment="top")

        axes[0].set_xticklabels([])
        axes[1].set_xlabel("time [KBJD]")

        return fig


class TargetPixelFile(_datafile):
    """
    A reference to a target pixel dataset on the MAST severs. Like the
    :class:`LightCurve` object, this object handles local caching of the file
    in a strict directory structure.

    """

    product = "target_pixel_files"
    suffixes = ["lpd-targ", "spd-targ"]
    filetype = ".fits.gz"

    def plot(self, normed=False):
        """
        Make a simple diagnostic plot of the target pixel light curves. The
        pixels in the optimal aperture are plotted in black and those outside
        are plotted in red.

        :return fig:
            A :class:`matplotlib.Figure` object containing the grid of axes.

        """
        super(TargetPixelFile, self).plot()

        # Load the data.
        with self.open() as f:
            data = f[1].data
            aperture = f[2].data
        time, flux = data["time"], data["flux"]
        nx, ny = flux[0].shape

        # Set up the figures.
        factor = 3.0
        dx = factor * nx
        dy = factor * ny
        fig, axes = pl.subplots(nx, ny, figsize=(dy, dx))
        fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0,
                            wspace=0.0, hspace=0.0)

        # Loop over the pixels.
        m = np.isfinite(time)
        xlim = [np.min(time[m]), np.max(time[m])]
        ylim = [0., 0.]
        for xi, yi in product(range(nx), range(ny)):
            ax = axes[xi, yi]
            f = flux[:, xi, yi]

            m = np.isfinite(f)
            if normed:
                if len(f[m]) > 0:
                    ax.plot(time[m], (f[m]/np.median(f[m]))-1.0, ".",
                            color="rk"[int(aperture[xi, yi] == 3)], ms=2)
                    minval = np.where(
                        np.min((f[m]/np.median(f[m]))-1.0) < ylim[0],
                        np.min((f[m]/np.median(f[m]))-1.0),
                        ylim[0])
                    maxval = np.where(
                        np.max((f[m]/np.median(f[m]))-1.0) > ylim[1],
                        np.max((f[m]/np.median(f[m]))-1.0),
                        ylim[1])
                    minmax = [minval, maxval]
            else:
                ax.plot(time[m], f[m], ".",
                        color="rk"[int(aperture[xi, yi] == 3)], ms=2)

            ax.set_xlim(xlim)
            ax.set_xticklabels([])
            ax.set_yticklabels([])

        # If normed, put everything on the same scale
        if normed:
            for ax in axes.flatten():
                try:
                    ax.set_ylim(minmax)
                except TypeError:
                    continue

        return fig


class K2TargetPixelFile(TargetPixelFile):
    """
    A reference to a K2 target pixel dataset on the MAST severs. Like the
    :class:`LightCurve` object, this object handles local caching of the file
    in a strict directory structure.

    """

    _id = "\"{sci_data_set_name}_{ktc_target_type}\""
    kepid_template = "{ktc_k2_id:09d}"
    product = "target_pixel_files"
    suffixes = ["lpd-targ", "spd-targ"]
    filetype = ".fits.gz"

    @property
    def base_dir(self):
        """
        The local base directory for these types of products.

        """
        return os.path.join(self.api.data_root, "data", "k2", self.product,
                            self.kepid)

    @property
    def filename(self):
        """
        The local filename of the data file. This file is only guaranteed to
        exist after ``fetch()`` has been called.

        """
        return os.path.join(self.base_dir, self._filename)

    @property
    def url(self):
        """
        The remote URL for the data file on the MAST servers.

        """
        base_url = "http://archive.stsci.edu/pub/k2/"
        if self.ktc_k2_id < 201000000:
            base_url += "{0}/c%d/200000000/{1:05d}/{2}" % self.sci_campaign
        else:
            base_url += ("{{0}}/c%d/{0}/{{1:05d}}/{{2}}" % self.sci_campaign) \
                .format(int(int(self.ktc_k2_id * 1e-5) * 1e5))

        return base_url.format(self.product,
                               int(int(int(self.kepid[-5:][-5:])*1e-3)*1e3),
                               self._filename)

class k2sff(object): 
    pass

def K2SFF(EPIC, version = 1, clobber = False, sci_campaign = None):
    '''
    
    '''
    
    # Local directory 
    base_dir = os.path.join(KPLR_ROOT, "data", "k2sff", str(EPIC))
    
    # Check for local copies
    file_exists = False
    if sci_campaign is None:
        for c in range(30):
            filename = "hlsp_k2sff_k2_lightcurve_%09d-c%02d_kepler_v%d_llc.fits" % \
                       (EPIC, c, version)
            if os.path.exists(os.path.join(base_dir, filename)):
                sci_campaign = c
                file_exists = True
                break
    else:
        filename = "hlsp_k2sff_k2_lightcurve_%09d-c%02d_kepler_v%d_llc.fits" % \
                   (EPIC, sci_campaign, version)
        if os.path.exists(os.path.join(base_dir, filename)):
            file_exists = True
    
    # Download the data
    if clobber or not file_exists:
      
      # Get the campaign number
      if sci_campaign is None:
        client = API()
        star = client.k2_star(EPIC)
        tpf = star.get_target_pixel_files(clobber = clobber)
        sci_campaign = tpf[0].sci_campaign
        filename = "hlsp_k2sff_k2_lightcurve_%09d-c%02d_kepler_v%d_llc.fits" % \
                   (EPIC, sci_campaign, version)
      
      # Get the url
      first_four = int(str(EPIC)[:4])
      last_five = int(str(EPIC)[-5:])
      url = "http://archive.stsci.edu/missions/hlsp/k2sff/"
      url += "c%02d/%04d00000/%05d/" % (sci_campaign, first_four, last_five)
      url += filename

      # Query the server
      r = urllib.request.Request(url)
      handler = urllib.request.urlopen(r)
      code = handler.getcode()
      if int(code) != 200:
          raise APIError(code, url, "")
      data = handler.read()
      
      # Make sure that the root directory exists.
      try:
          os.makedirs(base_dir)
      except os.error:
          pass

      # Atomically write to disk.
      f = NamedTemporaryFile("wb", delete=False)
      f.write(data)
      f.flush()
      os.fsync(f.fileno())
      f.close()
      shutil.move(f.name, os.path.join(base_dir, filename))
  
    # Now open the fits file 
    res = k2sff() 
    res._file = os.path.join(base_dir, filename)
    with pyfits.open(os.path.join(base_dir, filename)) as f:  
        res.time = f[1].data['T']
        res.fcor = f[1].data['FCOR']
        res.apertures = np.array(np.vstack([f[22].data, f[23].data]), dtype = int)
    
    return res

class k2varcat(object): 
    pass

def K2VARCAT(EPIC, version = 2, clobber = False, sci_campaign = None):
    '''
    
    '''
    
    # Local directory 
    base_dir = os.path.join(KPLR_ROOT, "data", "k2varcat", str(EPIC))
    
    # Check for local copies
    file_exists = False
    if sci_campaign is None:
        for c in range(30):
            filename = "hlsp_k2varcat_k2_lightcurve_%09d-c%02d_kepler_v%d_llc.fits" % \
                       (EPIC, c, version)
            if os.path.exists(os.path.join(base_dir, filename)):
                sci_campaign = c
                file_exists = True
                break
    else:
        filename = "hlsp_k2varcat_k2_lightcurve_%09d-c%02d_kepler_v%d_llc.fits" % \
                   (EPIC, sci_campaign, version)
        if os.path.exists(os.path.join(base_dir, filename)):
            file_exists = True
    
    # Download the data
    if clobber or not file_exists:
      
      # Get the campaign number
      if sci_campaign is None:
        client = API()
        star = client.k2_star(EPIC)
        tpf = star.get_target_pixel_files(clobber = clobber)
        sci_campaign = tpf[0].sci_campaign
        filename = "hlsp_k2varcat_k2_lightcurve_%09d-c%02d_kepler_v%d_llc.fits" % \
                   (EPIC, sci_campaign, version)
      
      # Get the url
      first_four = int(str(EPIC)[:4])
      next_two = int(str(EPIC)[-5:-3] + '000')
      url = "https://archive.stsci.edu/missions/hlsp/k2varcat/"
      url += "c%02d/%04d00000/%05d/" % (sci_campaign, first_four, next_two)
      url += filename

      # Query the server
      r = urllib.request.Request(url)
      handler = urllib.request.urlopen(r)
      code = handler.getcode()
      if int(code) != 200:
          raise APIError(code, url, "")
      data = handler.read()
      
      # Make sure that the root directory exists.
      try:
          os.makedirs(base_dir)
      except os.error:
          pass

      # Atomically write to disk.
      f = NamedTemporaryFile("wb", delete=False)
      f.write(data)
      f.flush()
      os.fsync(f.fileno())
      f.close()
      shutil.move(f.name, os.path.join(base_dir, filename))
  
    # Now open the fits file 
    res = k2varcat() 
    res._file = os.path.join(base_dir, filename)
    with pyfits.open(os.path.join(base_dir, filename)) as f:  
        res.time = f[1].data['TIME']
        res.flux = f[1].data['DETFLUX']
    
    return res

class k2sc(object): 
    pass

def K2SC(EPIC, version = 1, clobber = False, sci_campaign = None):
    '''
    
    '''
    
    # Local directory 
    base_dir = os.path.join(KPLR_ROOT, "data", "k2sc", str(EPIC))
    
    # Check for local copies
    file_exists = False
    if sci_campaign is None:
        for c in range(30):
            filename = "hlsp_k2sc_k2_llc_%09d-c%02d_kepler_v%d_lc.fits" % \
                       (EPIC, c, version)
            if os.path.exists(os.path.join(base_dir, filename)):
                sci_campaign = c
                file_exists = True
                break
    else:
        filename = "hlsp_k2sc_k2_llc_%09d-c%02d_kepler_v%d_lc.fits" % \
                   (EPIC, sci_campaign, version)
        if os.path.exists(os.path.join(base_dir, filename)):
            file_exists = True
    
    # Download the data
    if clobber or not file_exists:
      
      # Get the campaign number
      if sci_campaign is None:
        client = API()
        star = client.k2_star(EPIC)
        tpf = star.get_target_pixel_files(clobber = clobber)
        sci_campaign = tpf[0].sci_campaign
        filename = "hlsp_k2sc_k2_llc_%09d-c%02d_kepler_v%d_lc.fits" % \
                   (EPIC, sci_campaign, version)
      
      # Get the url
      first_four = int(str(EPIC)[:4])
      next_two = int(str(EPIC)[-5:-3] + '000')
      url = "https://archive.stsci.edu/missions/hlsp/k2sc/"
      url += "c%02d/%04d00000/" % (sci_campaign, first_four)
      url += filename

      # Query the server
      r = urllib.request.Request(url)
      handler = urllib.request.urlopen(r)
      code = handler.getcode()
      if int(code) != 200:
          raise APIError(code, url, "")
      data = handler.read()
      
      # Make sure that the root directory exists.
      try:
          os.makedirs(base_dir)
      except os.error:
          pass

      # Atomically write to disk.
      f = NamedTemporaryFile("wb", delete=False)
      f.write(data)
      f.flush()
      os.fsync(f.fileno())
      f.close()
      shutil.move(f.name, os.path.join(base_dir, filename))
  
    # Now open the fits file 
    res = k2sc() 
    res._file = os.path.join(base_dir, filename)
    with pyfits.open(os.path.join(base_dir, filename)) as f:  
        res.time = f[1].data['TIME']
        
        # Get the data minus the instrumental trends only
        res.pdcflux = f[1].data['FLUX'] + f[1].data['TREND_T'] - np.median(f[1].data['TREND_T'])
        res.sapflux = f[2].data['FLUX'] + f[2].data['TREND_T'] - np.median(f[2].data['TREND_T'])
        
    return res