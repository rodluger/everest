# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["get_quad_coeffs"]

import os
import shutil
import sqlite3
import logging
from tempfile import NamedTemporaryFile

from six.moves import urllib

from .config import KPLR_ROOT

DB_FILENAME = "ldcoeffs.db"


def get_quad_coeffs(teff=5778, logg=None, feh=None, data_root=None,
                    clobber=False):
    """
    Get the quadratic coefficients for the standard Kepler limb-darkening
    profile.

    :param teff: (optional)
        The effective temperature in degrees K.

    :param logg: (optional)
        The log10 surface gravity in cm/s/s.

    :param feh: (optional)
        The metallicity [Fe/H].

    :param data_root: (optional)
        The local base directory where the grids will be downloaded to. This
        can also be set using the ``KPLR_ROOT`` environment variable. The
        default value is ``~/.kplr``.

    :param clobber: (optional)
        Should the database file be overwritten even if it already exists?
        (default: False)

    """
    assert teff is not None

    # Make sure that the database is saved locally.
    filename = download_database(data_root=data_root, clobber=clobber)

    # Construct the SQL query.
    q = """
    SELECT mu1,mu2 FROM claret11 WHERE
    teff=(SELECT teff FROM claret11 ORDER BY abs(teff-?) LIMIT 1)
    ORDER BY (logg-?) * (logg-?) + (feh-?) * (feh-?) LIMIT 1
    """
    pars = [teff, logg, logg, feh, feh]

    # Execute the command.
    with sqlite3.connect(filename) as conn:
        c = conn.cursor()
        rows = c.execute(q, pars)
        mu1, mu2 = rows.fetchone()

    return mu1, mu2


def download_database(data_root=None, clobber=False):
    """
    Download a SQLite database containing the limb darkening coefficients
    computed by `Claret & Bloemen (2011)
    <http://adsabs.harvard.edu/abs/2011A%26A...529A..75C>`_. The table is
    available online on `Vizier
    <http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+A/529/A75>`_.
    Using the ASCII data table, the SQLite database was generated with the
    following Python commands:

    .. code-block:: python

        import sqlite3
        import numpy as np

        with sqlite3.connect("ldcoeffs.db") as conn:
            c = conn.cursor()
            c.execute("CREATE TABLE IF NOT EXISTS claret11 "
                    "(teff REAL, logg REAL, feh REAL, veloc REAL, mu1 REAL, "
                    "mu2 REAL)")
            data = np.loadtxt("claret11.txt", skiprows=59, delimiter="|",
                            usecols=range(1, 7))
            c.executemany("INSERT INTO claret11 (logg,teff,feh,veloc,mu1,mu2) "
                        "VALUES (?,?,?,?,?,?)", data)

    """
    # Figure out the local filename for the database.
    if data_root is None:
        data_root = KPLR_ROOT
    filename = os.path.join(data_root, DB_FILENAME)

    if not clobber and os.path.exists(filename):
        return filename

    # Make sure that the target directory exists.
    try:
        os.makedirs(data_root)
    except os.error:
        pass

    # MAGIC: specify the URL for the remote file.
    url = "http://bbq.dfm.io/~dfm/ldcoeffs.db"

    # Fetch the database from the server.
    logging.info("Downloading file from: '{0}'".format(url))
    r = urllib.request.Request(url)
    handler = urllib.request.urlopen(r)
    code = handler.getcode()
    if int(code) != 200:
        raise RuntimeError("Couldn't download file from {0}. Returned: {1}"
                           .format(url, code))

    # Save the contents of the file.
    logging.info("Saving file to: '{0}'".format(filename))

    # Atomically write to disk.
    # http://stackoverflow.com/questions/2333872/ \
    #        atomic-writing-to-file-with-python
    f = NamedTemporaryFile("wb", delete=False)
    f.write(handler.read())
    f.flush()
    os.fsync(f.fileno())
    f.close()
    shutil.move(f.name, filename)

    return filename
