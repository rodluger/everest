# -*- coding: utf-8 -*-

from __future__ import division, print_function

__all__ = ["download"]

import os
import shutil
import logging
from tempfile import NamedTemporaryFile

from six.moves import urllib

try:
    import pandas as pd
except ImportError:
    pd = None

from .config import KPLR_ROOT

_URL = "https://zenodo.org/record/13297/files/huber-kic-join-2014-12-16.tsv.gz"
_FILENAME = os.path.join(KPLR_ROOT, "huber-kic-join-2014-12-16.tsv.gz")


def download(clobber=False):
    if os.path.exists(_FILENAME) and not clobber:
        print("File '{0}' already exists".format(_FILENAME))
        return

    try:
        os.makedirs(os.path.dirname(_FILENAME))
    except os.error:
        pass

    # Fetch the remote file.
    logging.info("Downloading file from: '{0}'".format(_URL))
    r = urllib.request.Request(_URL)
    handler = urllib.request.urlopen(r)
    code = handler.getcode()
    assert int(code) == 200, "Download returned HTTP error {0}".format(code)

    # Atomically write to disk.
    # http://stackoverflow.com/questions/2333872/ \
    #        atomic-writing-to-file-with-python
    logging.info("Saving file to: '{0}'".format(_FILENAME))
    f = NamedTemporaryFile("wb", delete=False)
    f.write(handler.read())
    f.flush()
    os.fsync(f.fileno())
    f.close()
    shutil.move(f.name, _FILENAME)


def get_catalog():
    if pd is None:
        raise ImportError("pandas is required to read the Huber catalog")
    if not os.path.exists(_FILENAME):
        download()
    return pd.read_csv(_FILENAME, sep="\t", compression="gzip")
