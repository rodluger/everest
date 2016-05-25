#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import requests
from bs4 import BeautifulSoup

url = "http://archive.stsci.edu/search_fields.php"

types = {"string": "unicode", "datetime": "unicode", "float": "float",
         "integer": "int", "numeric": "unicode", "ra": "float",
         "dec": "float", "real": "float", "double": "float",
         "ustring": "unicode"}


def gen_adaptor(mission):
    r = requests.get(url, params={"mission": mission})
    if r.status_code != requests.codes.ok:
        r.raise_for_status()

    # Parse the HTML.
    tree = BeautifulSoup(r.content)

    # Find the table that describes the parameters.
    table_body = tree.find_all("tbody")
    assert table_body is not None

    for row in table_body[0].find_all("tr"):
        short_name, long_name, desc, ex, t = row.find_all("td")
        print("\"{0}\": (\"{1}\", {2}),".format(long_name.text.strip(),
                                                short_name.text.strip(),
                                                types[t.text.strip()]))


if __name__ == "__main__":
    print("Confirmed planets")
    gen_adaptor("kepler_cp")

    print("\nStars")
    gen_adaptor("kic10")

    print("\nData Search")
    gen_adaptor("kepler")

    print("\nK2 EPIC Search")
    gen_adaptor("epic")

    print("\nK2 TPFs")
    gen_adaptor("k2")
