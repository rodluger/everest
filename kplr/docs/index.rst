.. module:: kplr

A Python interface to the Kepler data
=====================================

If you're here, then you probably already know about the `Kepler mission
<http://kepler.nasa.gov/>`_. You probably also know that it can be a bit of a
pain to get access to this public dataset. As I understand things, the
canonical source for a catalog of planet candidates—or more precisely Kepler
Objects of Interest (KOIs)—is the `NASA Exoplanet Archive
<http://exoplanetarchive.ipac.caltech.edu/>`_ but the data is available
through `MAST <http://archive.stsci.edu/>`_. There are programmatic interfaces
(APIs) available for both of these services but it can still be frustrating to
interact with them in an automated way. That's why I made **kplr**.

**kplr** provides a lightweight Pythonic interface to the catalogs and data.
Below, I'll describe the features provided by kplr but to get things
started, let's see an example of how you would go about finding the published
parameters of a KOI and download the light curve data.

.. code-block:: python

    import kplr
    client = kplr.API()

    # Find a KOI.
    koi = client.koi(952.01)
    print(koi.koi_period)

    # This KOI has an associated star.
    star = koi.star
    print(star.kic_teff)

    # Download the lightcurves for this KOI.
    lightcurves = koi.get_light_curves()
    for lc in lightcurves:
        print(lc.filename)


Table of Contents
-----------------

* `Installation`_
* `API Interface`_
    * `Kepler Objects of Interest`_
    * `Confirmed Planets`_
    * `Kepler Input Catalog Targets`_
* `Data Access`_

For a detailed description of the Python bindings see the Python API
documentation:

.. toctree::
   :maxdepth: 2

   api


Installation
------------

You can install kplr using the standard Python packaging tool `pip
<http://www.pip-installer.org/>`_:

.. code-block:: bash

    pip install kplr

or (if you must) `easy_install <https://pypi.python.org/pypi/setuptools>`_:

.. code-block:: bash

    easy_install kplr

The development version can be installed using pip:

.. code-block:: bash

    pip install -e git+https://github.com/dfm/kplr#egg=kplr-dev

or by cloning the `GitHub repository <https://github.com/dfm/kplr>`_:

.. code-block:: bash

    git clone https://github.com/dfm/kplr.git
    cd kplr
    python setup.py install


API Interface
-------------

The basic work flow for interacting with the APIs starts by initializing an
:class:`API` object:

.. code-block:: python

    import kplr
    client = kplr.API()

Then, this object provides methods for constructing various queries to find

* `Kepler objects of interest`_,
* `confirmed planets`_, and
* `Kepler input catalog targets`_.

Kepler Objects of Interest
^^^^^^^^^^^^^^^^^^^^^^^^^^

The kplr KOI search interfaces with `The Exoplanet Archive API
<http://exoplanetarchive.ipac.caltech.edu/docs/program_interfaces.html>`_ to
return the most up to date information possible. In particular, it searches
the `cumulative table
<http://exoplanetarchive.ipac.caltech.edu/docs/API_kepcandidate_columns.html>`_.
As shown in the sample code at the top of this page, it is very easy to
retrieve the listing for a single :class:`KOI`:

.. code-block:: python

    koi = client.koi(952.01)

Note the ``.01`` in the KOI ID. This is required because a KOI is specified by
the full number (not just ``952`` which will fail).
The object will have an attribute for each column listed in the `Exoplanet
Archive documentation
<http://exoplanetarchive.ipac.caltech.edu/docs/API_kepcandidate_columns.html>`_.
For example, the period and error bars (positive and negative respectively)
are given by

.. code-block:: python

    print(koi.koi_period, koi.koi_period_err1, koi.koi_period_err2)

For KOI 952.01, this result will print ``5.901269, 1.7e-05, -1.7e-05``.

Finding a set of KOIs that satisfy search criteria is a little more
complicated because you must provide syntax that is understood by the
Exoplanet Archive. For example, to find all the KOIs with period longer than
200 days, you would run

.. code-block:: python

    kois = client.kois(where="koi_period>200")

At the time of writing, this should return 224 :class:`KOI` objects. If you
then wanted to sort by period, you could include the ``sort`` keyword
argument:

.. code-block:: python

    kois = client.kois(where="koi_period>200", sort="koi_period")

or, equivalently,

.. code-block:: python

    kois = client.kois(where="koi_period>200", sort=("koi_period", 1))

You can specify the sort order to be descending by using

.. code-block:: python

    kois = client.kois(where="koi_period>200", sort=("koi_period", -1))


Confirmed Planets
^^^^^^^^^^^^^^^^^

The confirmed planet interface queries the `confirmed planets
<http://archive.stsci.edu/search_fields.php?mission=kepler_cp>`_ table using
the `MAST API <http://archive.stsci.edu/vo/mast_services.html>`_. To find a
specific planet using this interface, you can use the :func:`API.planet`
function

.. code-block:: python

    planet = client.planet("32b")

or equivalently

.. code-block:: python

    planet = client.planet("Kepler-32b")

This object has attributes for each column given in the `table in the MAST
documentation
<http://archive.stsci.edu/search_fields.php?mission=kepler_cp>`_. For example,
the corresponding KOI name for this planet is given by

.. code-block:: python

    print(planet.kepoi_name)

In this case, you should see ``952.01``.

The query syntax on MAST is a little different than on the Exoplanet Archive.
For example, to find planets with estimated radii less than 2 Earth radii,
you would run

.. code-block:: python

    planets = client.planets(koi_prad="<2")

where ``koi_prad`` is the name of a column in the `MAST documentation table
<http://archive.stsci.edu/search_fields.php?mission=kepler_cp>`_.

The syntax for sorting the results is the same as described above for the
KOIs. To sort the above search by period, you would run

.. code-block:: python

    planets = client.planets(koi_prad="<2", sort="koi_period")


Kepler Input Catalog Targets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Access to the Kepler Input Catalog (KIC) is also provided by the `MAST API's
kic10 table <http://archive.stsci.edu/search_fields.php?mission=kic10>`_. It
can, therefore, be queried using syntax similar to the confirmed planet table.
For example, a particular star can be found as follows

.. code-block:: python

    star = client.star(9787239)

Similarly, a query can be run on the table using the following syntax:

.. code-block:: python

    stars = client.stars(kic_teff="5700..5800")

To select a set of stars in a 2MASS color range with (non-NULL) estimated
temperatures, you would run something like:

.. code-block:: python

    stars = client.stars(kic_jkcolor="0.3..0.4", kic_teff="!\\null")

*Note: by default, the* :func:`API.stars` *endpoint is limited to 100 results
because it's very easy to time out the MAST server if you're not careful. To
change this behavior, you can specify the* ``max_records`` *keyword argument:*

.. code-block:: python

    stars = client.stars(kic_jkcolor="0.3..0.4", kic_teff="!\\null", max_records=500)


Data Access
-----------

*Note: to interact with the Kepler data, you will need to be able to read the
FITS files. kplr automatically supports loading the data using* `pyfits
<http://pythonhosted.org/pyfits/>`_ *so it's probably easiest to make sure
that you have that installed before trying the examples in this section.*

The MAST servers are the main source for the Kepler data products. kplr
supports two types of data: light curves and target pixel files. These
products are described in detail in the `Kepler Archive Manual
<http://archive.stsci.edu/kepler/manuals/archive_manual.pdf>`_ but in summary:

* the target pixel files contain the lightly-processed CCD readouts from small
  fields around the telemetered Kepler targets, and
* the light curve files contain the results of the aperture photometric
  pipeline applied to the pixel files and various housekeeping columns.

All of the objects described above (:class:`KOI`, :class:`Planet` and
:class:`Star`) have a :func:`get_light_curves` method and a
:func:`get_target_pixel_files` method. These methods return (possibly empty)
lists of :class:`LightCurve` and :class:`TargetPixelFile` objects,
respectively. Both of the above methods take three keyword arguments:
``short_cadence``, ``fetch`` and ``clobber``. ``short_cadence`` defaults to
``True`` and it decides whether or not the "short cadence" data should be
included. If it is ``False``, only the "long cadence" data are returned. If
``fetch`` is ``True``, the data are automatically downloaded from the MAST
server if they don't already exist locally. Otherwise, if ``fetch`` is
``False`` (default) the data aren't downloaded until the first time they are
opened. Finally, ``clobber`` sets the behavior when a local copy of the file
already exists. If a data file has been corrupted, it can be useful to set
``clobber=True`` to make sure that the bad file is overwritten.

Below is an example for the best practice for loading a set of light curves
for a particular object:

.. code-block:: python

    # Find the target KOI.
    koi = client.koi(952.01)

    # Get a list of light curve datasets.
    lcs = koi.get_light_curves(short_cadence=False)

    # Loop over the datasets and read in the data.
    time, flux, ferr, quality = [], [], [], []
    for lc in lcs:
        with lc.open() as f:
            # The lightcurve data are in the first FITS HDU.
            hdu_data = f[1].data
            time.append(hdu_data["time"])
            flux.append(hdu_data["sap_flux"])
            ferr.append(hdu_data["sap_flux_err"])
            quality.append(hdu_data["sap_quality"])

