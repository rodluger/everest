Using Everest
-------------

.. role:: python(code)
   :language: python

There are two ways of interacting with the :py:obj:`everest` catalog: via the command line and 
through the Python interface. For quick visualization, check out the :doc:`everest <everest>` and
:doc:`estats <estats>` command line tools.
For customized access to de-trended light curves in Python, keep reading!

A Simple Example
================

Once you've installed :py:obj:`everest`, you can easily import it in Python:

.. code-block :: python
   
   import everest

Say we're interested in `EPIC 201367065`, a :py:obj:`K2` transiting exoplanet host star.
Let's instantiate the :py:class:`Everest <everest.user.Everest>` class for this target:

.. code-block :: python
   
   star = everest.Everest(201367065)

You will see the following printed to the screen:

.. code-block :: bash
  
   INFO  [everest.user.DownloadFile()]: Downloading the file...
   INFO  [everest.user.load_fits()]: Loading FITS file for 201367065.

Everest automatically downloaded the light curve and created an object containing all of
the de-trending information. For more information on the various methods and attributes 
of :py:obj:`star`, check out the 
:py:class:`everest.Detrender <everest.detrender.Detrender>` documentation (from which 
:py:class:`Everest <everest.user.Everest>` inherits a bunch of parameters) and the
:py:class:`Everest <everest.user.Everest>` docstring.

To bring up the DVS (:doc:`data validation summary <dvsfigs>`) for the target, execute

.. code-block :: python
   
   star.dvs()

You can also plot it interactively:

.. code-block :: python
  
   star.plot()

.. figure:: everest_plot.jpeg
   :width: 600px
   :align: center
   :figclass: align-center

The raw light curve is shown at the top and the de-trended light curve at the bottom.
The 6 hr CDPP (a photometric precision metric) is shown at the top of each plot in
red. Since this light curve was de-trended with a break point, which divides it into
two separate segments, the CDPP is shown for each one. At the top, below the title,
we indicate the CDPP for the entire light curve (raw → de-trended). Outliers are 
indicated in red, and arrows indicate points that are beyond the limits of the plot
(zoom out to see them). You can read more about these plots :doc:`here <dvsfigs>`.

Finally, if you want to manipulate the light curve yourself, the timeseries is stored
in :python:`star.time` and :python:`star.flux` (PLD-de-trended flux) or :python:`star.fcor` (de-trended
flux with CBV correction). The indices of all outliers are stored in :python:`star.mask`.

Masking Transits
================

CBV Corrections
===============

Tuning the Model
================

Pipeline Comparison
===================

Folded Transits
===============

If the time of first transit and period of an exoplanet/EB are known, plotting the
folded transit/eclipse is easy:

.. code-block :: python
  
   star.plot_folded(1980.42, 10.054)

.. figure:: everest_folded.jpeg
   :width: 400px
   :align: center
   :figclass: align-center