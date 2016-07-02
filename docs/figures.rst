Figures
=======

Below are descriptions of the figures produced by the :py:mod:`everest` pipeline.
These are available as JPEGs in the same directory as the FITS files on MAST.

.. contents::
   :local:

Aperture
~~~~~~~~

.. figure:: 205071984_aper.jpeg
 :width: 600px
 :align: center
 :height: 100px
 :figclass: align-center

The thumbnails at the left show each of the 20 apertures derived in the
`K2SFF <https://archive.stsci.edu/prepds/k2sff/>`_ pipeline. The first ten
are increasingly larger circular apertures centered on the target; the 
last ten are derived from PSF fitting. For all :py:mod:`everest 0.1` 
runs, aperture 15 was used; this aperture is enlarged at the right. The
color scale is logarithmic, and the image shown is the flux integrated
*over the entire campaign*, so that the total motion of the sources is
visible. EPIC sources are overlaid based on positions obtained from MAST;
the green circle is the target, and red circles are neighboring sources.
Their `Kepler` magnitudes are also indicated.

Contamination
~~~~~~~~~~~~~

.. figure:: 205071984_contamination.jpeg
 :width: 600px
 :align: center
 :height: 100px
 :figclass: align-center

The left four panels show the contamination statistics for the first
timestamp: the actual stellar image (top left), the optimized error function
PSF model (top center), the model binned to the pixels (bottom left), and the
difference between the data and the model (bottom center). As in the previous
figure, the thick black line is the aperture, and the circles are nearby EPIC targets.

The top right panel shows the absolute value of the difference between the data
and the model divided by the flux in the brightest pixel within the aperture. The
largest value in this image is the contamination metric for this timestamp:
it is a measure of how well the aperture can be modeled with a single symmetrical
PSF. The bottom right panel shows the contamination metric evaluated 50 times throughout 
the campaign. The median value is indicated with a dashed line: this is the value
reported in the ``CONTAM`` field of the FITS file. Typically, :py:mod:`everest` 
performs well for targets with contamination less than about 0.3.

Note that sources with asymmetrical PSFs (due to distortions in the `Kepler` PRF)
or very saturated stars with bleeding trails will generally have large contamination
metrics, even if there are no other sources in the aperture.

Outliers
~~~~~~~~

.. figure:: 205071984_outliers.jpeg
 :width: 600px
 :align: center
 :height: 100px
 :figclass: align-center

The top panel shows the background flux in the aperture. For campaigns 0, 1, and 2,
this was computed by taking the median flux in the region of the postage stamp
outside the aperture at each timestamp. For campaigns 3 onwards, the `Kepler` team
pre-subtracts the background in the target pixel files, so no background subtraction
was performed.

The bottom panel shows the background-corrected SAP flux in the aperture. Black points
were used to compute the PLD fit; red points were identified as outliers by iterative
sigma-clipping and did not inform the modeling process. Note, however, that these
points are still de-trended and are present in the FITS files.

Autocorrelation
~~~~~~~~~~~~~~~

.. figure:: 205071984_acor.jpeg
 :width: 600px
 :align: center
 :height: 100px
 :figclass: align-center

These plots show the GP fitting process. The top panel is the PLD-de-trended data after
a full iteration of the GP optimizer. The center panel is the Lomb-Scargle periodogram
of the de-trended data. Three main periods are highlighted; these are used as initial
guesses when fitting periodic kernels to the autocorrelation function (bottom panel).

In the bottom panel, the autocorrelation function is the black line, 
and the best-fit kernel is the blue line. Its functional form is indicated in the legend
at the top right.
The standard deviation of the autocorrelation function is indicated by the grey envelope,
and was arbitrarily selected to downweight points at large lags. The chi-squared of the
fit is indicated at the top left, along with the number of iterations performed.
Finally, the numbers at the bottom left indicate the amplitude of the white noise component
*W* in units of the data error bars and the amplitude of the red noise component *A* in
units of the standard deviation of the de-trended light curve.

Cross-validation
~~~~~~~~~~~~~~~~

.. figure:: 205071984_crossvalidation.jpeg
 :width: 600px
 :align: center
 :height: 100px
 :figclass: align-center

This plot shows how :py:mod:`everest` selected the optimal number of principal components
to use in the de-trending step. The CDPP in the training set as a function of the number
of components is shown as blue points; the blue curve is a GP fit to the points. The CDPP
in the validation set (which was masked during de-trending) is shown as red points, with the
corresponding GP fit indicated by the red curve. The minimum in the red curve is chosen
as the optimal number of components. Beyond this point, the validation CDPP increases because
of overfitting.

Note that some cross-validation plots may look very different from this one. In particular, if
the red curve falls sufficiently below the blue curve, or if the red curve monotonically
decreases past 250 components, this may be a sign that the target is saturated/contaminated.
Highly oscillating cross-validation curves may also indicate poor :py:mod:`everest` performance.

De-trended
~~~~~~~~~~

.. figure:: 205071984_detrended.jpeg
 :width: 600px
 :align: center
 :height: 100px
 :figclass: align-center

The top plot shows the raw background-subtracted SAP flux in red (outliers in black);
the bottom plot shows the optimum de-trended light curve  in blue. The legends at the
top left indicate the CDPP before (left) and after (right) applying a 2-day, quadratic 
Savitsky-Golay filter. The approximate photon limit is also indicated in the legend of
the top panel.


.. raw:: html

  <script>
    (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
    (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
    m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
    })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

    ga('create', 'UA-47070068-2', 'auto');
    ga('send', 'pageview');
  </script>