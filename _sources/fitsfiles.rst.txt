FITS Files
==========

The FITS files available in the :py:mod:`everest` catalog contain the de-trended
light curves, as well as information used by :py:mod:`everest` to do
customized post-processing. Each FITS file is composed of six extensions:

.. contents::
   :local:

.. note:: It is **highly recommended** that users access the catalog through the \
          :py:mod:`everest` :doc:`user interface <ui>`, as this allows for post-processing of \
          the light curves with custom masks.


``[0]`` Primary HDU
~~~~~~~~~~~~~~~~~~~

A header unit containing general information about the data collection for this target;
this info was copied directly from the original *K2* target pixel files. The ``VERSION``
and ``DATE`` tags were added to indicate the :py:mod:`everest` pipeline version and
creation date, respectively.

``[1]`` Lightcurve HDU
~~~~~~~~~~~~~~~~~~~~~~

The main data unit containing the de-trended light curve. Again, the header keywords
were taken from the equivalent extension in the original *K2* target pixel files. The
following :py:mod:`everest` keywords were added:
  
================  =================================================================================
 **Keyword**       **Description**
----------------  ---------------------------------------------------------------------------------
MISSION           The photometry mission name
VERSION           The EVEREST pipeline version
DATE              FITS file creation date (YYYY-MM-DD)
MODEL             The name of the :py:obj:`everest` model used to de-trend
APNAME            The name of the aperture used
BPAD              Light curve segment overlap in cadences for mending at breakpoints
BRKPT\ **XX**     Index/indices of light curve breakpoint(s)
CBVNUM            Number of CBV signals to regress on
CBVNITER          Number of CBV SysRem iterations
CBVWIN            Window size for smoothing CBVs
CBVORD            Filter order when smoothing CBVs
CDIVS             Number of cross-validation subdivisions
CDPP              Average de-trended CDPP
CDPPR             Raw light curve CDPP
CDPPV             Estimated average validation CDPP
CDPPG             Average GP-de-trended CDPP
CDPP\ **XX**      De-trended CDPP in light curve segment **XX**
CDPPR\ **XX**     Raw CDPP in light curve segment **XX**
CDPPV\ **XX**     Estimated validation CDPP in light curve segment **XX**
CVMIN             Cross-validation objective function
GITER             Number of GP optimiziation iterations
GPFACTOR          GP amplitude initialization factor
GPWHITE           GP white noise amplitude (e-/s)
GPRED             GP red noise amplitude (e-/s)
GPTAU             GP red noise timescale (days)
LAM\ **XXYY**     Cross-validation parameter for segment **XX** and PLD order **YY**
RECL\ **XXYY**    Cross-validation parameter for segment **XX** and PLD order **YY** (:py:class:`iPLD` only)
LEPS              Cross-validation minimum-finding CDPP tolerance
MAXPIX            Maximum size of TPF aperture in pixels
NRBY\ **XX**\ ID  Nearby source **XX** ID number
NRBY\ **XX**\ X   Nearby source **XX** X position (pixels)
NRBY\ **XX**\ Y   Nearby source **XX** Y position (pixels)
NRBY\ **XX**\ M   Nearby source **XX** magnitude
NRBY\ **XX**\ X0  Nearby source reference X (pixels)
NRBY\ **XX**\ Y0  Nearby source reference Y (pixels)
NEIGH\ **XX**     Neighboring star ID used to de-trend (:py:class:`nPLD` only)
OITER             Number of outlier clipping iterations
OPTGP             GP optimization performed?
OSIGMA            Outlier tolerance (standard deviations)
P\ **XX**\ T0     Masked planet transit time (days)
P\ **XX**\ PER    Masked planet period (days)
P\ **XX**\ DUR    Masked planet transit duration (days)
PLDORDER          PLD de-trending order
SATUR             Is target saturated?
SATTOL            Fractional saturation tolerance
================  =================================================================================

The ``data`` container of this extension contains the following arrays:

==============  =================================================================================
  **Keyword**     **Description**
--------------  ---------------------------------------------------------------------------------
TIME            The original timestamp. For :py:obj:`K2`, this is :py:obj:`(BJD - 2454833)`
CADN            The original cadence number
FLUX            The :py:mod:`everest` de-trended flux, same units as original SAP flux (e-/s)
FCOR            The CBV-corrected de-trended flux (e-/s)
FRAW            The original (raw) SAP flux
FRAW_ERR        The observing errors on the raw flux
QUALITY         An :py:obj:`int64` array of quality flags for each cadence (see note below)
BKG             If present, the background flux subtracted from each cadence
CBV\ **XX**     The co-trending basis vectors used to produce the corrected flux :py:obj:`FCOR`
==============  =================================================================================

.. note:: The :py:obj:`QUALITY` array uses the same bit flags as `K2`, with the addition of \
          four :py:mod:`everest` flags that indicate a data point was masked when computing the model:
            
            ====== =================================================
            **23** Data point is flagged in the raw `K2` TPF
            **24** Data point is a :py:obj:`NaN`
            **25** Data point was determined to be an outlier
            **26** *Not used*
            **27** Data point is during a transit/eclipse
            ====== =================================================
          
          As an example, to get the indices of all points that were flagged with bit **23** in
          Python, run
          
          .. code-block:: python
             
             inds = np.where(QUALITY & 2 ** (23 - 1))[0]
          
          Note, however, that the :py:obj:`everest` :doc:`user interface <ui>` provides easy access to
          these masks via the :py:attr:`Everest.badmask`, :py:attr:`Everest.nanmask`,
          :py:attr:`Everest.outmask`, and :py:attr:`Everest.transitmask` attributes.
          
``[2]`` Pixels HDU
~~~~~~~~~~~~~~~~~~

An extension containing the pixel-level light curve.
The ``data`` container stores two arrays:

==============  =================================================================================
 **Keyword**      **Description**
--------------  ---------------------------------------------------------------------------------
FPIX            The flux in each of the pixels in the aperture
X1N             The first order PLD vectors for the neighbors (:py:class:`nPLD` only)
==============  =================================================================================

``[3]`` Aperture HDU
~~~~~~~~~~~~~~~~~~~~

An extension containing the aperture mask used for PLD de-trending. Ones correspond to pixels
that were included in the de-trending; zeros correspond to pixels that were ignored.


``[4]`` Images HDU
~~~~~~~~~~~~~~~~~~

Stores images of the full target postage stamp at three points in the light curve, for
plotting purposes only.

==============  =================================================================================
 **Keyword**      **Description**
--------------  ---------------------------------------------------------------------------------
STAMP1          The postage stamp at the first cadence
STAMP2          The postage stamp at the midpoint
STAMP3          The postage stamp at the last cadence
==============  =================================================================================


``[5]`` HiRes HDU
~~~~~~~~~~~~~~~~~

An image HDU containing a higher resolution image of the target. For :py:obj:`K2`, this is
obtained from the Palomar Observatory Sky Survey.

.. raw:: html

  <script>
    (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
    (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
    m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
    })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

    ga('create', 'UA-47070068-3', 'auto');
    ga('send', 'pageview');

  </script>
