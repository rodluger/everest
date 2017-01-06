FITS Files
==========

The FITS files available in the :py:mod:`everest` catalog contain the de-trended
light curves, as well as information used by :py:mod:`everest` to do
customized post-processing. Each FITS file is composed of six extensions:

.. contents::
   :local:

.. note:: It is **highly recommended** that users access the catalog through the \
          :py:mod:`everest` interface, as this allows for post-processing of \
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
  
==============  =================================================================================
 **Keyword**      **Description**
--------------  ---------------------------------------------------------------------------------
MISSION         The photometry mission name
VERSION         The EVEREST pipeline version
DATE            FITS file creation date (YYYY-MM-DD)
MODEL           The name of the :py:obj:`everest` model used to de-trend
APNAME          The name of the aperture used
BPAD            Light curve segment overlap in cadences for mending at breakpoints
BRKPTXX         Index/indices of light curve breakpoint(s)
CBVNUM          Number of CBV signals to regress on
CBVNITER        Number of CBV SysRem iterations
CBVWIN          Window size for smoothing CBVs
CBVORD          Filter order when smoothing CBVs
CDIVS           Number of cross-validation subdivisions
CDPP            Average de-trended CDPP
CDPPR           Raw light curve CDPP
CDPPV           Estimated average validation CDPP
CDPPG           Average GP-de-trended CDPP
CDPPXX          De-trended CDPP in light curve segment :py:obj:`XX`
CDPPRXX         Raw CDPP in light curve segment :py:obj:`XX`
CDPPVXX         Estimated validation CDPP in light curve segment :py:obj:`XX`
CVMIN           Cross-validation objective function
GITER           Number of GP optimiziation iterations
GPFACTOR        GP amplitude initialization factor
GPWHITE         GP white noise amplitude (e-/s)
GPRED           GP red noise amplitude (e-/s)
GPTAU           GP red noise timescale (days)
LAMXXYY         Cross-validation parameter for segment XX and PLD order YY
RECLXXYY        Cross-validation parameter for segment XX and PLD order YY (:py:class:`iPLD` only)
LEPS            Cross-validation minimum-finding CDPP tolerance
MAXPIX          Maximum size of TPF aperture in pixels
NRBYXXID        Nearby source XX ID number
NRBYXXX         Nearby source XX X position (pixels)
NRBYXXY         Nearby source XX Y position (pixels)
NRBYXXM         Nearby source XX magnitude
NRBYXXX0        Nearby source reference X (pixels)
NRBYXXY0        Nearby source reference Y (pixels)
NEIGHXX         Neighboring star ID used to de-trend (:py:class:`nPLD` only)
OITER           Number of outlier clipping iterations
OPTGP           GP optimization performed?
OSIGMA          Outlier tolerance (standard deviations)
PXXT0           Masked planet transit time (days)
PXXPER          Masked planet period (days)
PXXDUR          Masked planet transit duration (days)
PLDORDER        PLD de-trending order
SATUR           Is target saturated?
SATTOL          Fractional saturation tolerance
==============  =================================================================================

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
QUALITY         An :py:obj:`int64` array of quality flags for each cadence
BKG             If present, the background flux subtracted from each cadence
==============  =================================================================================

    
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
