FITS Files
==========

The FITS files available in the :py:mod:`everest` catalog contain the de-trended
light curves, as well as information used by :py:mod:`everest.usertools` to do
customized post-processing. Each FITS file is composed of 5 extensions:

.. contents::
   :local:

.. note:: It is **highly recommended** that users access the catalog through the \
          :py:mod:`everest.usertools` module, as this allows for post-processing of \
          the light curves with custom masks. For more information, check out the \
          `Running Everest <running_everest.html>`_ section.

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
   Keyword        Description 
--------------  ---------------------------------------------------------------------------------
VERSION         The EVEREST pipeline version
DATE            FITS file creation date (YYYY-MM-DD)
PLDORDER        Order of PLD used (3 for all *K2* targets)
NPC             Number of principal components used in de-trending
NSMAT           Number of sub-matrices (chunks) used. 2 for campaign 1; 1 for all other campaigns
APNUM           Number of the *K2SFF* aperture used (15 for all *K2* targets)
CDPP6RAW        6-hr CDPP estimate for the raw SAP flux (ppm)
CDPP6           6-hr CDPP estimate for the de-trended flux (ppm)
SATFLAG         Saturation flag (0-5). Beware of stars with ``SATFLAG > 2`` 
CRWDFLAG        Saturation flag (0-5). Beware of stars with ``CRWDFLAG > 2`` 
CONTAM          Median contamination metric used to compute ``CRWDFLAG``
ACRCHSQ         Autocorrelation fit chi squared. Beware of stars with ``ACRCHSQ > 30``
GPITER          Number of GP fitting iterations (default is 2)
GITHASH         :py:mod:`everest` git repository commit hash when file was generated
BRKPT0          Time at which light curve was split
==============  =================================================================================

The ``data`` container of this extension contains the following arrays:

==============  =================================================================================
   Keyword        Description 
--------------  ---------------------------------------------------------------------------------
TIME            The original `K2` timestamp (BJD - 2454833)
FLUX            The :py:mod:`everest` de-trended flux, same units as original SAP flux (e-/s)
OUTLIER         An integer array of outlier masks. 1 = outlier (not used in fit); 0 = used in fit
BKG_FLUX        The background flux subtracted from the SAP data
RAW_FLUX        The original SAP flux from the `K2` TPF
RAW_FERR        The error bars on the original SAP flux
==============  =================================================================================

``[2]`` PLD weights HDU
~~~~~~~~~~~~~~~~~~~~~~~

An extension containing the values of the PLD weights computed by regression. Recall that the PLD
model is given by

  :math:`m = \mathbf{X\cdot w}`

where :math:`w` is the vector of weights and :math:`X` is the design matrix (see below).

The ``data`` container stores a single array:

==============  =================================================================================
   Keyword        Description 
--------------  ---------------------------------------------------------------------------------
C               The PLD coefficients (weights) array
==============  =================================================================================

``[3]`` Design Matrix HDU
~~~~~~~~~~~~~~~~~~~~~~~~~

An extension containing the design matrix :math:`X` and some GP information. The ``header``
contains the following:

+------------+---------------------------------------------------------------------------------+
| Keyword    | Description                                                                     |
+============+=================================================================================+
| KNUM       | | The number of the kernel used; this is the number of the element in the       |
|            | | :py:obj:`everest.kernels.KernelModels` list. See `kernels.py <kernels.html>`_.|
+------------+---------------------------------------------------------------------------------+
| KPARXX     | | The value of the XXth element of the kernel, obtained during the              |
|            | | optimization step.                                                            |
+------------+---------------------------------------------------------------------------------+

If the FITS file handle is ``hdulist``, the :py:mod:`george` kernel used during de-trending may
be re-constructed by typing

.. code-block:: python

  knum = hdulist[3].header['KNUM']
  kpars = [hdulist[3].header['KPAR%02d' % n] for n in range(10)]
  kpars = [k for k in kpars if k != '']
  kernel = everest.kernels.KernelModels[knum]
  kernel[:] = kpars
  george_kernel = kernel.george_kernel()  

The ``data`` container stores a single array:

==============  =================================================================================
   Keyword        Description 
--------------  ---------------------------------------------------------------------------------
X               The PLD design matrix
==============  =================================================================================

``[4]`` Aperture HDU
~~~~~~~~~~~~~~~~~~~~

An extension containing the aperture mask used for PLD de-trending. Ones correspond to pixels
that were included in the de-trending; zeros correspond to pixels that were ignored.