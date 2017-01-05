:py:mod:`everest` - Command line tools
--------------------------------------

The :py:mod:`everest` command accepts several options, which we list below.

======================  =============================================================================================
:py:obj:`ID`            The target :py:obj:`ID` number (`required`)
:py:obj:`-1`            Plot the :py:obj:`everest 1.0` de-trended light curve (:py:obj:`K2` `only`)
:py:obj:`-2`            Plot the :py:obj:`everest 2.0` de-trended light curve (`default`)
:py:obj:`-a`            Plot the target aperture/postage stamp
:py:obj:`-c` `cadence`  Cadence type (:py:obj:`lc` | :py:obj:`sc`). Default :py:obj:`lc`
:py:obj:`-d`            Show the :py:obj:`everest 2.0` data validation summary
:py:obj:`-f`            Plot the :py:obj:`k2sff` light curve for the target (:py:obj:`K2` `only`)
:py:obj:`-h`            Show the help message and exit
:py:obj:`-m` `mission`  The mission name (:py:obj:`k2` | :py:obj:`kepler` | :py:obj:`tess`). Default :py:obj:`k2`
:py:obj:`-n`            Do not apply the CBV correction
:py:obj:`-r`            Plot the raw light curve
:py:obj:`-s`            Plot the :py:obj:`k2sc` light curve for the target (:py:obj:`K2` `only`)
======================  =============================================================================================