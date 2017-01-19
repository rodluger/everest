Known Issues
============

Below we list some issues known to crop up here and there in the :py:obj:`everest`
light curves. This list is not comprehensive, so if you think you've found an issue
with the de-trending, please let us know by `opening an issue <https://github.com/rodluger/everest/issues>`_
on the github page.

.. contents::
   :local:

Missing eclipses
~~~~~~~~~~~~~~~~
Some very deep transits and eclipses may be missing from the plots in the
:doc:`data validation summaries <dvsfigs>`. This is almost always a plotting
issue, and not an issue with the pipeline. We recommend plotting the light curve interactively
using the :doc:`user interface <user>` and zooming out -- the eclipses are likely
just below the plot limits!

Stars in crowded apertures
~~~~~~~~~~~~~~~~~~~~~~~~~~
Everest 2.0 light curves are far more robust to crowding than the previous ones,
but stars in crowded fields remain an issue for the algorithm. PLD implicitly
assumes that the target star is the only source in the aperture and may overfit
astrophysical signals if bright contaminants are present. It is important to
inspect the apertures and high resolution images provided in the
:doc:`data validation summaries <dvsfigs>` to ensure that crowding is not affecting
the light curves. In general, crowding is only a problem when the contaminant source
is of a similar magnitude (or brighter) than the target and is inside the optimal
aperture. Future versions of the pipeline will improve aperture selection and
adapt PLD to work for these stars.

RR Lyrae
~~~~~~~~
Everest 2.0 again improves the performance of PLD on extremely variable stars
relative to Everest 1.0. Nonetheless, about 40% of long cadence RR Lyrae seem to be 
overfitted by the algorithm. We recommend inspecting the raw and de-trended light
curves in the :doc:`data validation summaries <dvsfigs>` to check whether PLD
has dampened the astrophysical signal.

Aperture losses
~~~~~~~~~~~~~~~
Extremely saturated stars whose flux bleeds out of the target pixel file postage
stamp are not properly de-trended, since a large portion of the astrophysical signal
is missing from the light curve. We recommend inspection of the target aperture to
ensure all of the stellar flux is enclosed by the chosen aperture.

Faint short cadence targets
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Extremely faint stars observed in short cadence mode are not properly de-trended by
PLD. This is because they are dominated by photon noise (instead of instrumental
noise) and the algorithm has trouble mending the light curve segments at the breakpoints,
resulting in discontinuous jumps for short cadence targets fainter than about 16th magnitude.


.. raw:: html

  <script>
    (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
    (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
    m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
    })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

    ga('create', 'UA-47070068-3', 'auto');
    ga('send', 'pageview');

  </script>