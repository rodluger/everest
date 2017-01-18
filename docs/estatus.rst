:py:mod:`everest-status` - De-trending Status
---------------------------------------------

The :py:mod:`everest-status` command prints the status of the batch de-trending runs
for a given mission, model, and season. The command accepts several options, which we list below.

+--------------------------+---------------------------------------------------------------------------------+
| :py:obj:`season`         | | The season number. For :py:obj:`K2`, this is the campaign number. Note that   |
|                          | | fractional seasons are allowed (i.e., 6.0). Default is 0                      |
+--------------------------+---------------------------------------------------------------------------------+
| :py:obj:`model`          | | The :py:obj:`everest` model name. Default :py:class:`nPLD`                    |
+--------------------------+---------------------------------------------------------------------------------+
| :py:obj:`-m` `mission`   | | The mission name (:py:obj:`k2` | :py:obj:`kepler` | :py:obj:`tess`).          |
|                          | | Default :py:obj:`k2`                                                          |
+--------------------------+---------------------------------------------------------------------------------+
| :py:obj:`-s`             | | Check the status of the short cadence runs.                                   |
+--------------------------+---------------------------------------------------------------------------------+
| :py:obj:`-i`             | | Check the status of the injection runs.                                       |
+--------------------------+---------------------------------------------------------------------------------+

.. raw:: html

  <script>
    (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
    (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
    m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
    })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

    ga('create', 'UA-47070068-2', 'auto');
    ga('send', 'pageview');

  </script>