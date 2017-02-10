Bulk Download
=============

For users interested in downloading de-trended light curves in bulk, tarballs
for each campaign are available on `MAST <http://archive.stsci.edu/missions/hlsp/everest/v2/bundles/>`_.
If you then wish to interact with these light curves using the :doc:`user interface <ui>`,
untar the light curves and place them in the directory

.. code-block:: bash

   ~/.everest2/k2/

organized by campaign and sub-campaign following the same directory structure as on 
`MAST <http://archive.stsci.edu/missions/hlsp/everest/v2/bundles/>`_. For instance,
the :py:obj:`.fits` and :py:obj:`.pdf` files for the campaign 1 star EPIC 201367065 
should be placed in

.. code-block:: bash

   ~/.everest2/k2/c01/201300000/67065/

Then, you can load the cached light curve in Python by running

.. code-block:: python

   import everest
   star = everest.Everest(201367065)

.. raw:: html

  <script>
    (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
    (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
    m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
    })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

    ga('create', 'UA-47070068-3', 'auto');
    ga('send', 'pageview');

  </script>
