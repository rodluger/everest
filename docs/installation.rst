Installation
============

As of version 2.0, :py:mod:`everest` is `pip <https://en.wikipedia.org/wiki/Pip_(package_manager)>`_ 
installable:

.. code-block:: bash

   pip install everest-pipeline

This should allow you to run the :doc:`command line utilities <ui>` as well
as to :py:obj:`import everest` in Python.

.. warning:: If you have the previous version of :py:obj:`everest` and/or \
             previous versions of :py:obj:`k2plr` and :py:obj:`pysyzygy` \
             installed, it might be a good idea to **completely** remove them before \
             attempting to install :py:obj:`everest 2.0`.

.. warning:: :py:obj:`everest` requires :py:obj:`george`, which in turn needs access \
             to the :py:obj:`Eigen3` linear algebra package. If you don't have that \
             installed, follow the instructions `here <http://dan.iel.fm/george/current/user/quickstart/>`_ 
             to get :py:obj:`george` set up.

Alternatively, the development version of :py:mod:`everest` is maintained on 
`github <https://github.com/rodluger/everest>`_.
You can install it by cloning the repository and running the setup script:

.. code-block:: bash

   git clone https://github.com/rodluger/everest
   cd everest
   python setup.py install --user

.. role:: python(code)
   :language: python

.. raw:: html

  <script>
    (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
    (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
    m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
    })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

    ga('create', 'UA-47070068-3', 'auto');
    ga('send', 'pageview');

  </script>