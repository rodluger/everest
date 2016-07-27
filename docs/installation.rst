Installation
============

Using :py:mod:`everest` is easy once you get it set up. The trickiest step
will likely be getting the :py:mod:`george` package installed, as it requires
the **Eigen3** distribution. If you don't have :py:mod:`george`, follow the
instructions on `Dan Foreman-Mackey's page <http://dan.iel.fm/george/current/user/quickstart/>`_.
Then, to install :py:mod:`everest`, all you need to do is run

.. code-block:: bash

   git clone --recursive https://github.com/rodluger/everest
   cd everest
   python setup.py install --user

.. note:: :py:mod:`everest` was coded using :py:mod:`george 0.2.1`. If you have the development \
          version of :py:mod:`george` installed, you'll get tons of errors. Consider installing \
          :py:mod:`george 0.2.1` in a separate directory and adding it to your :py:obj:`$PATH` \
          before calling :py:mod:`everest`.

.. raw:: html

  <script>
    (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
    (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
    m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
    })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

    ga('create', 'UA-47070068-2', 'auto');
    ga('send', 'pageview');
  </script>