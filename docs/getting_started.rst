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