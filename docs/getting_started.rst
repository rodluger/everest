Installation
============

Using :py:mod:`everest` is easy once you get it set up.
Clone the git repository with

.. code-block:: bash

   git clone --recursive https://github.com/rodluger/everest

to ensure you clone all the submodules. In the :py:mod:`everest.pysyzygy` subfolder, compile the
transit code by running `make`.

You'll also have to have the *astroML* and *george* Python packages installed. You can
do this easily with `pip <https://pypi.python.org/pypi/pip>`_:

.. code-block:: bash

   pip install astroML
   pip install george --global-option=build_ext --global-option=-I/path/to/eigen3

where :obj:`/path/to/eigen3` should be replaced with the path to the *parent directory*
containing the **Eigen3** distribution. You can download **Eigen3** from 
`here <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ (no installation
necessary!) More information on installing :py:mod:`george` at
`Dan Foreman-Mackey's page <http://dan.iel.fm/george/current/user/quickstart/>`_.

You also need *scipy*, *numpy*, *matplotlib*, and some other standard Python packages --
these should be included in all of the major Python distributions (I'm using
`Anaconda <https://www.continuum.io/downloads>`_). Note that :py:mod:`everest` was developed
using **Python 3**.