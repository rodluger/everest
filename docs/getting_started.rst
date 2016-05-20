Installation
============

Using :py:mod:`everest` is easy once you get it set up.
Clone the git repository with

.. code-block:: bash

   git clone --recursive https://github.com/rodluger/everest
   cd everest
   python setup.py install --user

If you don't have the :py:mod:`george` package installed, you might get an error
complaining that the installer can't find the **Eigen3** distribution. 
You can download the latest **Eigen3** tarball from 
`here <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_. Just untar it
in the directory of your choice and install :py:mod:`george` by telling it where
to look:

.. code-block:: bash

   pip install george --global-option=build_ext --global-option=-I/path/to/eigen3

where :obj:`/path/to/eigen3` should be replaced with the path to the *parent directory*
containing the **Eigen3** distribution. More information on installing :py:mod:`george` at
`Dan Foreman-Mackey's page <http://dan.iel.fm/george/current/user/quickstart/>`_.