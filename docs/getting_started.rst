Getting Started
===============

Project Description
-------------------
:py:mod:`everest` is a suite of Python routines to detrend **K2** light curves
using pixel level decorrelation (PLD) and gaussian processes (GPs).

Installation
------------
Clone the git repository with

.. code-block:: bash

   git clone --recursive https://github.com/rodluger/everest

to ensure you clone all the submodules. In the :py:mod:`everest.pysyzygy` subfolder, compile the
transit code by running `make`.

You'll also have to have the *astroML* and *george* Python packages installed. It is
recommended that you run :py:mod:`everest` under Python 3.