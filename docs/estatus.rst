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