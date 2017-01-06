:py:mod:`estats` - De-trending Statistics
-----------------------------------------

The :py:mod:`estats` command accepts several options, which we list below.

+--------------------------+---------------------------------------------------------------------------------+
| :py:obj:`season`         | | The season number. For :py:obj:`K2`, this is the campaign number. Default 0   |
+--------------------------+---------------------------------------------------------------------------------+
| :py:obj:`model`          | | The :py:obj:`everest` model name. Default :py:class:`nPLD`                    |
+--------------------------+---------------------------------------------------------------------------------+
| :py:obj:`compare_to`     | | The :py:obj:`everest` model or pipeline to compare against.                   |
|                          | | Default :py:obj:`everest1`                                                    |
+--------------------------+---------------------------------------------------------------------------------+
| :py:obj:`-s`             | | Plot the short cadence versus long cadence CDPP statistics.                   |
|                          | | If no campaign is specified, shows all campaigns                              |
+--------------------------+---------------------------------------------------------------------------------+
| :py:obj:`-i`             | | Plot the transit injection/recovery results.                                  |
+--------------------------+---------------------------------------------------------------------------------+
| :py:obj:`-p`             | | Plot the CDPP comparison for all confirmed planet hosts                       |
+--------------------------+---------------------------------------------------------------------------------+