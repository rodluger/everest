The Catalog
===========

.. raw:: html

  <div class="admonition note">
  <p class="first admonition-title">Search the EVEREST Catalog</p>
  <p></p>
  <div class="last"><div role="search">
  <form id="target-search-form" class="wy-form" action="target_search.html" method="get">
      <div style="margin-bottom:1em;">
      Target ID: &nbsp;<input type="text" name="id" placeholder="e.g., 201367065" style = "height: 1.75em;">
      </div>
      Mission:  &nbsp;<input type="radio" name="mission" value="k2" checked>K2
                &nbsp;<input type="radio" name="mission" value="kepler" disabled><span style="color: #cccccc;">Kepler</span>
                &nbsp;<input type="radio" name="mission" value="tess" disabled><span style="color: #cccccc;">TESS</span>
  </form><br>
  </div></div></div>

The **EVEREST** catalog is available `here <https://archive.stsci.edu/prepds/everest/>`_
and includes ``.fits`` files with the de-trended light curve data as well as ``.pdf``
data validation summaries associated with each EPIC target. Follow the links 
below for detailed information about the catalog. Note that it is **highly recommended** 
that users access the catalog through the
:doc:`everest interface <user>`, as this allows for post-processing of
the light curves with custom masks.

.. toctree::
   :maxdepth: 2   
   
   fitsfiles
   dvsfigs
   issues