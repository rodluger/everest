<p align="center">
  <img width = "300" src="http://staff.washington.edu/rodluger/everest/_images/everest.png"/>
</p>
<p align="center">
  <a href="https://travis-ci.org/rodluger/everest/"><img src="https://travis-ci.org/rodluger/everest.svg?branch=master"/></a>
  <a href="http://arxiv.org/abs/1607.00524"><img src="https://img.shields.io/badge/arXiv-1607.00524-blue.svg?style=flat"/></a>
  <a href="https://raw.githubusercontent.com/rodluger/everest/master/LICENSE"><img src="https://img.shields.io/badge/license-MIT-brightgreen.svg"/></a>
  <a href="http://staff.washington.edu/rodluger/everest"><img src="https://img.shields.io/badge/read-the_docs-blue.svg?style=flat"/></a>
  <a href="https://archive.stsci.edu/prepds/everest/"><img src="https://img.shields.io/badge/MAST-lightcurves-brightgreen.svg?style=flat"/></a>
</p>

<div align="justify">
<b>E</b>PIC <b>V</b>ariability <b>E</b>xtraction and <b>R</b>emoval for <b>E</b>xoplanet <b>S</b>cience <b>T</b>argets: A pipeline for de-trending <b>K2</b> light curves with pixel level decorrelation and Gaussian processes. Here you'll find the Python code used to generate the <b>EVEREST</b> catalog, as well as tools for accessing and interacting with the de-trended light curves.
<br/><br/>

To install the latest <b>EVEREST</b> release (2.0):
<br/><br/>
<pre><code>pip install everest-pipeline</code></pre>
Note that <b>EVEREST</b> depends on <b>george</b>, which requires the <b>Eigen3</b> package. If you don't have <b>Eigen3</b>, follow the instructions <a href="http://dan.iel.fm/george/current/user/quickstart/">here</a> to get <b>george</b> set up.
<br/><br/>
You can also install the current development version of <b>EVEREST</b> from source:
<br/><br/>
<pre><code>git clone https://github.com/rodluger/everest
cd everest
python setup.py install --user</code></pre>
For more information on how to install and use <b>EVEREST</b>, please refer to the <a href="http://staff.washington.edu/rodluger/everest">documentation</a>. And if you're looking for the previous version of <b>EVEREST</b>, check out the latest <a href="https://github.com/rodluger/everest/tree/1.0.5">v1</a> release.
</div>
<br>
