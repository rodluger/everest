EVEREST
-------

Pixel Level Decorrelation of K2 Light curves


To-do
=====

- Easy EVEREST calling for masked transits.
- Take median of unsaturated background!
- Make Eigen3 installation easier (copy george?)
- Import astroML
- MacOS backend doesn't work

- Re-run injections; re-run precision plots.
- Plot the autocorrelation function and scatter in the binned residuals 
  as a function of bin size.  Of course stellar variability that remains 
  unmodeled can cause this scatter, but it will be interesting to compare 
  with the other pipelines.
- It might be fair to show an example transit light curve in which EVEREST 
  does worse than K2SFF (and perhaps instructive!)
- It might be good to show an example of all of the steps applied to a single EPIC.  
  
- Run a quick version of the pipeline, no optimization!
- Flag outliers in the output file
- Cite `this paper <https://arxiv.org/abs/1604.07442>`_ on asteroseismology
- Better saturation and crowding metrics
- Test **EVEREST** on short-cadence data
- Compute gain curve to correct saturated pixels
- Interpolate over each module with basis vectors from bright stars (paper 2)