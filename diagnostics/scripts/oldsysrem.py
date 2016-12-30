def _OldGetCBVs(campaign, module = None, model = 'nPLD', nrec = 5, clobber = False, plot = True, **kwargs):
  '''
  DEPRECATED!
  
  '''
  
  # Initialize logging?
  if len(logging.getLogger().handlers) == 0:
    InitLog(file_name = None, screen_level = logging.DEBUG)
  
  # All modules?
  if module is None:
    
    if not plot:
      for module in range(2, 25):
        X = GetCBVs(campaign, module = module, model = model, clobber = clobber, **kwargs)
    
    else:
    
      # We're going to plot the CBVs on the CCD
      fig = [None] + [None for n in range(1, 1 + nrec)]
      ax = [None] + [None for n in range(1, 1 + nrec)]
      for n in range(1, nrec + 1):
        fig[n], ax[n] = pl.subplots(5, 5, figsize = (9, 9))
        fig[n].subplots_adjust(wspace = 0.025, hspace = 0.025)
        ax[n] = [None] + list(ax[n].flatten())
        for axis in [ax[n][1], ax[n][5], ax[n][21], ax[n][25]]:
          axis.set_visible(False)
        for i in range(1, 25):
          ax[n][i].set_xticks([])
          ax[n][i].set_yticks([])
          ax[n][i].annotate('%02d' % i, (0.5, 0.5), 
                            va = 'center', ha = 'center',
                            xycoords = 'axes fraction',
                            color = 'k', fontsize = 60, alpha = 0.05)
          ax[n][i].margins(0.1, 0.1)
        
      # Get the CBVs
      for module in range(2, 25):
        X = GetCBVs(campaign, module = module, model = model, clobber = clobber, **kwargs)
        if X is not None:
        
          # Get the timeseries info
          infofile = os.path.join(EVEREST_DAT, 'k2', 'cbv', 'c%02d' % campaign, str(module), model, 'info.npz')
          info = np.load(infofile)
          time = info['time']
          nstars = info['nstars']
          breakpoints = info['breakpoints']
          
          # Plot the CBVs
          for b in range(len(breakpoints)):
            inds = GetChunk(time, breakpoints, b)
            for n in range(1, min(6, X.shape[1])):
              ax[n][module].plot(time[inds], X[inds,n])
              if b == 0:
                ax[n][module].annotate(nstars, (0.01, 0.02), va = 'bottom', ha = 'left',
                                       xycoords = 'axes fraction', color = 'k', fontsize = 8,
                                       alpha = 0.5)
                
      # Save the figures
      for n in range(1, 1 + nrec):
        figname = os.path.join(EVEREST_DAT, 'k2', 'cbv', 'c%02d' % campaign, model + '_%02d.pdf' % n)
        fig[n].suptitle('CBV #%02d' % n, fontsize = 18, y = 0.94)
        fig[n].savefig(figname, bbox_inches = 'tight')
        pl.close(fig[n])
    
    return
  
  log.info('Computing CBVs for campaign %d, module %d...' % (campaign, module))
    
  # Output path
  path = os.path.join(EVEREST_DAT, 'k2', 'cbv', 'c%02d' % campaign, str(module), model)
  if not os.path.exists(path):
    os.makedirs(path)
  
  # Get the design matrix
  xfile = os.path.join(path, 'X.npz')
  if clobber or not os.path.exists(xfile):
    
    # Get the light curves
    log.info('Obtaining light curves...')
    lcfile = os.path.join(path, 'lcs.npz')
    infofile = os.path.join(path, 'info.npz')
    if clobber or not os.path.exists(lcfile):
      try:
        time, breakpoints, fluxes, errors, kpars = GetStars(campaign, module, model = model, **kwargs)
      except AssertionError:
        np.savez(lcfile, time = None, breakpoints = None, fluxes = None, errors = None, kpars = None)
        np.savez(infofile, time = None, breakpoints = None, nstars = None)
        np.savez(xfile, X = None)
        return None
      np.savez(lcfile, time = time, breakpoints = breakpoints, fluxes = fluxes, errors = errors, kpars = kpars)
      np.savez(infofile, time = time, breakpoints = breakpoints, nstars = len(fluxes))
    else:
      lcs = np.load(lcfile)
      time = lcs['time']
      breakpoints = lcs['breakpoints']
      fluxes = lcs['fluxes']
      errors = lcs['errors']
      kpars = lcs['kpars']
    
    # Compute the design matrix  
    log.info('Running SysRem...')
    X = np.ones((len(time), 1 + nrec))
    
    # Loop over the segments
    new_fluxes = np.zeros_like(fluxes)
    for b in range(len(breakpoints)):
      
      # Get the current segment's indices
      inds = GetChunk(time, breakpoints, b)
    
      # Update the error arrays with the white GP component
      for j in range(len(errors)):
        errors[j] = np.sqrt(errors[j] ** 2 + kpars[j][0] ** 2)
    
      # Get de-trended fluxes
      X[inds,1:] = SysRem(time[inds], fluxes[:,inds], errors[:,inds], nrec = nrec, **kwargs).T
      
    # Save
    np.savez(xfile, X = X)
  
  else:
    
    # Load from disk
    X = np.load(xfile)['X'][()]
  
  if X is not None: 
    # Ensure we only return as many as we asked for 
    return X[:,:nrec + 1]
  else:
    return X
