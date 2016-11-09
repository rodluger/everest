from everest2 import StandardPLD
import matplotlib.pyplot as pl


c2 = [204611221, 204353339, 203008341, 204856827]
c3 = [205998445, 206137557, 206143176, 206066800]

fig, ax = pl.subplots(4,2, figsize = (4,8))
for i, EPIC in enumerate(c3):
  StandardPLD(EPIC, clobber = True, ax = ax[i])
fig.savefig('debug.pdf')

# - # - #

def debug_poss(self):
  '''
  
  '''
  
  # Get colormap
  plasma = pl.get_cmap('plasma')
  plasma.set_bad('w')
  
  # Get aperture contour
  def PadWithZeros(vector, pad_width, iaxis, kwargs):
    vector[:pad_width[0]] = 0
    vector[-pad_width[1]:] = 0
    return vector
  ny, nx = self.pixel_images[0].shape
  contour = np.zeros((ny,nx))
  contour[np.where(self.aperture)] = 1
  contour = np.lib.pad(contour, 1, PadWithZeros)
  highres = zoom(contour, 100, order = 0, mode='nearest') 
  extent = np.array([-1, nx, -1, ny])
  
  # Plot first, mid, and last TPF image
  image = self.pixel_images[1]
  ax = self.debug_ax[0]
  ax.imshow(image, aspect = 'auto', interpolation = 'nearest', cmap = plasma)
  ax.contour(highres, levels=[0.5], extent=extent, origin='lower', colors='r', linewidths=1)
  ax.axis('off')       
  ax.set_xlim(-0.7, nx - 0.3)
  ax.set_ylim(-0.7, ny - 0.3)
  for source in self.nearby:
    ax.annotate('%.1f' % source['mag'], 
                xy = (source['x'] - source['x0'], source['y'] - source['y0']), 
                ha = 'center', va = 'center', size = 6, color = 'w', fontweight = 'bold')    
    
  # Plot external image
  ax = self.debug_ax[1]
  from .missions.k2 import GetHiResImage
  image = GetHiResImage(self.ID)
  ax.imshow(image, aspect = 'auto', extent = (-0.5, nx - 0.5, -0.5, ny - 0.5), interpolation = 'bicubic', cmap = plasma)
  ax.contour(highres, levels=[0.5], extent=extent, origin='lower', colors='r', linewidths=1)
  ax.axis('off')
  ax.set_xlim(-0.7, nx - 0.3)
  ax.set_ylim(-0.7, ny - 0.3)
