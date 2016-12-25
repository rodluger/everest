from everest.missions.k2.aux import GetK2Campaign
from everest.missions.k2 import TargetDirectory
import os
import shutil

keep = GetK2Campaign(6.0, epics_only = True)

for campaign in [5, 6]:
  all = GetK2Campaign(campaign, epics_only = True)
  sc = GetK2Campaign(campaign, epics_only = True, cadence = 'sc')
  lc = list(set(all) - set(sc) - set(keep))
  for i, star in enumerate(lc):
    print("C%02d: %d/%d..." % (campaign, i + 1, len(lc)))
    infile = os.path.join(TargetDirectory(star, campaign), 'data.npz')
    outpath = os.path.join('/lolo/archive/hyak/vsm/rodluger/everest2/', 'c%02d' % campaign, 
              ('%09d' % star)[:4] + '00000', 
              ('%09d' % star)[4:])
    if not os.path.exists(outpath):
      os.makedirs(outpath)
    outfile = os.path.join(outpath, 'data.npz')
    shutil.move(infile, outfile)