import numpy as np
import os, glob

path = '/Users/rodrigo/src/everest/paper/scripts/CDPP'

# K2SC
for c in range(8):
  try:
    epic, _, cdpp = np.loadtxt(os.path.join(path, 'k2sc_C%02d.tsv' % c), unpack = True)
    epic = np.array(epic, dtype = int)
    with open('/Users/rodrigo/src/everest2/everest2/missions/k2/tables/c%02d_k2sc.tsv' % c, 'w') as outfile:
      for i in range(len(epic)):
        print("{:>09d} {:>15.3f}".format(epic[i], cdpp[i]), file = outfile)
  except:
    continue
    
# K2SFF
for c in range(8):
  epic, _, cdpp = np.loadtxt(os.path.join(path, 'k2sff_C%02d.tsv' % c), unpack = True)
  epic = np.array(epic, dtype = int)
  with open('/Users/rodrigo/src/everest2/everest2/missions/k2/tables/c%02d_k2sff.tsv' % c, 'w') as outfile:
    for i in range(len(epic)):
      print("{:>09d} {:>15.3f}".format(epic[i], cdpp[i]), file = outfile)

# Everest 1
for c in range(8):
  epic, _, cdpp = np.loadtxt(os.path.join(path, 'everest_C%02d.tsv' % c), unpack = True)
  epic = np.array(epic, dtype = int)
  with open('/Users/rodrigo/src/everest2/everest2/missions/k2/tables/c%02d_everest1.tsv' % c, 'w') as outfile:
    for i in range(len(epic)):
      print("{:>09d} {:>15.3f}".format(epic[i], cdpp[i]), file = outfile)