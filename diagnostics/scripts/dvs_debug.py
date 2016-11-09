import sys
sys.path.insert(1, '/Users/rodrigo/src/everest2')
from everest2.dvs import DVS
import matplotlib.pyplot as pl
dvs = DVS()
dvs.title()
dvs.footer()
dvs.top_left()
[dvs.top_right() for i in range(4)]
[dvs.left() for i in range(4)]
[dvs.right() for i in range(3)]
pl.show()