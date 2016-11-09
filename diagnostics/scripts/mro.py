from types import MethodType

class special(object):
  def x(self):
    print(self.y)

class PrimaryModel(object):
  def __init__(self, *args, **kwargs):
    print("Init Primary")

class SecondaryModel(PrimaryModel):
  def __init__(self, *args, **kwargs):
    print("Init Secondary")

class StandardPLD(special):   
  def __init__(self, *args, **kwargs):
    super(StandardPLD, self).__init__(self, *args, **kwargs)
    print("Init Standard")
    self.y = 10
    
class NeighborPLD(object):
  def __init__(self, *args, **kwargs):
    super(NeighborPLD, self).__init__(self, *args, **kwargs)
    print("Init Neighbor")
    
# User-facing

def Model(method = StandardPLD, parent = None):
  '''
  
  '''
  
  if parent is None:
    model = PrimaryModel
    cname = method.__name__[0].lower()
    pname = None
  else:
    model = SecondaryModel
    cname = "%s(%s)" % (method.__name__[0].lower(), parent.__name__[0].lower())
    pname = parent.__name__
  
  class PLD(method, model):
    @property
    def name(self):
      return cname
    @property
    def parent(self):
      return pname
  return PLD()

Model(StandardPLD, StandardPLD).x()