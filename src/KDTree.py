import math

def distsq((a,b,c),(d,e,f)):
  return (a-d)*(a-d)+(b-e)*(b-e)+(c-f)*(c-f)

# based on Andrew Moore's tutorial on kd-trees:
#   http://www.autonlab.org/autonweb/14665/version/2/part/5/data/moore-tutorial.pdf?branch=main&language=en

# example of how to use this:
#  kdtree = KDTree([x.co for x in receptor])
#  (dsq,coord) = kdtree.nearest(prot[j].co)
#  j = kdtree.index_of(coord)

class KDTree:

  # members: leaf (0/1), dim, val (float), cutpt (coord)
  #          left, right (KDTrees)
  #          mins,maxs (coords)

  def __init__(self,coords):
    self.INFY = 999999
    self.inverse = {}
    for i in xrange(len(coords)):
      self.inverse[coords[i][0],coords[i][1],coords[i][2]] = i
    if len(coords)==0: self.leaf = 1; return
    else: self.leaf = 0
    (self.mins,self.maxs) = self.compute_extremes(coords)
    diffs = [x-y for x,y in zip(self.maxs,self.mins)]
    if diffs[0]>=diffs[1] and diffs[0]>=diffs[2]: self.dim = 0
    elif diffs[1]>=diffs[0] and diffs[1]>=diffs[2]: self.dim = 1
    else: self.dim = 2
    self.val = (self.maxs[self.dim]+self.mins[self.dim])/2.0

    best = None
    for co in coords:
      Del = math.fabs(co[self.dim]-self.val)
      if best==None or Del<delbest: (best,delbest) = (co,Del)
    self.cutpt = best
    self.val = self.cutpt[self.dim]
    (lset,rset)=([],[])
    for co in coords:
      if co==best: continue # exclude cutpt from either sub-tree
      if co[self.dim]<=self.val: lset.append(co)
      else: rset.append(co)
    self.left = KDTree(lset)
    self.right = KDTree(rset)

  def num_nodes(self):
    if self.leaf==1: return 1
    return 1+self.left.num_nodes()+self.right.num_nodes()

  def num_leaves(self):
    if self.leaf==1: return 1
    return self.left.num_leaves()+self.right.num_leaves()

  def index_of(self,co): return self.inverse[co[0],co[1],co[2]]

  # returns INFY if no pt found within dmax, else returns distance-squared

  def nearest(self,targ):
    return self.nearest_(targ,self.mins,self.maxs,self.INFY)

  def nearest_(self,targ,mins,maxs,dmax):
    if self.leaf==1: return (self.INFY,None)

    (minN,minF) = (mins[0:3],mins[0:3])
    (maxN,maxF) = (maxs[0:3],maxs[0:3])
    if targ[self.dim]<=self.val: 
      (nearer,further) = (self.left,self.right)
      (maxN[self.dim],minF[self.dim]) = (self.val,self.val)
    else: 
      (nearer,further) = (self.right,self.left)
      (minN[self.dim],maxF[self.dim]) = (self.val,self.val)

    (dsq,output) = nearer.nearest_(targ,minN,maxN,dmax)
    if dsq<dmax: dmax = dsq
    if dmax<(targ[self.dim]-self.val)*(targ[self.dim]-self.val): 
      return (dsq,output)
    d2 = distsq(targ,self.cutpt)
    if d2<dsq:
      output = self.cutpt
      dmax = dsq = d2

    (d3,out) = further.nearest_(targ,minF,maxF,dmax)
    if d3<dsq: (dsq,output) = (d3,out)
    return (dsq,output)

  def near(self,coord):
    (d,pt) = self.nearest(coord)
    if d<DMESH*DMESH: return 1
    return 0

  def compute_extremes(self,coords):
    xs = map(lambda x: x[0],coords)
    ys = map(lambda x: x[1],coords)
    zs = map(lambda x: x[2],coords)
    mins = [min(xs),min(ys),min(zs)]
    maxs = [max(xs),max(ys),max(zs)]
    return (mins,maxs)

  def range_query_radius(self,co,rad):
    lows = map(lambda x: x-rad,co)
    highs = map(lambda x: x+rad,co)
    hits = self.range_query(lows,highs)
    return filter(lambda x: distsq(x,co)<rad*rad,hits)

  def range_query(self,lows,highs):
    if self.leaf==1: return []
    for i in range(3):
      if lows[i]>self.maxs[i]: return []
      if highs[i]<self.mins[i]: return []
    # there must be some overlap volume
    a = self.left.range_query(lows,highs)
    b = self.right.range_query(lows,highs)
    c = [self.cutpt]
    for i in range(3):
      if self.cutpt[i]<lows[i] or self.cutpt[i]>highs[i]: c = []
    return a+b+c

