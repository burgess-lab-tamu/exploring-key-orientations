import math,random

def stats(vals):
  sum,ss = 0,0
  for x in vals: sum += x; ss += x*x
  N = float(len(vals))
  mean = sum/N
  var = ss/N-mean*mean
  if var<0: var = 0
  stdev = math.sqrt(var)
  return mean,stdev

def corr(X,Y):
  muX,sdX = stats(X)
  muY,sdY = stats(Y)
  cX = [x-muX for x in X]
  cY = [y-muY for y in Y]
  s = sum([x*y for (x,y) in zip(cX,cY)])
  return s/(float(len(X))*sdX*sdY)

# return list of (item,cnt) sorted by counts

def popularity(lst):
  hash = {}
  for x in lst:
    if x not in hash: hash[x] = 0
    hash[x] += 1
  data = [(hash[x],x) for x in hash.keys()]
  data.sort(reverse=True)
  data = [(y,x) for (x,y) in data]
  return data

def frand_pm1(): return 2.*random.random()-1.

def dist((a,b,c),(d,e,f)): return math.sqrt(distsq((a,b,c),(d,e,f)))

def distsq((a,b,c),(d,e,f)):
  return (a-d)*(a-d)+(b-e)*(b-e)+(c-f)*(c-f)

def remove_duplicates(list):
  unique = []
  for x in list:
    if not x in unique: unique.append(x)
  return unique

def intersection(a,b):
  subset = []
  for i in a:
    if i in b: subset.append(i)
  return subset

def pdb_rec(anum,anam,rnam,rnum,chn,x,y,z,occ,bf):
  return "ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" % \
      (anum,anam,rnam,chn,rnum%10000,x,y,z,occ,bf)
    
def diff(a,b): return [a[k]-b[k] for k in xrange(3)]

def vsum(a,b): return [a[k]+b[k] for k in xrange(3)]
  
def prod(mat,v):
  w = [0 for i in xrange(3)]
  for i in xrange(3):
    for j in xrange(3): w[i] += mat[i][j]*v[j]
  return w

def cross(u,v):
  a = u[1]*v[2]-u[2]*v[1]
  b = u[2]*v[0]-u[0]*v[2]
  c = u[0]*v[1]-u[1]*v[0]
  return [a,b,c]

def my_dot(a,b): return sum([x*y for x,y in zip(a,b)])

def mult_mat_vec(M,v):
  w = []
  for i in range(3):
    z = 0
    for j in range(3): z += M[i][j]*v[j]
    w.append(z)
  return w

def angle(u,v):
  a2 = my_dot(u,u)
  b2 = my_dot(v,v)
  c2 = my_dot(u,v)
  if a2<=0 or b2<=0: return 0
  d = c2/math.sqrt(a2*b2)
  if d<-1: d = -1
  if d>1: d = 1
  d = math.acos(d)
  return 180*d/3.14159

def angle3(a,b,c):
  u = diff(a,b)
  v = diff(c,b)
  return angle(u,v)

# from p. 75 in Altmann

def rotation_matrix(dir,ang): # in degrees

  (nx,ny,nz) = (dir[0],dir[1],dir[2])
  mag = math.sqrt(nx*nx+ny*ny+nz*nz)
  (nx,ny,nz) = (nx/mag,ny/mag,nz/mag) # normalize to unit vector
  ang *= 3.1415927/180.0
  sn = math.sin(ang)
  sq = math.sin(0.5*ang)
  sq *= sq

  m = [[0 for i in xrange(3)] for j in xrange(3)]

  m[0][0] = 1-2*(ny*ny+nz*nz)*sq
  m[0][1] = -nz*sn+2*nx*ny*sq
  m[0][2] = ny*sn+2*nx*nz*sq

  m[1][0] = nz*sn+2*nx*ny*sq
  m[1][1] = 1-2*(nx*nx+nz*nz)*sq
  m[1][2] = -nx*sn+2*ny*nz*sq

  m[2][0] = -ny*sn+2*nx*nz*sq
  m[2][1] = nx*sn+2*ny*nz*sq
  m[2][2] = 1-2*(nx*nx+ny*ny)*sq

  return m

