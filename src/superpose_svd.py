import sys,math,os
import scipy,scipy.linalg

def dot(X,Y): return scipy.dot(X,Y)

# returns transform of second onto first

def get_transform(co1,co2):
  co1,co2 = co2,co1 # this code actually transforms first onto second
  N = len(co1)
  if N!=len(co2): print "error: superpose() requires same number of coordinates"; sys.exit(0)

  center1,center2 = [],[]
  for i in range(3): center1.append(sum([x[i] for x in co1])/float(N))
  for i in range(3): center2.append(sum([x[i] for x in co2])/float(N))

  A = [[x[j]-center1[j] for j in range(3)] for x in co1] # centered
  B = [[x[j]-center2[j] for j in range(3)] for x in co2] # centered

  U = superpose_matrix_svd(A,B)
  return (U,center1,center2)

def apply_transform(coords,U,center1,center2):
  A = [[x[j]-center1[j] for j in range(3)] for x in coords]
  return [(dot(U,x)+center2).tolist() for x in A]

# input 2 sets of parallel coords (lists of float-triples)
# output list of co2 transformed onto co1

def superpose_svd(co1,co2):
  co1,co2 = co2,co1 # this code actually transforms first onto second
  N = len(co1)
  if N!=len(co2): print "error: superpose() requires same number of coordinates"; sys.exit(0)

  center1,center2 = [],[]
  for i in range(3): center1.append(sum([x[i] for x in co1])/float(N))
  for i in range(3): center2.append(sum([x[i] for x in co2])/float(N))

  A = [[x[j]-center1[j] for j in range(3)] for x in co1] # centered
  B = [[x[j]-center2[j] for j in range(3)] for x in co2] # centered

  U = superpose_matrix_svd(A,B)
  return [(dot(U,x)+center2).tolist() for x in A]

# assume coords are already centered
# return rotation matrix
# see http://en.wikipedia.org/wiki/Kabsch_algorithm

def superpose_matrix_svd(A,B):
  N = len(A)
  cov = scipy.zeros((3,3))
  for i in range(3):
    for j in range(3):
      cov[i,j] = 0
      for k in range(N): cov[i,j] += A[k][i]*B[k][j]
  V,S,Wt = scipy.linalg.svd(cov) # S is vec, use diag(S)
  #print "svd = V.S.Wt\n",dot(V,dot(diag(S),Wt)); print "cov = \n",cov
  Vt = scipy.transpose(V)
  W = scipy.transpose(Wt)
  D = scipy.eye(3)
  if scipy.linalg.det(cov)<0: D[2,2] = -1. # shouldn't be necessary if symmetric
  U = dot(W,dot(D,Vt)) # matmult(W,matmult(D,Vt)) # W*D*Vt
  return U

