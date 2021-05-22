#include "math.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_cblas.h"
#include "gsl/gsl_eigen.h"
#include "gsl_utils.h"

void invert_matrix(double* arr,double* out,int N)
{
  gsl_matrix *mat = gsl_matrix_alloc(N,N);
  int i,j;
  for (i=0 ; i<N ; i++) 
    for (j=0 ; j<N ; j++)
      gsl_matrix_set(mat,i,j,arr[i*N+j]);

  int signum;
  gsl_matrix *inv = gsl_matrix_alloc(3,3);
  gsl_permutation *perm = gsl_permutation_alloc(3);

  gsl_linalg_LU_decomp(mat,perm,&signum);
  gsl_linalg_LU_invert(mat,perm,inv);

  for (i=0 ; i<N ; i++) 
    for (j=0 ; j<N ; j++)
      out[i*N+j] = gsl_matrix_get(inv,i,j);

  gsl_matrix_free(mat);
  gsl_matrix_free(inv);
}

// A*B = C

void matrix_product(double* A,double* B,double* C,int N)
{
  int i,j,k;
  for (i=0 ; i<N ; i++)
    for (j=0 ; j<N ; j++)
    {
      C[i*N+j] = 0;
      for (k=0 ; k<N ; k++)
        C[i*N+j] += A[i*N+k]*B[k*N+j];
    }
}

void eigen_solve(double* arr,double* evals,double* evecs,int N)
{
  int i,j;
  gsl_matrix *eigen = gsl_matrix_alloc(N,N);
  gsl_matrix_view mat = gsl_matrix_view_array((double*)arr,N,N);
  gsl_matrix_memcpy(eigen,&mat.matrix);
  gsl_vector_complex *evalsC = gsl_vector_complex_alloc(N);
  gsl_matrix_complex *evecsC = gsl_matrix_complex_alloc(N,N);
     
  gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(N);
  gsl_eigen_nonsymmv(eigen,evalsC,evecsC,w); // results are complex matrix, vector
  gsl_eigen_nonsymmv_free(w);
  gsl_eigen_nonsymmv_sort(evalsC,evecsC,GSL_EIGEN_SORT_ABS_DESC); // largest mag first
       
  for (i=0 ; i<N ; i++)
  {
    evals[i] = GSL_REAL(gsl_vector_complex_get(evalsC,i));
    for (j=0 ; j<N ; j++) evecs[i*N+j] = GSL_REAL(gsl_matrix_complex_get(evecsC,i,j));
  }
}

void mul_mat_vec(double* mat,double* vec,double* out,int N)
{
  int i,j;
  for (i=0 ; i<3 ; i++)
  {
    out[i] = 0;
    for (j=0 ; j<3 ; j++) out[i] += mat[i*N+j]*vec[j];
  }
}

void transpose(double* X,double* Xt,int N)
{
  int i,j;
  for (i=0 ; i<N ; i++)
    for (j=0 ; j<N ; j++)
      Xt[i*N+j] = X[j*N+i];
}

double det33(double* a)
{
  double s=0;
  for (int i=0 ; i<3 ; i++)
  {
    double p=1,q=1;
    for (int j=0 ; j<3 ; j++) 
    {
      p *= a[3*j+(i+j)%3];
      q *= a[3*j+(3+i-j)%3];
    }
    s += p-q;
  }
  return s;
}

double det(double* arr,int N)
{
  gsl_matrix *mat = gsl_matrix_alloc(N,N);
  int i,j;
  for (i=0 ; i<N ; i++) 
    for (j=0 ; j<N ; j++)
      gsl_matrix_set(mat,i,j,arr[i*N+j]);

  int signum;
  gsl_matrix *inv = gsl_matrix_alloc(3,3);
  gsl_permutation *perm = gsl_permutation_alloc(3);

  gsl_linalg_LU_decomp(mat,perm,&signum);
  double x=gsl_linalg_LU_det(mat,signum);

  gsl_matrix_free(mat);
  gsl_matrix_free(inv);

  return x;
}

void eye(double* A,int N)
{
  int i,j;
  for (i=0 ; i<N ; i++)
    for (j=0 ; j<N ; j++)
    {
      if (i==j) A[i*N+j] = 1;
      else A[i*N+j] = 0;
    }
}

// A = U.S.Vt

void SVD(double* _A,double* _U,double* _S,double* _Vt,int N)
{
  gsl_matrix *A=gsl_matrix_alloc(N,N);
  gsl_matrix *V=gsl_matrix_alloc(N,N);
  gsl_vector *S=gsl_vector_alloc(N);
  gsl_vector *work=gsl_vector_alloc(N);

  int i,j;
  for (i=0 ; i<N ; i++)
    for (j=0 ; j<N ; j++) gsl_matrix_set(A,i,j,_A[i*N+j]);
  gsl_linalg_SV_decomp(A,V,S,work); // A = U.S.Vt, A replaced with U

  for (i=0 ; i<N ; i++)
    for (j=0 ; j<N ; j++) 
    {
      _U[i*N+j] = gsl_matrix_get(A,i,j);
      if (i==j) _S[i*N+j] = gsl_vector_get(S,i);
      else _S[i*N+j] = 0;
      _Vt[i*N+j] = gsl_matrix_get(V,j,i); // transpose!
    }

  gsl_matrix_free(V);
  gsl_vector_free(S);
  gsl_matrix_free(A);
  gsl_vector_free(work);
}

void vec_diff(double* A,double* B,double* C)
{
  for (int i=0 ; i<3 ; i++) C[i] = A[i]-B[i];
}

void vec_add(double* A,double* B,double* C)
{
  for (int i=0 ; i<3 ; i++) C[i] = A[i]+B[i];
}

void get_center(double* X,int N,double* cen)
{
  for (int i=0 ; i<3 ; i++) cen[i] = 0;
  for (int i=0 ; i<N ; i++)
    for (int j=0 ; j<3 ; j++) cen[j] += X[i*3+j];
  for (int i=0 ; i<3 ; i++) cen[i] /= float(N);
}

double dist2(double* A,double* B)
{
  double d=0;
  for (int i=0 ; i<3 ; i++)
    d += (A[i]-B[i])*(A[i]-B[i]);
  return d;
}

void covariance(double cov[3][3],double* A,double* B,int N)
{
  for (int i=0 ; i<3 ; i++)
    for (int j=0 ; j<3 ; j++) 
    {
      cov[i][j] = 0;
      for (int k=0 ; k<N ; k++) cov[i][j] += A[k*3+i]*B[k*3+j];
    }
}

// assume coords are pre-centered
// rotation of first onto second

void superpose(double *X,double *Y,int N,double* rot)
{
  double cov[3][3],U[3][3],V[3][3],Vt[3][3],S[3][3],W[3][3],Wt[3][3],D[3][3],P[3][3];
  covariance(cov,X,Y,N);
  SVD((double*)cov,(double*)V,(double*)S,(double*)Wt,3); // M = V.S.Wt
  transpose((double*)Wt,(double*)W,3);
  transpose((double*)V,(double*)Vt,3);
  eye((double*)D,3);
  //if (det((double*)cov,3)<0) D[2][2] = -1;
  if (det33((double*)cov)<0) D[2][2] = -1;
  matrix_product((double*)W,(double*)D,(double*)P,3); // P is an intermediate
  matrix_product((double*)P,(double*)Vt,(double*)rot,3); // rot = W.D.Vt, rotation matrix
}

double superpose_rms(double *X,double *Y,int N,double* rot)
{
  superpose(X,Y,N,rot);
  double r=0,transformed[3];
  for (int k=0 ; k<N ; k++)
  {
    mul_mat_vec(rot,&X[3*k],transformed,3);
    r += dist2(transformed,&Y[3*k]);
  }
  r = sqrt(r/(float)N);
  return r;  
}


