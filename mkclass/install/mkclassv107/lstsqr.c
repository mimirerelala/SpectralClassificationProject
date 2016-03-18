// Linear Least Squares routine
#include <stdio.h>
#include <math.h>
#include "util.h"

void lstsqr(float *X, float *Y, int nq, float wlow, float whigh, float *a, float *b)
{
  int i,j,k,nlow,nhigh;
  float *x,*y;
  float S,Sx,Sy,Sxx,Sxy,D;

  nlow = 0;
  
  for(i=0;i<nq;i++) {
    if(X[i] <= wlow) nlow = i;
    if(X[i] <= whigh) nhigh = i;
  }

  k = nhigh - nlow + 1;

  x = vector(0,k);
  y = vector(0,k);

  j = 0;
  for(i=nlow;i<=nhigh;i++) {
    x[j] = X[i];
    y[j] = Y[i];
    j++;
  }

  S = Sx = Sy = Sxx = Sxy = 0.0;
  for(i=0;i<=k;i++) {
    S += 1.0;
    Sx += x[i];
    Sy += y[i];
    Sxx += x[i]*x[i];
    Sxy += x[i]*y[i];
  }

  D = S*Sxx - Sx*Sx;

  *a = (Sxx*Sy-Sx*Sxy)/D;
  *b = (S*Sxy-Sx*Sy)/D;

  free_vector(x,0,k);
  free_vector(y,0,k);
  return;
}
