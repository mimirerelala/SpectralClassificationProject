// Determines a spectral type based on the metallic-line spectrum.  Mostly useful for F and G-type stars
#include <stdio.h>
#include <math.h>
#include "util.h"
void sp2class(float spcode, float lumcode, float metcode, 
	      float *x, float *y, int *n);
void lstsqr(float *X, float *Y, int nq, float wlow, float whigh, float *a, float *b);
float CHband2(float spcode, float lumcode);

float sptmetal(float spcode)
{
  extern float *X,*Y;
  extern float Lumcode;
  extern int N;
  extern int Nq;
  float *x,*y;
  float chi,sumy,sumY,avey,aveY,sigy,sigY;
  float sa,sb,pa,pb;
  float chi2a = 0.0;
  float chi2b = 0.0;
  float chi2 = 0.0;
  float chi2b1,chi2b2,chi2b3;
  int i,n,nq;
  int flag1,flag2,flag3;
  extern FILE *Log;

  if(spcode > 36.0) return(1.0e+10);

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,Lumcode,0.0,x,y,&n);


  flag1 = flag2 = flag3 = 0;
  nq = 0;
  sumy = sumY = 0.0;
  chi2a = 0.0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4135.0 && x[i] <= 4315.0) {
      if(flag1 == 0) {
        lstsqr(X,Y,n,4135.0,4315.0,&pa,&pb);
        lstsqr(x,y,n,4135.0,4315.0,&sa,&sb);
        flag1 = 1;
      }
      chi = y[i]/(sa+x[i]*sb) - Y[i]/(pa+x[i]*pb);
      chi2a += chi*chi;
      nq++;
    }
    else if(x[i] >= 4410.0 && x[i] <= 4600.0) {
      if(flag2 == 0) {
	lstsqr(X,Y,n,4410.0,4600.0,&pa,&pb);
	lstsqr(x,y,n,4410.0,4600.0,&sa,&sb);
	flag2 = 1;
      }
      chi = y[i]/(sa+x[i]*sb) - Y[i]/(pa+x[i]*pb);
      chi2a += chi*chi;
      nq++;
    } else if(x[i] >= 4900.0 && x[i] <= 5400.0) {
      if(flag3 == 0) {
        lstsqr(X,Y,n,4900.0,5400.0,&pa,&pb);
        lstsqr(x,y,n,4900.0,5400.0,&sa,&sb);
        flag3 = 1;
      }
      chi = y[i]/(sa+x[i]*sb) - Y[i]/(pa+x[i]*pb);
      // chi2a += chi*chi;  removed because was skewing metallic-line type.
      // nq++;
    }
  }
  chi2a /= (float)nq;


  // added 5/26/15
  chi = CHband2(spcode,Lumcode);

  fprintf(Log,"CH band dif = %f\n",chi);

  chi2 = chi2a + 3.0*chi*chi;

  /*
  chi2b = 0.0;
  flag1 = flag2 = flag3 = 0;
  nq = 0;
  sumy = sumY = 0.0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4150.0 && x[i] <= 4280.0) {
      if(flag1 == 0) {
        lstsqr(X,Y,n,4150.0,4280.0,&pa,&pb);
        lstsqr(x,y,n,4150.0,4280.0,&sa,&sb);
        flag1 = 1;
      }
      chi = y[i]/(sa+x[i]*sb) - Y[i]/(pa+x[i]*pb);
      chi2b += chi*chi;
      chi2b1 = chi2b;
      nq++;
    }
    else if(x[i] >= 4410.0 && x[i] <= 4600.0) {
      if(flag2 == 0) {
        lstsqr(X,Y,n,4410.0,4600.0,&pa,&pb);
        lstsqr(x,y,n,4410.0,4600.0,&sa,&sb);
        flag2 = 1;
      }
      chi = y[i]/(sa+x[i]*sb) - Y[i]/(pa+x[i]*pb);
      chi2b += chi*chi;
      chi2b2 = chi2b;
      nq++;
    } else if(x[i] >= 4900.0 && x[i] <= 5400.0) {
      if(flag3 == 0) {
        lstsqr(X,Y,n,4900.0,5400.0,&pa,&pb);
        lstsqr(x,y,n,4900.0,5400.0,&sa,&sb);
        flag3 = 1;
      }
      chi = y[i]/(sa+x[i]*sb) - Y[i]/(pa+x[i]*pb);
      // chi2b += chi*chi;
      // nq++;
    }
  }
  chi2b /= (float)nq;

  /* Note: the below ensures that the change in criteria between F0 and
     F5 does not result in a discontinuity in chi2 */

  /*
  if(spcode >= 26.0) chi2 = chi2a;
  else if(spcode <= 23.0) chi2 = chi2b;
  else chi2 = chi2b + (chi2a-chi2b)*(spcode-23.0)/(3.0);
*/

  // chi2 = chi2a;

  free_vector(x,0,N);
  free_vector(y,0,N);

  // printf("spcode = %f  chi2 = %e\n",spcode,chi2);

  return(chi2);
}
