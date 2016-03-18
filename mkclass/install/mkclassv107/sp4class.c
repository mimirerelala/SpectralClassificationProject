// Interpolation function for pulling spectra from spectral libraries
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "util.h"
#define I 37
#define J 3
#define K 1
FILE *ffopen();

void sp2class(float spcode, float lumcode, float metcode, 
	      float *x, float *y, int *n)
{
  double t[I] = {3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.0,13.0,14.0,15.0,16.0,17.0,19.0,20.0,21.0,23.0,24.0,25.0,
		 26.0,27.0,28.0,29.0,30.0,31.0,32.0,33.0,34.0,36.0,37.0,39.0,40.0,40.7,42.5,44.0,45.5};
  double l[J] = {1.0,3.0,5.0};
  double m[K] = {0.0};
  float *Y1,*y2,*y3,*y4,q,T,U;
  int i,k,k1,k2,l1,l2,nq;
  char name1[80],name2[80],name3[80],name4[80];
  extern char lib[40];
  extern int N;
  extern char MKLIB[300];
  FILE *in1,*in2,*in3,*in4;

  Y1 = vector(0,N);
  y2 = vector(0,N);
  y3 = vector(0,N);
  y4 = vector(0,N);

  for(i=0;i<I-1;i++) {
    if(spcode < t[0]) {
      k1 = 0;
      k2 = 1;
      break;
    }
    if(spcode > t[I-1]) {
      k1 = I-2;
      k2 = I-1;
      break;
    }
    if(spcode >= t[i] && spcode <= t[i+1]) {
      k1 = i;
      k2 = i+1;
      break;
    }
  }


  /* Allows extrapolation beyond limits of luminosity classes V -> Ib */
  if(lumcode < l[0]) {
    l1 = 0;
    l2 = 1;
  } else if(lumcode > l[J-1]) {
    l1 = J-2;
    l2 = J-1;
  } else {
    for(i=0;i<J-1;i++) {
      if(lumcode >= l[i] && lumcode <= l[i+1]) {
        l1 = i;
        l2 = i+1;
        break;
      }
    }
  }
  T = (spcode - t[k1])/(t[k2]-t[k1]);
  U = (lumcode - l[l1])/(l[l2]-l[l1]);

  sprintf(name1,"%s/%s/t%03dl%2dp00.rbn",
          MKLIB,lib,(int)(10.0*t[k1]),(int)(10.0*l[l1]));
  sprintf(name2,"%s/%s/t%03dl%2dp00.rbn",
          MKLIB,lib,(int)(10.0*t[k2]),(int)(10.0*l[l1]));
  sprintf(name3,"%s/%s/t%03dl%2dp00.rbn",
          MKLIB,lib,(int)(10.0*t[k2]),(int)(10.0*l[l2]));
  sprintf(name4,"%s/%s/t%03dl%2dp00.rbn",
          MKLIB,lib,(int)(10.0*t[k1]),(int)(10.0*l[l2]));

  in1 = ffopen(name1,"r");
  in2 = ffopen(name2,"r");
  in3 = ffopen(name3,"r");
  in4 = ffopen(name4,"r");

  i = 0;
  while(fscanf(in1,"%f %f",&x[i],&Y1[i]) != EOF) i++;
  i=0;
  while(fscanf(in2,"%f %f",&q,&y2[i]) != EOF) i++;
  i=0;
  while(fscanf(in3,"%f %f",&q,&y3[i]) != EOF) i++;
  i=0;
  while(fscanf(in4,"%f %f",&q,&y4[i]) != EOF) i++;

  nq = i;

  for(i=0;i<nq;i++) {
    y[i] = (1.0-T)*(1.0-U)*Y1[i] + T*(1-U)*y2[i] + T*U*y3[i] + (1.0-T)*U*y4[i];
  }
  *n = nq;
  free_vector(Y1,0,N);
  free_vector(y2,0,N);
  free_vector(y3,0,N);
  free_vector(y4,0,N);
  fclose(in1);
  fclose(in2);
  fclose(in3);
  fclose(in4);
  return;
}

