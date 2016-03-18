// Functions to identify various common chemical peculiarities
#include <stdio.h>
#include <math.h>
#include "util.h"
void sp2class(float spcode, float lumcode, float metcode, 
	      float *x, float *y, int *n);

int peculiarity(float spcode, float lumcode, int *he, int *sr, int *si, 
                 int *eu, int *cr, int *ba, int *ch, int *cn )
{
  int flag = 0;
  extern float *X,*Y;
  extern int N;
  float *x,*y;
  float sum1,sum2,ratio1,ratio2,dif;
  int i,n;
  extern FILE *Log;

  x = vector(0,N);
  y = vector(0,N);

  *he = *sr = *si = *eu = *cr = *ba = *ch = *cn = 0;

  sp2class(spcode,lumcode,0.0,x,y,&n);

  /* Strontium Peculiarity: from A0 to K0 */
  if(spcode >= 16.0 && spcode <= 34.0) {
    sum1 = sum2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4068.0 && x[i] <= 4074.5) sum2 += y[i];
      if(x[i] >= 4074.5 && x[i] <= 4081.0) sum1 += y[i];
      if(x[i] >= 4081.0 && x[i] <= 4087.5) sum2 += y[i];
    }
    ratio1 = sum1/(0.5*sum2);
    sum1 = sum2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4068.0 && x[i] <= 4074.5) sum2 += Y[i];
      if(x[i] >= 4074.5 && x[i] <= 4081.0) sum1 += Y[i];
      if(x[i] >= 4081.0 && x[i] <= 4087.5) sum2 += Y[i];
    }
    ratio2 = sum1/(0.5*sum2);
    dif = (ratio1-ratio2)/ratio1;
    fprintf(Log,"Sr II: dif = %f\n",dif);
    if(dif >= 0.015 && dif <= 0.024) {
      *sr = 1;
      flag = 1;
    }
    if(dif > 0.024) {
      *sr = 2;
      flag = 1;
    }
  }

  /* Europium Peculiarity: from A0 to F5 */
  if(spcode >= 16.0 && spcode <= 26.0) {
    sum1 = sum2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4199.0 && x[i] <= 4203.0) sum2 += y[i];
      if(x[i] >= 4203.0 && x[i] <= 4207.0) sum1 += y[i];
      if(x[i] >= 4207.0 && x[i] <= 4211.0) sum2 += y[i];
    }
    ratio1 = sum1/(0.5*sum2);
    sum1 = sum2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4199.0 && x[i] <= 4203.0) sum2 += Y[i];
      if(x[i] >= 4203.0 && x[i] <= 4207.0) sum1 += Y[i];
      if(x[i] >= 4207.0 && x[i] <= 4211.0) sum2 += Y[i];
    }
    ratio2 = sum1/(0.5*sum2);
    dif = (ratio1-ratio2)/ratio1;
    // printf("Eu II: dif = %f\n",dif);
    if(dif >= 0.015) {
      *eu = 1;
      flag = 1;
    }
  }

  /* Silicon Peculiarity: from B5 to F0 */
  if(spcode >= 12.0 && spcode <= 23.0) {
    sum1 = sum2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4118.0 && x[i] <= 4126.0) sum2 += y[i];
      if(x[i] >= 4126.0 && x[i] <= 4134.0) sum1 += y[i];
      if(x[i] >= 4134.0 && x[i] <= 4142.0) sum2 += y[i];
    }
    ratio1 = sum1/(0.5*sum2);
    sum1 = sum2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4118.0 && x[i] <= 4126.0) sum2 += Y[i];
      if(x[i] >= 4126.0 && x[i] <= 4134.0) sum1 += Y[i];
      if(x[i] >= 4134.0 && x[i] <= 4142.0) sum2 += Y[i];
    }
    ratio2 = sum1/(0.5*sum2);
    dif = (ratio1-ratio2)/ratio1;
    fprintf(Log,"Si II: dif = %f\n",dif);
    if(dif >= 0.015) {
      *si = 1;
      flag = 1;
    }
  }

  free_vector(x,0,N);
  free_vector(y,0,N);
  return(flag);
}


int barium(float spcode, float lumcode)
{
  float *x,*y;
  extern float *X,*Y;
  extern int N;
  extern FILE *Log;
  float ratio1,ratio2,sum1,sum2,sum3,sum4,dif;
  int i,n;
  int ba = 0;

  /* In F and G-type stars with a Barium enhancement, the Ba II line
doesn't look strong in classification spectra, but Sr II 4077, especially,
does look strong.  So, for stars earlier than G5, check to see if Sr is 
enhanced.  For later than G5, check the actual Ba line. */

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,lumcode,0.0,x,y,&n);

  if(spcode <= 32.0) {
    sum1 = sum2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4068.0 && x[i] <= 4074.5) sum2 += y[i];
      if(x[i] >= 4074.5 && x[i] <= 4081.0) sum1 += y[i];
      if(x[i] >= 4081.0 && x[i] <= 4087.5) sum2 += y[i];
    }
    ratio1 = sum1/(0.5*sum2);
    sum1 = sum2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4068.0 && x[i] <= 4074.5) sum2 += Y[i];
      if(x[i] >= 4074.5 && x[i] <= 4081.0) sum1 += Y[i];
      if(x[i] >= 4081.0 && x[i] <= 4087.5) sum2 += Y[i];
    }
    ratio2 = sum1/(0.5*sum2);
    dif = (ratio1-ratio2)/ratio1;
    fprintf(Log,"Ba routine: Sr II: dif = %f\n",dif);
    if(dif >= 0.05) ba = 1;
  } else {
    sum1 = sum2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4520.0 && x[i] <= 4539.0) sum2 += y[i];
      if(x[i] >= 4546.0 && x[i] <= 4560.0) sum1 += y[i];
    }
    ratio1 = sum1/sum2;
    sum1 = sum2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4520.0 && x[i] <= 4539.0) sum2 += Y[i];
      if(x[i] >= 4546.0 && x[i] <= 4560.0) sum1 += Y[i];
    }
    ratio2 = sum1/sum2;
    dif = (ratio1-ratio2)/ratio1;
    fprintf(Log,"Ba routine: Ba index: %f\n",dif);
    if(dif >= 0.05) ba = 1;
  }

  free_vector(x,0,N);
  free_vector(y,0,N);

  return(ba);
}

/* This function measures a Mg II index based on the Mg II 4481 line, and
   can be used to indicate a possible Lambda Boo star */

float MgII(float spcode, float lumcode)
{
  float *x,*y;
  extern float *X,*Y;
  extern int N;
  float ratio1,ratio2,sum1,sum2,sum3,sum4,dif;
  int i,n;
  float mg = 0;

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,lumcode,0.0,x,y,&n);

  if(spcode <= 34.0) {
    sum1 = sum2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4478.0 && x[i] <= 4484.0) sum1 += y[i];
      if(x[i] >= 4484.0 && x[i] <= 4490.0) sum2 += y[i];
    }
    ratio1 = sum1/sum2;
    sum1 = sum2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4478.0 && x[i] <= 4484.0) sum1 += Y[i];
      if(x[i] >= 4484.0 && x[i] <= 4490.0) sum2 += Y[i];
    }
    ratio2 = sum1/sum2;
    mg = (ratio1-ratio2)/ratio1;
    return(mg);
  }

  free_vector(x,0,N);
  free_vector(y,0,N);

  return(mg);
}

/* See if the C2 4737 band is enhanced */
float carbon4737(float spcode, float lumcode)
{
  float *x,*y;
  extern float *X,*Y;
  extern int N;
  float ratio1,ratio2,sum1,sum2,sum3,sum4,dif;
  int i,n;
  float C = 0;

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,lumcode,0.0,x,y,&n);

  sum1 = sum2 = sum3 = sum4 = 0.0;
   for(i=0;i<n;i++) {
     if(x[i] >= 4585.0 && x[i] <= 4630.0) sum1 += y[i];
     if(x[i] >= 4630.0 && x[i] <= 4720.0) sum2 += y[i];
     if(x[i] >= 4720.0 && x[i] <= 4765.0) sum1 += y[i];
   }
   if(sum2 != 0) ratio1 = sum1/sum2;
   else {
     ratio1 = 1;
     C = -1;
   }
   sum1 = sum2 = 0.0;
   for(i=0;i<n;i++) {
     if(x[i] >= 4585.0 && x[i] <= 4630.0) sum1 += Y[i];
     if(x[i] >= 4630.0 && x[i] <= 4720.0) sum2 += Y[i];
     if(x[i] >= 4720.0 && x[i] <= 4765.0) sum1 += Y[i]; 
   }
   if(sum2 != 0) ratio2 = sum1/sum2;
   else {
     ratio2 = 1;
     C = -1;
   }

   if(C != -1) {
     dif = fabs(ratio1-ratio2)/ratio2;
     /* printf("Carbon: diff = %f\n",dif); */
   }
   
  free_vector(x,0,N);
  free_vector(y,0,N);

  return(dif);
}

/* See if the CN 4215 band is enhanced */
float CN4215(float spcode, float lumcode)
{
  float *x,*y;
  extern float *X,*Y;
  extern int N;
  float ratio1,ratio2,sum1,sum2,sum3,sum4,dif;
  int i,n;
  float C = 0;

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,lumcode,0.0,x,y,&n);

  sum1 = sum2 = sum3 = sum4 = 0.0;
   for(i=0;i<n;i++) {
     if(x[i] >= 4043.0 && x[i] <= 4088.0) sum1 += y[i];
     if(x[i] >= 4140.0 && x[i] <= 4210.0) sum2 += y[i];
     if(x[i] >= 4219.0 && x[i] <= 4264.0) sum1 += y[i];
   }
   if(sum2 != 0) ratio1 = sum1/sum2;
   else {
     ratio1 = 1;
     C = -1;
   }
   sum1 = sum2 = 0.0;
   for(i=0;i<n;i++) {
     if(x[i] >= 4043.0 && x[i] <= 4088.0) sum1 += Y[i];
     if(x[i] >= 4140.0 && x[i] <= 4210.0) sum2 += Y[i];
     if(x[i] >= 4219.0 && x[i] <= 4264.0) sum1 += Y[i]; 
   }
   if(sum2 != 0) ratio2 = sum1/sum2;
   else {
     ratio2 = 1;
     C = -1;
   }

   if(C != -1) {
     dif = fabs(ratio1-ratio2)/ratio2;
     /* printf("CN: diff = %f\n",dif); */
   }
   
  free_vector(x,0,N);
  free_vector(y,0,N);

  return(dif);
}

/* See if the G-band is enhanced or weak */
float CHband(float spcode, float lumcode)
{
  float *x,*y;
  extern float *X,*Y;
  extern int N;
  float ratio1,ratio2,sum1,sum2,sum3,sum4,dif;
  int i,n;
  float C = 0;

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,lumcode,0.0,x,y,&n);

  sum1 = sum2 = sum3 = sum4 = 0.0;
   for(i=0;i<n;i++) {
     if(x[i] >= 4243.0 && x[i] <= 4265.0) sum1 += y[i];
     if(x[i] >= 4294.0 && x[i] <= 4314.0) sum2 += y[i];
     if(x[i] >= 4347.0 && x[i] <= 4377.0) sum1 += y[i];
   }
   if(sum2 != 0) ratio1 = sum1/sum2;
   else {
     ratio1 = 1;
     C = -1;
   }
   sum1 = sum2 = 0.0;
   for(i=0;i<n;i++) {
     if(x[i] >= 4243.0 && x[i] <= 4265.0) sum1 += Y[i];
     if(x[i] >= 4294.0 && x[i] <= 4314.0) sum2 += Y[i];
     if(x[i] >= 4347.0 && x[i] <= 4377.0) sum1 += Y[i]; 
   }
   if(sum2 != 0) ratio2 = sum1/sum2;
   else {
     ratio2 = 1;
     C = -1;
   }

   if(C != -1) {
     dif = fabs(ratio1-ratio2)/ratio2;
     /* printf("CH: diff = %f\n",dif); */
   }
   
  free_vector(x,0,N);
  free_vector(y,0,N);

  return(dif);
}

float CHband2(float spcode, float lumcode)
{
  float *x,*y;
  extern float *X,*Y;
  extern int N;
  float ratio1,ratio2,sum1,sum2,sum3,sum4,dif;
  int i,n;
  float C = 0;
  float G1,G2;
  float cont1,cont2,cont,pos1,pos2,ratio;
  float flx;

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,lumcode,0.0,x,y,&n);

  /* first find continuum points */
  cont1 = cont2 = pos1 = pos2 = 0.0;

  for(i=0;i<n;i++) {
    if(x[i] >= 4261 && x[i] <= 4268) {
      if(y[i] > cont1) {
	cont1 = y[i];
	pos1 = x[i];
      }
    }
    if(x[i] >= 4316 && x[i] <= 4320) {
      if(y[i] > cont2) {
	cont2 = y[i];
	pos2 = x[i];
      }
    }
  }
  /* Now calculate G-band index */
  G1 = 0.0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4298.0 && x[i] <= 4309.0) {
      cont = cont1 + (cont2-cont1)*(x[i]-pos1)/(pos2-pos1);
      G1 += (1.0 - y[i]/cont);
    }
  }

  /* first find continuum points */
  cont1 = cont2 = pos1 = pos2 = 0.0;

  for(i=0;i<n;i++) {
    if(X[i] >= 4261 && X[i] <= 4268) {
      if(Y[i] > cont1) {
	cont1 = Y[i];
	pos1 = X[i];
      }
    }
    if(X[i] >= 4316 && X[i] <= 4320) {
      if(Y[i] > cont2) {
	cont2 = Y[i];
	pos2 = X[i];
      }
    }
  }
  /* Now calculate G-band index */
  G2 = 0.0;
  for(i=0;i<n;i++) {
    if(X[i] >= 4298.0 && X[i] <= 4309.0) {
      cont = cont1 + (cont2-cont1)*(X[i]-pos1)/(pos2-pos1);
      G2 += (1.0 - Y[i]/cont);
    }
  }

  /* printf("G1 = %f\n",G1);
     printf("G2 = %f\n",G2); */

  ratio = (G2-G1)/G1;
   
  free_vector(x,0,N);
  free_vector(y,0,N);

  return(ratio);
}

