// This file contains functions for identifying carbon stars, white dwarfs, etc.
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "util.h"

void hydrat(float *wave, float *flx, int N, float *width, float *ratio);
void getspectrum();
void rebin();
float lateM(float *wave, float *flx, int N);
float DOB(float *wave, float *flx, int N);
float DB(float *wave, float *flx, int N);
float DZ(float *wave, float *flx, int N, float *Cratio);
float DO(float *wave, float *flx, int N);
float Carbon(float *wave, float *flx, int N);
int emission(float *wave, float *flx, int N);
int lowSN(float *wave, float *flx, int N);
extern FILE *Log;

int DetectNN2(float *x, float *y, int n)
{
  int k,e;
  float width,ratio,index;
  extern float wlow;
  extern float whigh;
  extern int N;
  extern float space;
  extern FILE *Log;
  int l1,l2;
  float Cratio;
  int SN = 0;

  /* Will return:
     0   normal star
     1   DA
     2   DB
     3   DO
     4   DZ
     5   DQ
     6   unknown
     7   M-type star
     8   Carbon star
     9   emission-line (type undetermined)
    10   WN
    11   Helium nova
    12   WC
    13   Low S/N or faulty reduction
    14   Unclassifiable
  */

  /* Check to see if the star is low S/N.  This is currently done simply
by checking to see if some of the flux points are less than zero.  This may
also mean a bad cosmic ray or a zeropoint problem */

  SN = lowSN(x,y,N);
  if(SN == 1) return(13);
  if(SN == 2) return(14); 

  /* Check to see if star is an emission-line star */


  e = emission(x,y,n);
  fprintf(Log,"emission e = %d\n",e);

  if(e == 1) return(9);
  if(e == 2) return(10);
  if(e == 3) return(11);
  if(e == 4) return(12);

  /* Determine Hydrogen-line index based on ratio of width to depth */

  if(whigh < 4700.0) {
    hydrat(x,y,n,&width,&ratio);
    fprintf(Log,"NN3: ratio = %f width = %f\n",ratio,width);
    if(ratio > 0 && ratio < 65 && width > 2 && width < 27.5) return(0);
    // The next stars are in a questionable region
    if(ratio >= 60 && ratio <= 125 && width >= 13.0 && width <= 42) return(13);
    // These are probably DAs
    if(width > 26.0 && ratio > 125) return(1);
    if(width < 0.0 || ratio < 0.0) return(6);
    else return(6);
  } else {
    ratio = Carbon(x,y,n);
    fprintf(Log,"Carbon ratio = %f\n",ratio);
    if(ratio > 3.0) return(8);
    hydrat(x,y,n,&width,&ratio);
    fprintf(Log,"NN3: ratio = %f width = %f\n",ratio,width);
  /* Check for Normal stars */
    // printf("hydrat: ratio = %f  width = %f\n",ratio,width);
    if(ratio >= 0 && ratio <= 65 && width >= 1.9 && width <= 27.5) return(0);
  /* Check to see if late M-type star */
    index = lateM(x,y,n);
    if(index > 1.0 && ratio < 20 && ratio > 0 || index > 2.0) return(0);
  /* Check to see if DA */
    // The next guys may have defective spectra
    if(ratio >= 60 && ratio <= 125 && width >= 13.0 && width <= 42) return(13);
    // These are probably DAs
    if(width > 26.0 && ratio > 125) {
      fprintf(Log,"NN3: ratio = %f width = %f\n",ratio,width);
      return(1);
    }
    if(width < 0.0 || ratio < 0.0) return(6);
  /* Check for other types of WD and Carbon stars */
    index = DB(x,y,n);
    fprintf(Log,"DB = %f\n",index);
    if(index > 1.20) return(2);  // changed from 1.07
    index = DZ(x,y,n,&Cratio);
    fprintf(Log,"DZ index = %f Cratio = %f\n",index,Cratio);
    // printf("DZ index = %f Cratio = %f\n",index,Cratio);
    /* This distinguishes DZ stars and CN-strong stars that were not caught
       with the Carbon function above */
    if(index > 1.05 && Cratio < 1.3) return(4);
    else return(0);
    index = DO(x,y,n);
    if(index > 1.02) return(3);
  }
  return(6);
}

int lowSN(float *wave, float *flx, int N)
{
  int count = 0;
  int i;

  // Are more than 3 points negative?

  for(i=0;i<N;i++) {
    if(flx[i] < 0.0) count++;
  }

  // printf("count = %d\n",count);

  // if more than 50 points negative, declare it unclassifiable
  if(count > 3 && count < 50) return(1);
  else if(count >= 50) return(2);
  else return(0);
}

void hydrat(float *wave, float *flx, int N, float *width, float *ratio)
{
  float cont1,cont2,wid1,wid2,depth,depth2;
  float sum1,sum2,sumH,slope,bot,b,height,diff,min,y;
  int i,j,k,m,n;

  sum1 = sum2 = 0.0;
  bot = 1.0e+30;
  m = n = 0;
  for(i=0;i<N;i++) {
    if(wave[i] >= 4190.0 && wave[i] <= 4230.0) {
       sum1 += flx[i];
       m++;
    }
    if(wave[i] >= 4450.0 && wave[i] <= 4500.0) {
       sum2 += flx[i];
       n++;
    }
    if(wave[i] >= 4330.0 && wave[i] <= 4350.0 && bot > flx[i]) bot = flx[i];
  }
  sum1 /= m;
  sum2 /= n;
  slope = (sum2-sum1)/265.0;
  b = sum1 - slope*4210.0;
  height = slope*4340.0 + b;
  depth = height - bot;
  depth2 = depth/2;
  depth /= height;
  // printf("sum1 = %f sum2 = %f b = %f\n",sum1,sum2,b);
  fprintf(Log,"NN3: slope = %f  height = %f depth = %f\n",slope, height,depth);

  /* Find wavelength point closest to center of H gamma */
  min = 10;
  for(i=0;i<N;i++) {
    diff = fabs(wave[i]-4340.4);
    if(diff <= min) {
      min = diff;
      j = i;
    }
  }

  /* Find blue wing midpoint */

  k = j;
  y = slope*wave[k] + (b-depth2) - flx[k];

  while(y >= 0) {
    y = slope*wave[k] + (b-depth2) - flx[k];
    k--;
  }
  // printf("k = %d\n",k);
  wid1 = wave[k+1];

  /* Find red wing midpoint */

  k = j;
  y = slope*wave[k] + (b-depth2) - flx[k];

  while(y >= 0) {
    y = slope*wave[k] + (b-depth2) - flx[k];
    k++;
  }
  wid2 = wave[k-1];
  // printf("k = %d\n",k);
  fprintf(Log,"wid1 = %f  wid2 = %f\n",wid1,wid2);
  *width = wid2-wid1;
  *ratio = *width/depth;
}

float lateM(float *wave, float *flx, int N)
{
  float sum1,sum2;
  int i;

  sum1 = sum2 = 0.0;

  for(i=0;i<N;i++) {
    if(wave[i] >= 4918.0 && wave[i] <= 4948.0) sum1 += flx[i];
    if(wave[i] >= 4958.0 && wave[i] <= 4988.0) sum2 += flx[i];
  }
  return(sum1/sum2);
}

float Carbon(float *wave, float *flx, int N)
{
  float sum1,sum2,sum3,sum4;
  int i;
  sum1 = sum2 = sum3 = sum4 = 0.0;

  for(i=0;i<N;i++) {
    if(wave[i] >= 4675.0 && wave[i] <= 4725.0) sum1 += flx[i];
    if(wave[i] >= 4755.0 && wave[i] <= 4805.0) sum2 += flx[i];
    if(wave[i] >= 5070.0 && wave[i] <= 5140.0) sum3 += flx[i];
    if(wave[i] >= 5176.0 && wave[i] <= 5246.0) sum4 += flx[i];
  }

  return(sum2/sum1 + sum4/sum3);
}

float DOB(float *wave, float *flx, int N)
{
  float sum1,sum2;
  int i;

  sum1 = sum2 = 0.0;

  for(i=0;i<N;i++) {
    if(wave[i] >= 4815.0 && wave[i] <= 4905.0) sum1 += flx[i];
    if(wave[i] >= 4905.0 && wave[i] <= 4960.0) sum2 += flx[i];
  }

  return((sum2*90)/(sum1*55));
}

float DO(float *wave, float *flx, int N)
{
  float cont,line,ratio;
  int i;

  cont = line = 0.0;

  for(i=0;i<N;i++) {
    if(wave[i] >= 4780.0 && wave[i] <= 4820.0) cont += flx[i];
    if(wave[i] >= 4820.0 && wave[i] <= 4900.0) line += flx[i];
    if(wave[i] >= 4900.0 && wave[i] <= 4940.0) cont += flx[i];
  }
  ratio = cont/line;
  return(ratio);
}

float DB(float *wave, float *flx, int N)
{
  float cont1,sum,cont2,cont,line,ratio;
  int i;

  cont1 = sum = cont2 = 0.0;

  for(i=0;i<N;i++) {
    if(wave[i] >= 4301.0 && wave[i] <= 4351.0) cont1 += flx[i];
    if(wave[i] >= 4421.0 && wave[i] <= 4521.0) sum += flx[i];
    if(wave[i] >= 4591.0 && wave[i] <= 4641.0) cont2 += flx[i];
  }
  cont = (cont1+cont2)/100.0;
  line = sum/100.0;
  ratio = cont/line;

  return(ratio);
}

float DZ(float *wave, float *flx, int N, float *Cratio)
{
  float cont,line,ratio;
  int i;
  float ratio1,ratio2,sum1,sum2,sum3,sum4,dif;

  cont = line = 0.0;

  for(i=0;i<N;i++) {
    if(wave[i] > 3850.0 && wave[i] <= 3900.0) cont += flx[i];
    if(wave[i] > 3925.0 && wave[i] <= 3940.0) line += flx[i];
    if(wave[i] > 3960.0 && wave[i] <= 3975.0) line += flx[i];
    if(wave[i] > 4000.0 && wave[i] <= 4050.0) cont += flx[i];
  }

  ratio = (30.0*cont)/(line*100);

  /* Check to see if there is a CN band -- that will distinguish DZ
     stars from carbon-enhanced stars */

  sum1 = sum2 = sum3 = sum4 = 0.0;
   for(i=0;i<N;i++) {
     if(wave[i] >= 4043.0 && wave[i] <= 4088.0) sum1 += flx[i];
     if(wave[i] >= 4140.0 && wave[i] <= 4210.0) sum2 += flx[i];
     if(wave[i] >= 4219.0 && wave[i] <= 4264.0) sum1 += flx[i];
   }
   if(sum2 != 0) *Cratio = sum1/sum2;
   else *Cratio = -9.99;

   /*
   printf("DZ: ratio = %f\n",ratio);
   printf("DZ: CN ratio = %f\n",ratio1);
   */

  return(ratio);
}

/* Detects most emission-line stars.  Tries to identify the type */
int emission(float *wave, float *flx, int N)
{
  float ave,adev,sdev,var,skew,curt,flxmax,wavmax;
  float ave4770,count4770,flxhalf;
  int i,j,k;
  float *Flx;
  int e = 0;
  extern FILE *Log;

  Flx = vector(1,10000);

  for(i=1;i<N-1;i++) Flx[i] = flx[i];


  flxmax = wavmax = 0.0;
  for(i=100;i<=N-20;i++) {
    if(flxmax < flx[i]) {
      flxmax = flx[i];
      wavmax = wave[i];
    }
  }
  fprintf(Log,"wavemax = %7.2f\n",wavmax);
  for(i=1;i<=N-2;i++) Flx[i] /= flxmax;

  moment(Flx,N-2,&ave,&adev,&sdev,&var,&skew,&curt);

  fprintf(Log,"Emission routine: ave = %f sdev = %f\n",ave,sdev);

  e = 0;
  //  printf("ave = %f  sdev = %f\n",ave,sdev);
  if((ave <= 0.27 && sdev < 0.185185185*ave + 0.1) || 
     (ave > 0.27 && sdev < -0.1*ave +0.130)) e = 1;
  // constant 0.165 replaced with 0.130 in above inequality

  /* Are there He II emission lines, or Carbon emission lines? */
  if(e == 1) {
    if(wavmax > 4680.0 && wavmax <= 4690.0) e = 2;
    if(wavmax > 4655.0 && wavmax <= 4665.0) e = 4;
  }

  /* If there are He II emission lines, is it a nova or a WN? */
  ave4770 = count4770 = 0.0;
  if(e == 2) {
    for(i=1;i<=N-2;i++) {
      if(wave[i] >= 4770.0 && wave[i] <= 4780.0) {
	ave4770 += Flx[i];
	count4770 += 1.0;
      }
    }
    ave4770 /= count4770;
    flxhalf = (1.0 + ave4770)/2.0;
    for(i=1;i<=N-2;i++) {
      if(wave[i] >= 4840.0 && wave[i] <= 4870.0 && Flx[i] > flxhalf) e = 3;
    }
  }

  /*
   e = 0 absorption-line spectrum
   e = 1 unidentified emission-line spectrum
   e = 2 WN
   e = 3 Helium nova
   e = 4 WC
  */

  free_vector(Flx,1,10000);

  return(e);
}
