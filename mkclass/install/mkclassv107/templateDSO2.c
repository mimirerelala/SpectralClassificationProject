// Applies a flux template (derived from a library spectrum) to a program spectrum to help deal with poorly flux-calibrated spectra
#include <stdio.h>
#include <stdlib.h>
double fflux();

void templateDSO(float *xin,float *yin,float *xt,float *yt,int kin)
{
  double *flxin,*flxt,*ratio,rat,flxmax = -0.01;
  int i,j,kt,k,nmax;
  double band[9] = {3875.0,4020.0,4211.0,4500.0,4570.0,4805.0,4940.0,
                    5100.0,5450.0};
  extern int N;

  flxin = (double *)calloc(9,sizeof(double));
  flxt = (double *)calloc(9,sizeof(double));
  ratio = (double *)calloc(9,sizeof(double));

  nmax = 9;

  for(i=0;i<9;i++) {
    if(band[i] + 10.0 > xin[kin-1]) {
      nmax = i;
      break;
    }
  }

  for(i=0;i<nmax;i++) {
    flxin[i] = fflux(xin,yin,band[i],kin);
    flxt[i] = fflux(xt,yt,band[i],N);
    ratio[i] = flxt[i]/flxin[i];
  }

  for(i=0;i<kin;i++) {
    if(xin[i] < band[0])
      rat = ratio[0] + (ratio[1]-ratio[0])*(xin[i]-band[0])/(band[1]-band[0]);
    else if(xin[i] >= band[nmax-1])
      rat = ratio[nmax-2] + (ratio[nmax-1]-ratio[nmax-2])*
        (xin[i]-band[nmax-2])/(band[nmax-1]-band[nmax-2]);
    else{
      k = 0;
      for(j=0;j<nmax;j++) {
	if(xin[i] >= band[j] && xin[i] < band[j+1]) k = j;
      }
      rat = ratio[k] + (ratio[k+1]-ratio[k])*(xin[i]-band[k])/
            (band[k+1]-band[k]);
    }
    yin[i] = yin[i]*rat;
  }

  for(i=0;i<kin;i++) {
    if(xin[i] >= 4502.0 && xin[i] <= 4508.0) {
      if(yin[i] > flxmax) flxmax = yin[i];
    }
  }

  for(i=0;i<kin;i++) yin[i] /= flxmax;
}


double fflux(x,y,band,kin)
float *x,*y;
double band;
int kin;
{
  static double wave[21] = {-10.0,-9.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,
                       -1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
  static double trans[21] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.0, 1.0, 1.0,
			1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.8,0.6,0.4,0.2, 0.0};

  double wav[21];
  int i,j,k;
  double Y1,y2,ft,A=0.0,integral=0.0,a1,a2;
  k = 1;

  Y1 = 0.0;
  a1 = 0.0;

  for(i=0;i<21;i++) wav[i] = wave[i] + band;

  j = 0;
  while(x[j] < wav[20] && j < kin) {
    if(x[j] <= wav[0]) {
      j++;
      continue;
    }
    if(x[j] >= wav[20]) break;
    for(i=k;i<=20;i++) {
      if(x[j] > wav[i-1] && x[j] <= wav[i]) {
	k = i;
	break;
      }
    } 
    ft = trans[k-1] + (trans[k]-trans[k-1])*(x[j]-wav[k-1])/(wav[k]-wav[k-1]);
    y2 = ft*y[j];
    a2 = ft;
    integral += 0.5*(Y1+y2)*(x[j]-x[j-1]);
    A += 0.5*(a1+a2)*(x[j]-x[j-1]);
    a1 = a2;
    Y1 = y2;
    j++;
  }
  return(integral/A);
}
