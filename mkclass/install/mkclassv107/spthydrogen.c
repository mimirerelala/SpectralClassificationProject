// Various functions having to do with hydrogen and helium lines
// Some changes to sptGlines made on March 14, 2014.
#include <stdio.h>
#include <math.h>
#include "util.h"
void sp2class(float spcode, float lumcode, float metcode, 
	      float *x, float *y,int *n);
float hydrogen_profile_hot(float spcode);

float spthydrogen(float spcode)
{
  extern float *X,*Y;
  extern float Lumcode;
  extern int N;
  extern int sf;
  float *x,*y;
  float chi;
  float chi2 = 0.0;
  float lHd[5] = {4097.0,4091.0,4091.0,4097.0,4097.0};
  float hHd[5] = {4107.0,4112.0,4112.0,4107.0,4107.0};
  float lHg[5] = {4335.0,4330.0,4330.0,4335.0,4335.0};
  float hHg[5] = {4345.0,4350.0,4350.0,4345.0,4345.0};
  float lHb[5] = {4855.0,4850.0,4850.0,4855.0,4855.0};
  float hHb[5] = {4865.0,4870.0,4870.0,4865.0,4865.0};
  int i,n;

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,Lumcode,0.0,x,y,&n);

  for(i=0;i<n;i++) {
   if((x[i] >= lHd[sf] && x[i] <= hHd[sf]) || 
         (x[i] >= lHg[sf] && x[i] <= hHg[sf]) ||
         (x[i] >= lHb[sf] && x[i] <= hHb[sf]))
     {
       chi = y[i] - Y[i];
       chi2 += chi*chi;
     }
  }
  chi2 /= 20.0;

  free_vector(x,0,N);
  free_vector(y,0,N);
  return(chi2);
}


 /* This function determines the width of the H gamma line at mid-depth.
    It can be used to distinguish between a normal star and a DA white dwarf */

float hydD2()
{
  extern float *X,*Y;
  extern int Nq;
  extern float space;
  float wavemin,flxmin,halfblue,halfred;
  float cont1,cont2,cont;
  float sum1,sum2;
  float b = 0;
  float r = 0;
  float D2,bD2,rD2,S2,bw,rw,sb,sr;
  int i;
  extern FILE *Log;

  sum1 = sum2 = 0.0;
  flxmin = 1.0e+30;
  for(i=0;i<Nq;i++) {
    if(X[i] >= 4200.0 && X[i] <= 4220.0) {
      b += 1.0;
      sum1 += Y[i];
    }
    if(X[i] > 4330.0 && X[i] < 4350.0) {
      if(Y[i] <= flxmin) {
	flxmin = Y[i];
	wavemin = X[i];
      }
    }
    if(X[i] >= 4510 && X[i] <= 4570.0) {
      r += 1.0;
      sum2 += Y[i];
    }
  }
  cont1 = sum1/b;
  cont2 = sum2/r;
  fprintf(Log,"sum1 = %f sum2 = %f wavemin = %f\n",sum1,sum2,wavemin);

  cont = cont1 + (cont2 - cont1)*(wavemin - 4210)/(4540.0-4210.0);
  D2 = (cont + flxmin)/2.0;
  fprintf(Log,"D2 = %f\n",D2);

  for(i=0;i<Nq;i++) {
    if(X[i] >= 4210) {
      S2 = Y[i] - (D2 + (cont2 - cont1)*(wavemin - 4210)/(4540.0-4210.0));
      if(S2 <= 0.0) {
	bw = X[i];
	break;
      }
    }
  }
  for(i=0;i<Nq;i++) {
    if(X[i] >= wavemin) {
      S2 = Y[i] - (D2 + (cont2 - cont1)*(wavemin - 4210)/(4540.0-4210.0));
      if(S2 >= 0.0) {
	rw = X[i];
	break;
      }
    }
  }
  fprintf(Log,"rw = %f  bw = %f\n",rw,bw);
  /* Check symmetry to make certain we are actually measuring the hydrogen 
     lines */
  sb = wavemin - bw;
  sr = rw - wavemin;
  if(fabs(sb - sr) > 2.0) return(0.0);
  else {
    D2 = rw - bw;
    return(D2);
  }
  return(0.0);
}

float spt2hyd(float pcode[])
{
  extern float *X,*Y;
  extern int N;
  extern int sf;
  extern float whigh;
  float lumcode,spcode;
  float *x,*y;
  int i,n;
  float chi;
  float chi2 = 0.0;
  float ssumc1,ssumI,ssumc2,psumc1,psumI,psumc2;
  float sindx,pindx;

  x = vector(0,N);
  y = vector(0,N);

  spcode = pcode[1];
  lumcode = pcode[2];

  sp2class(spcode,lumcode,0.0,x,y,&n);
 
  ssumc1 = ssumI = ssumc2 = psumc1 = psumI = psumc2 = 0.0;

  for(i=0;i<n;i++) {
    if(x[i] >= 4053.0 && x[i] <= 4073.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4094.0 && x[i] <= 4110.0) {
      ssumI += y[i];
      psumI += Y[i];
    }
    if(x[i] >= 4131.0 && x[i] <= 4151.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  ssumc1 /= 20.0;
  psumc1 /= 20.0;
  ssumI /= 16.0;
  psumI /= 16.0;
  ssumc2 /= 20.0;
  psumc2 /= 20.0;
  sindx = ssumI/(0.5*ssumc1 + 0.5*ssumc2);
  pindx = psumI/(0.5*psumc1 + 0.5*psumc2);
  chi2 += (sindx - pindx)*(sindx-pindx);

  ssumc1 = ssumI = ssumc2 = psumc1 = psumI = psumc2 = 0.0;

  for(i=0;i<n;i++) {
    if(x[i] >= 4233.0 && x[i] <= 4248.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4332.0 && x[i] <= 4348.0) {
      ssumI += y[i];
      psumI += Y[i];
    }
    if(x[i] >= 4355.0 && x[i] <= 4378.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  ssumc1 /= 15.0;
  psumc1 /= 15.0;
  ssumI /= 16.0;
  psumI /= 16.0;
  ssumc2 /= 23.0;
  psumc2 /= 23.0;
  sindx = ssumI/(0.206*ssumc1 + 0.794*ssumc2);
  pindx = psumI/(0.206*psumc1 + 0.794*psumc2);
  chi2 += (sindx - pindx)*(sindx-pindx);

  if(whigh > 4920) {

    for(i=0;i<n;i++) {
      if(x[i] >= 4805.0 && x[i] <= 4825.0) {
        ssumc1 += y[i];
        psumc1 += Y[i];
      }
      if(x[i] >= 4853.0 && x[i] <= 4869.0) {
        ssumI += y[i];
        psumI += Y[i];
      }
      if(x[i] >= 4897.0 && x[i] <= 4917.0) {
        ssumc2 += y[i];
        psumc2 += Y[i];
      }
    }
    ssumc1 /= 20.0;
    psumc1 /= 20.0;
    ssumI /= 16.0;
    psumI /= 16.0;
    ssumc2 /= 20.0;
    psumc2 /= 20.0;
    sindx = ssumI/(0.5*ssumc1 + 0.5*ssumc2);
    pindx = psumI/(0.206*psumc1 + 0.794*psumc2);
    chi2 += (sindx - pindx)*(sindx-pindx);
  }

  free_vector(x,0,N);
  free_vector(y,0,N);
  return(chi2);
}

float hydrogen_index(float spcode)
{
  extern float *X,*Y;
  extern int N;
  extern int sf;
  extern float Lumcode;
  extern float whigh;
  float *x,*y;
  int i,n;
  float chi;
  float chi2 = 0.0;
  float ssumc1,ssumI,ssumc2,psumc1,psumI,psumc2;
  float sindx,pindx;

  if(spcode > 36.0) return(1.0e+10);

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,Lumcode,0.0,x,y,&n);
 
  ssumc1 = ssumI = ssumc2 = psumc1 = psumI = psumc2 = 0.0;

  if(spcode <= 19.0) return(1.0e+30); /* only for A stars and later */
  for(i=0;i<n;i++) {
    if(x[i] >= 4052.0 && x[i] <= 4072.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4082.0 && x[i] <= 4122.0) {
      ssumI += y[i];
      psumI += Y[i];
    }
    if(x[i] >= 4132.0 && x[i] <= 4152.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  ssumc1 /= 20.0;
  psumc1 /= 20.0;
  ssumI /= 40.0;
  psumI /= 40.0;
  ssumc2 /= 20.0;
  psumc2 /= 20.0;
  sindx = ssumI/(0.5*ssumc1 + 0.5*ssumc2);
  pindx = psumI/(0.5*psumc1 + 0.5*psumc2);
  chi2 += (sindx - pindx)*(sindx-pindx);

  ssumc1 = ssumI = ssumc2 = psumc1 = psumI = psumc2 = 0.0;

  for(i=0;i<n;i++) {
    if(x[i] >= 4233.0 && x[i] <= 4248.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4320.0 && x[i] <= 4360.0) {
      ssumI += y[i];
      psumI += Y[i];
    }
    if(x[i] >= 4355.0 && x[i] <= 4378.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  ssumc1 /= 15.0;
  psumc1 /= 15.0;
  ssumI /= 40.0;
  psumI /= 40.0;
  ssumc2 /= 23.0;
  psumc2 /= 23.0;
  sindx = ssumI/(0.206*ssumc1 + 0.794*ssumc2);
  pindx = psumI/(0.206*psumc1 + 0.794*psumc2);
  chi2 += (sindx - pindx)*(sindx-pindx);

  if(whigh > 4920) {
    for(i=0;i<n;i++) {
      if(x[i] >= 4805.0 && x[i] <= 4825.0) {
        ssumc1 += y[i];
        psumc1 += Y[i];
      }
      if(x[i] >= 4841.0 && x[i] <= 4881.0) {
        ssumI += y[i];
        psumI += Y[i];
      }
      if(x[i] >= 4897.0 && x[i] <= 4917.0) {
        ssumc2 += y[i];
        psumc2 += Y[i];
      }
    }
    ssumc1 /= 20.0;
    psumc1 /= 20.0;
    ssumI /= 40.0;
    psumI /= 40.0;
    ssumc2 /= 20.0;
    psumc2 /= 20.0;
    sindx = ssumI/(0.5*ssumc1 + 0.5*ssumc2);
    pindx = psumI/(0.206*psumc1 + 0.794*psumc2);
    chi2 += (sindx - pindx)*(sindx-pindx);
  }

  free_vector(x,0,N);
  free_vector(y,0,N);
  return(chi2);
}

/* The following code compares hydrogen-line profiles.  To be most general,
it includes a rough rectification routine for each profile, so that flux-
calibrated spectra can be accomodated */

float hydrogen_profile(float spcode)
{
  extern float *X,*Y;
  extern int N;
  extern int sf;
  extern float Lumcode;
  extern float whigh;
  float *x,*y;
  float yr,Yr,scont,pcont;
  int i,n,l,m;
  float chi;
  float chi2 = 0.0;
  float ssumc1,ssumc2,psumc1,psumc2;
  float sindx,pindx;

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,Lumcode,0.0,x,y,&n);

  /* H-delta line */

  if(spcode <= 19.0) return(1.0e+30); /* Only for A-stars and later */
  chi2 = 0.0;
  ssumc1 = ssumc2 = psumc1 = psumc2 = 0.0;
  l = m = 0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4012.0 && x[i] <= 4052.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
      l++;
    }
    if(x[i] >= 4152.0 && x[i] <= 4192.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
      m++;
    }
  }
  ssumc1 /= (double)l;
  psumc1 /= (double)l;
  ssumc2 /= (double)m;
  psumc2 /= (double)m;
  for(i=0;i<n;i++) {
    if(x[i] >= 4062.0 && x[i] <= 4142.0) {
      scont = ssumc1 + (ssumc2-ssumc1)*(x[i] - 4032.0)/(4172.0-4032.0);
      pcont = psumc1 + (psumc2-psumc1)*(x[i] - 4032.0)/(4172.0-4032.0);
      yr = y[i]/scont;
      Yr = Y[i]/pcont;
      chi2 += (yr-Yr)*(yr-Yr);
    }
  }

  /* H-gamma line */

  ssumc1 = ssumc2 = psumc1 = psumc2 = 0.0;
  l = m = 0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4250.0 && x[i] <= 4290.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
      l++;
    }
    if(x[i] >= 4390.0 && x[i] <= 4430.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
      m++;
    }
  }
  ssumc1 /= (double)l;
  psumc1 /= (double)l;
  ssumc2 /= (double)m;
  psumc2 /= (double)m;
  for(i=0;i<n;i++) {
    if(x[i] >= 4300.0 && x[i] <= 4380.0) {
      scont = ssumc1 + (ssumc2-ssumc1)*(x[i] - 4270.0)/(4410.0-4270.0);
      pcont = psumc1 + (psumc2-psumc1)*(x[i] - 4270.0)/(4410.0-4270.0);
      yr = y[i]/scont;
      Yr = Y[i]/pcont;
      chi2 += (yr-Yr)*(yr-Yr);
    }
  }

  if(whigh > 4960) {
  /* H Beta */
    ssumc1 = ssumc2 = psumc1 = psumc2 = 0.0;
    l = m = 0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4791.0 && x[i] <= 4831.0) {
        ssumc1 += y[i];
        psumc1 += Y[i];
	l++;
      }
      if(x[i] >= 4911.0 && x[i] <= 4951.0) {
        ssumc2 += y[i];
        psumc2 += Y[i];
	m++;
      }
    }
    ssumc1 /= (double)l;
    psumc1 /= (double)l;
    ssumc2 /= (double)m;
    psumc2 /= (double)m;
    for(i=0;i<n;i++) {
      if(x[i] >= 4841.0 && x[i] <= 4901.0) {
        scont = ssumc1 + (ssumc2-ssumc1)*(x[i] - 4811.0)/(4931.0-4811.0);
        pcont = psumc1 + (psumc2-psumc1)*(x[i] - 4811.0)/(4931.0-4811.0);
        yr = y[i]/scont;
        Yr = Y[i]/pcont;
        chi2 += (yr-Yr)*(yr-Yr);
      }
    }
  }

  free_vector(x,0,N);
  free_vector(y,0,N);
  return(chi2);
}

/* This function is designed to be used with B-type stars but is similar
   to the above */

float hydrogen_profile_hot(float spcode)
{
  extern float *X,*Y;
  extern int N;
  extern int sf;
  extern float Lumcode;
  extern float whigh;
  float lumcode;
  float *x,*y;
  float yr,Yr,scont,pcont;
  int i,n,l,m;
  float chi;
  float chi2 = 0.0;
  float ssumc1,ssumc2,psumc1,psumc2;
  float sindx,pindx;

  lumcode = Lumcode;

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,lumcode,0.0,x,y,&n);

  /* H-delta line */

  if(spcode >= 17.0) return(1.0e+30); /* Only for B-stars and earlier */
  chi2 = 0.0;
  ssumc1 = ssumc2 = psumc1 = psumc2 = 0.0;
  l = m = 0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4012.0 && x[i] <= 4052.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
      l++;
    }
    if(x[i] >= 4152.0 && x[i] <= 4192.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
      m++;
    }
  }
  ssumc1 /= (double)l;
  psumc1 /= (double)l;
  ssumc2 /= (double)m;
  psumc2 /= (double)m;
  for(i=0;i<n;i++) {
    if(x[i] >= 4062.0 && x[i] <= 4142.0) {
      scont = ssumc1 + (ssumc2-ssumc1)*(x[i] - 4032.0)/(4172.0-4032.0);
      pcont = psumc1 + (psumc2-psumc1)*(x[i] - 4032.0)/(4172.0-4032.0);
      yr = y[i]/scont;
      Yr = Y[i]/pcont;
      chi2 += (yr-Yr)*(yr-Yr);
    }
  }

  /* H-gamma line */

  ssumc1 = ssumc2 = psumc1 = psumc2 = 0.0;
  l = m = 0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4250.0 && x[i] <= 4290.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
      l++;
    }
    if(x[i] >= 4390.0 && x[i] <= 4430.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
      m++;
    }
  }
  ssumc1 /= (double)l;
  psumc1 /= (double)l;
  ssumc2 /= (double)m;
  psumc2 /= (double)m;
  for(i=0;i<n;i++) {
    if(x[i] >= 4300.0 && x[i] <= 4380.0) {
      scont = ssumc1 + (ssumc2-ssumc1)*(x[i] - 4270.0)/(4410.0-4270.0);
      pcont = psumc1 + (psumc2-psumc1)*(x[i] - 4270.0)/(4410.0-4270.0);
      yr = y[i]/scont;
      Yr = Y[i]/pcont;
      chi2 += (yr-Yr)*(yr-Yr);
    }
  }

  if(whigh > 4960) {
  /* H Beta */
    ssumc1 = ssumc2 = psumc1 = psumc2 = 0.0;
    l = m = 0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4791.0 && x[i] <= 4831.0) {
        ssumc1 += y[i];
        psumc1 += Y[i];
	l++;
      }
      if(x[i] >= 4911.0 && x[i] <= 4951.0) {
        ssumc2 += y[i];
        psumc2 += Y[i];
	m++;
      }
    }
    ssumc1 /= (double)l;
    psumc1 /= (double)l;
    ssumc2 /= (double)m;
    psumc2 /= (double)m;
    for(i=0;i<n;i++) {
      if(x[i] >= 4841.0 && x[i] <= 4901.0) {
        scont = ssumc1 + (ssumc2-ssumc1)*(x[i] - 4811.0)/(4931.0-4811.0);
        pcont = psumc1 + (psumc2-psumc1)*(x[i] - 4811.0)/(4931.0-4811.0);
        yr = y[i]/scont;
        Yr = Y[i]/pcont;
        chi2 += (yr-Yr)*(yr-Yr);
      }
    }
  }

  free_vector(x,0,N);
  free_vector(y,0,N);
  return(chi2);
}

/* This function checks for Helium peculiarities in B-type stars.  It
forms indices & then returns a value dif, such that if dif < 0, the
helium is weak for that spectral type, if dif ~0 it is normal, and
if dif > 0, it is strong */

float heIpec(float spcode, float lumcode)
{
  extern float *X,*Y;
  extern int N;
  extern int sf;
  extern float whigh;
  float *x,*y;
  int i,n;
  float chi;
  float chi2 = 0.0;
  float ssumc1,ssumI,ssumc2,psumc1,psumI,psumc2;
  float sindx,pindx;
  float dif;

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,lumcode,0.0,x,y,&n);
 
  ssumc1 = ssumI = ssumc2 = psumc1 = psumI = psumc2 = 0.0;

  if(spcode >= 19.0) return(1.0e+10); /* only for B stars and earlier */
  for(i=0;i<n;i++) {
    if(x[i] >= 4014.0 && x[i] <= 4019.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4019.0 && x[i] <= 4033.0) {
      ssumI += y[i];
      psumI += Y[i];
    }
    if(x[i] >= 4033.0 && x[i] <= 4038.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  ssumc1 /= 5.0;
  ssumc2 /= 5.0;
  psumc1 /= 5.0;
  psumc2 /= 5.0;
  ssumI /= 14.0;
  psumI /= 14.0;
  dif = (ssumI-psumI)/ssumI;

  return(dif);
}


/* The following code looks at typical hydrogen to iron ratios in late F and
G-type stars to determine a temperature type */
float sptGlines(float spcode)
{
  extern float *X,*Y;
  extern float Lumcode;
  extern int N;
  extern int sf;
  extern float whigh;
  float *x,*y;
  float pro1,pro2,std1,std2,ratpro,ratstd;
  float chi;
  float chi2 = 0.0;
  int i,n;

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,Lumcode,0.0,x,y,&n);

  if(spcode < 28.0) return(0.0);

  /* H delta and Fe I 4046 */
  pro1 = pro2 = std1 = std2 = 0.0;
  chi2 = 0.0;

  for(i=0;i<n;i++) {
    if(x[i] > 4038.0 && x[i] < 4050.0) {
      std1 += y[i];
      pro1 += Y[i];
    }
    if(x[i] > 4096.0 && x[i] < 4106.0) { //Narrowed March 16, 2014
      std2 += y[i];
      pro2 += Y[i];
    }
  }
  if(std2 != 0.0 && pro2 != 0.0) {
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    chi2 += pow(fabs(1.0-ratstd/ratpro),2.0);
  }

  /* H gamma and Fe I 4383 */
  pro1 = pro2 = std1 = std2 = 0.0;
  // Charged comparison from 4324 to 4383  March 16, 2014.  Also changed weight.
  for(i=0;i<n;i++) {
    if(x[i] > 4378.0 && x[i] < 4388.0) {
      std1 += y[i];
      pro1 += Y[i];
    }
    if(x[i] > 4332.0 && x[i] < 4347.0) {
      std2 += y[i];
      pro2 += Y[i];
    }
  }
  if(std2 != 0.0 && pro2 != 0.0) {
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    chi2 += 3.0*pow(fabs(1.0-ratstd/ratpro),2.0);
  }

  /* H beta and metal 4888 */
  if(whigh > 4900.0) {
    pro1 = pro2 = std1 = std2 = 0.0;

    for(i=0;i<n;i++) {
      if(x[i] > 4880.0 && x[i] < 4894.0) {
        std1 += y[i];
        pro1 += Y[i];
      }
      if(x[i] > 4851 && x[i] < 4866.0) {
        std2 += y[i];
        pro2 += Y[i];
      }
    }
    if(std2 != 0.0 && pro2 != 0.0) {
      ratstd = std1/std2;
      ratpro = pro1/pro2;
      chi2 += pow(fabs(1.0-ratstd/ratpro),2.0);
    }
  }

  free_vector(x,0,N);
  free_vector(y,0,N);
  return(chi2);
}

float sptHeII(float spcode)
{
  extern float *X,*Y;
  extern int N;
  extern int sf;
  extern float Lumcode;
  extern float whigh;
  extern double shot;
  float *x,*y;
  int i,n;
  float chi;
  float chi2 = 0.0;
  float ssumc1,ssumI,ssumc2,psumc1,psumI,psumc2;
  float sindx,pindx;

  x = vector(0,N);
  y = vector(0,N);

  if(spcode <= shot) return(1.0e+05);

  sp2class(spcode,Lumcode,0.0,x,y,&n);
 
  ssumc1 = ssumI = ssumc2 = psumc1 = psumI = psumc2 = 0.0;

  if(spcode >= 7.5) return(1.0e+30); /* only for early B stars and earlier */
  for(i=0;i<n;i++) {
    if(x[i] >= 4511.0 && x[i] <= 4531.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4531.0 && x[i] <= 4551.0) {
      ssumI += y[i];
      psumI += Y[i];
    }
    if(x[i] >= 4551.0 && x[i] <= 4571.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  sindx = ssumI/(0.5*ssumc1 + 0.5*ssumc2);
  pindx = psumI/(0.5*psumc1 + 0.5*psumc2);
  chi2 += (sindx - pindx)*(sindx-pindx);

  return(chi2);
}


/* This function calculates the spectral type based on He I strengths in
comparison with standards.  To help bridge the He I maximum at B2, it also
considers other temperature sensitive criteria, such as the Si II 4128-30/
He I 4144 ratio, the presence of the feature near 4071A, the strength
of the C II line at 4267, and the strength of the Mg II 4481 line. 
It also includes a contribution to chi2 based on the hydrogen lines */   
float sptHeImet(float spcode)
{
  extern float *X,*Y;
  extern int N;
  extern int sf;
  extern float whigh;
  extern float Lumcode;
  float *x,*y;
  int i,n;
  float chi;
  float chi2 = 0.0;
  float ssumc1,ssumI,ssumc2,psumc1,psumI,psumc2,ratio1,ratio2;
  float sindx,pindx;
  float dif;

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,Lumcode,0.0,x,y,&n);
 
  ssumc1 = ssumI = ssumc2 = psumc1 = psumI = psumc2 = 0.0;

  if(spcode >= 20.0) return(1.0e+10); /* only for B stars and earlier */
  chi2 = hydrogen_profile_hot(spcode);
  /* He I 4026 line */
  for(i=0;i<n;i++) {
    if(x[i] >= 4014.0 && x[i] <= 4019.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4019.0 && x[i] <= 4033.0) {
      ssumI += y[i];
      psumI += Y[i];
    }
    if(x[i] >= 4033.0 && x[i] <= 4038.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  ssumc1 /= 5.0;
  ssumc2 /= 5.0;
  psumc1 /= 5.0;
  psumc2 /= 5.0;
  ssumI /= 14.0;
  psumI /= 14.0;
  sindx = ssumI/(0.5*ssumc1+0.5*ssumc2);
  pindx = psumI/(0.5*psumc1+0.5*psumc2);
  chi2 += (sindx-pindx)*(sindx-pindx);

  /* He I 4387 */
ssumc1 = ssumI = ssumc2 = psumc1 = psumI = psumc2 = 0.0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4365.0 && x[i] <= 4377.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4377.0 && x[i] <= 4397.0) {
      ssumI += y[i];
      psumI += Y[i];
    }
    if(x[i] >= 4397.0 && x[i] <= 4409.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  ssumc1 /= 12.0;
  ssumc2 /= 12.0;
  psumc1 /= 12.0;
  psumc2 /= 12.0;
  ssumI /= 20.0;
  psumI /= 20.0;
  sindx = ssumI/(0.5*ssumc1+0.5*ssumc2);
  pindx = psumI/(0.5*psumc1+0.5*psumc2);
  chi2 += (sindx-pindx)*(sindx-pindx);

  /* He I 4471 */
ssumc1 = ssumI = ssumc2 = psumc1 = psumI = psumc2 = 0.0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4439.0 && x[i] <= 4454.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4465.0 && x[i] <= 4476.0) {
      ssumI += y[i];
      psumI += Y[i];
    }
    if(x[i] >= 4486.0 && x[i] <= 4501.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  ssumc1 /= 15.0;
  ssumc2 /= 15.0;
  psumc1 /= 15.0;
  psumc2 /= 15.0;
  ssumI /= 11.0;
  psumI /= 11.0;
  sindx = ssumI/(0.5*ssumc1+0.5*ssumc2);
  pindx = psumI/(0.5*psumc1+0.5*psumc2);
  chi2 += (sindx-pindx)*(sindx-pindx);

  /* Si II 4128-30 ratioed with He I 4144 */
ssumc1 = ssumI = ssumc2 = psumc1 = psumI = psumc2 = 0.0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4124.0 && x[i] <= 4135.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4136.0 && x[i] <= 4150.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  ratio1 = ssumc1/ssumc2;
  ratio2 = psumc1/psumc2;
  chi2 += 2.0*(ratio1-ratio2)*(ratio1-ratio2);

  /* The feature at 4071 in ratio with nearby continuum */
ssumc1 = ssumI = ssumc2 = psumc1 = psumI = psumc2 = 0.0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4053 && x[i] <= 4089.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4065.0 && x[i] <= 4080.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  ssumc1 /= 36.0;
  psumc1 /= 36.0;
  ssumc2 /= 15.0;
  psumc2 /= 15.0;
  ratio1 = ssumc1/ssumc2;
  ratio2 = psumc1/psumc2;
  chi2 += 2.0*(ratio1-ratio2)*(ratio1-ratio2);

  /* The C II 4267 line in ratio with nearby continua */
ssumc1 = ssumI = ssumc2 = psumc1 = psumI = psumc2 = 0.0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4248 && x[i] <= 4260.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4260.0 && x[i] <= 4272.0) {
      ssumI += y[i];
      psumI += Y[i];
    }
    if(x[i] >= 4272.0 && x[i] <= 4284.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  sindx = ssumI/(0.5*ssumc1+0.5*ssumc2);
  pindx = psumI/(0.5*psumc1+0.5*psumc2);
  chi2 += 2.0*(sindx-pindx)*(sindx-pindx);

  /* The ratio of He I 4471 to Mg II 4481 */
ssumc1 = ssumI = ssumc2 = psumc1 = psumI = psumc2 = 0.0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4463.0 && x[i] <= 4475.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4475.0 && x[i] <= 4486.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  ratio1 = ssumc1/ssumc2;
  ratio2 = psumc1/psumc2;
  chi2 += 5.0*(ratio1-ratio2)*(ratio1-ratio2);

  return(chi2);
}

/* This function judges the spectral type based only on the helium I lines */

float sptHeI(float spcode, float lumcode)
{
  extern float *X,*Y;
  extern int N;
  extern int sf;
  extern float whigh;
  float *x,*y;
  int i,n;
  float chi;
  float chi2 = 0.0;
  float ssumc1,ssumI,ssumc2,psumc1,psumI,psumc2,ratio1,ratio2;
  float sindx,pindx;
  float dif;

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,lumcode,0.0,x,y,&n);
 
  ssumc1 = ssumI = ssumc2 = psumc1 = psumI = psumc2 = 0.0;

  if(spcode >= 19.0) return(1.0e+10); /* only for B stars and earlier */
  /* He I 4026 line */
  for(i=0;i<n;i++) {
    if(x[i] >= 4014.0 && x[i] <= 4019.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4019.0 && x[i] <= 4033.0) {
      ssumI += y[i];
      psumI += Y[i];
    }
    if(x[i] >= 4033.0 && x[i] <= 4038.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  ssumc1 /= 5.0;
  ssumc2 /= 5.0;
  psumc1 /= 5.0;
  psumc2 /= 5.0;
  ssumI /= 14.0;
  psumI /= 14.0;
  sindx = ssumI/(0.5*ssumc1+0.5*ssumc2);
  pindx = psumI/(0.5*psumc1+0.5*psumc2);
  chi2 += (sindx-pindx)*(sindx-pindx);

  /* He I 4387 */
ssumc1 = ssumI = ssumc2 = psumc1 = psumI = psumc2 = 0.0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4365.0 && x[i] <= 4377.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4377.0 && x[i] <= 4397.0) {
      ssumI += y[i];
      psumI += Y[i];
    }
    if(x[i] >= 4397.0 && x[i] <= 4409.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  ssumc1 /= 12.0;
  ssumc2 /= 12.0;
  psumc1 /= 12.0;
  psumc2 /= 12.0;
  ssumI /= 20.0;
  psumI /= 20.0;
  sindx = ssumI/(0.5*ssumc1+0.5*ssumc2);
  pindx = psumI/(0.5*psumc1+0.5*psumc2);
  chi2 += (sindx-pindx)*(sindx-pindx);

  /* He I 4471 */
ssumc1 = ssumI = ssumc2 = psumc1 = psumI = psumc2 = 0.0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4439.0 && x[i] <= 4454.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4465.0 && x[i] <= 4476.0) {
      ssumI += y[i];
      psumI += Y[i];
    }
    if(x[i] >= 4486.0 && x[i] <= 4501.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  ssumc1 /= 15.0;
  ssumc2 /= 15.0;
  psumc1 /= 15.0;
  psumc2 /= 15.0;
  ssumI /= 11.0;
  psumI /= 11.0;
  sindx = ssumI/(0.5*ssumc1+0.5*ssumc2);
  pindx = psumI/(0.5*psumc1+0.5*psumc2);
  chi2 += (sindx-pindx)*(sindx-pindx);

  return(chi2);
}
