// This file contains functions for applying luminosity criteria
#include <stdio.h>
#include <math.h>
#include "util.h"
void sp2class(float spcode, float lumcode, float metcode, 
	      float *x, float *y, int *n);
float hydrogen_luminosity(float *x, float *y, int n);

float lummin(float lumcode)
{
  extern float *X,*Y;
  extern float Spcode;
  extern int N;
  float *x,*y;
  float chi;
  float chi2 = 0.0;
  int i,n,nq;

  x = vector(0,N);
  y = vector(0,N);

  sp2class(Spcode,lumcode,0.0,x,y,&n);

  if(Spcode < 21.0) {
    /* A-type star luminosity criteria -- Hydrogen lines */
    for(i=0;i<n;i++) {
      if((x[i] > 4091.0 && x[i] <= 4110.0) || (x[i] >= 4330.0 && x[i] <= 4350.0) || (x[i] >= 4850.0 && x[i] <= 4870.0)) {
	chi = y[i] - Y[i];
	chi2 += chi*chi;
      }
    }
    chi2 /= 60.0;
  } else if(Spcode < 27.5) {
    /* Early F-star luminosity criteria -- FeII TiII lines */
    for(i=0;i<n;i++) {
      if((x[i] > 4170.0 && x[i] < 4180.0) || (x[i] > 4390 && x[i] < 4600)) {
	chi = y[i] - Y[i];
	chi2 += chi*chi;
      }
    } 
    chi2 /= 220.0;
  } else if(Spcode < 31.5) {
  /* Late F and early G luminosity criteria -- FeII, TiII, SrII lines */
    for(i=0;i<n;i++) {
      if((x[i] > 4075 && x[i] < 4080) || (x[i] > 4213 && x[i] < 4219) || (x[i] > 4170.0 && x[i] < 4180.0) || (x[i] > 4390 && x[i] < 4600)) {
	chi = y[i] - Y[i];
	chi2 += chi*chi;
      }
    }
    chi2 /= 231.0;
  } else if(Spcode < 36.0) {
      /* Late G and early K-type luminosity criteria -- CN band & Sr II */
    for(i=0;i<n;i++) {
      if((x[i] > 4075 && x[i] < 4080) || (x[i] > 4110.0 && x[i] < 4216.0)) {
	chi = y[i] - Y[i];
	chi2 += chi*chi;
      }
    }
    chi2 /= 111.0;
  } else if(Spcode <= 40.5) {
  /* Late K-type luminosity criteria -- CN band, Ca I 4226 line, CaH
     band, and 5250/5296 */
  for(i=0;i<n;i++) {
    nq = 0;
    if(x[i] > 4000 && x[i] < 5600) {
      chi = y[i]-Y[i];
      chi2 += chi*chi;
      nq++;
    }
    if(x[i] > 4086 && x[i] < 4243) {
      chi = y[i] - Y[i];
      chi2 += 3.0*chi*chi;
      nq++;
    }
    if(x[i] >= 4700.0 && x[i] <= 4950.0) {
	chi = y[i] - Y[i];
	chi2 += 3.0*chi*chi;
	nq++;
    }
    if(x[i] >= 5100.0 && x[i] <= 5330.0) {
	chi = y[i] - Y[i];
	chi2 += 3.0*chi*chi;
	nq++;
    }
  }
  } else if(Spcode <= 45.0) {
    /* M-type luminosity criteria -- Ca I 4226 line, MgH band, 
       CaOH band and 5250/5296 */
    if(x[i] > 4086 && x[i] < 4243) {
      chi = y[i] - Y[i];
      chi2 += 3.0*chi*chi;
      nq++;
    }
    if(x[i] >= 4700.0 && x[i] <= 4950.0) {
	chi = y[i] - Y[i];
	chi2 += 3.0*chi*chi;
	nq++;
    }
    if(x[i] >= 5100.0 && x[i] <= 5330.0) {
	chi = y[i] - Y[i];
	chi2 += 3.0*chi*chi;
	nq++;
    }
    if(x[i] >= 5532.0 && x[i] <= 5548.0) {
      chi = y[i] - Y[i];
      chi2 += 3.0*chi*chi;
      nq++;
    }
  }

  chi2 /= (float)nq;


  free_vector(x,0,N);
  free_vector(y,0,N);
  return(chi2);
} 

float lumratiomin(float lumcode)
{
  extern float *X,*Y;
  extern float Spcode;
  extern int Ba;
  extern int N;
  extern FILE *Log;
  extern float whigh;
  float *x,*y;
  float pro1,pro2,std1,std2,ratpro,ratstd;
  float stdc1,stdc2,proc1,proc2,stdMgI,proMgI;
  float istdMgI,iproMgI;
  float chi;
  float chi2 = 0.0;
  float cont1p,cont1s,cont2p,cont2s,contp,conts;
  float sum1,sum2;
  int i,n,k,k1,k2;

  x = vector(0,N);
  y = vector(0,N);

  sp2class(Spcode,lumcode,0.0,x,y,&n);

  chi2 = 0.0;
  if(Spcode  >= 7.0 && Spcode < 10.0) {
    chi2 = hydrogen_luminosity(x,y,n);
    /* We also use with heavier weight luminosity criteria based on
       O II lines and Si III lines ratioed to He I lines */
    pro1 = pro2 = std1 = std2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4379.0 && x[i] <= 4395.0) {
	std1 += y[i];
	pro1 += Y[i];
      }
      if(x[i] >= 4407.0 && x[i] <= 4422.0) {
	std2 += y[i];
	pro2 += y[i];
      }
    }
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    chi2 += 4.0*pow(fabs(1.0 - ratstd/ratpro),2.0);

    pro1 = pro2 = std1 = std2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4463.0 && x[i] <= 4476.0) {
	std1 += y[i];
	pro1 += Y[i];
      }
      if(x[i] >= 4545.0 && x[i] <= 4581.0) {
	std2 += y[i];
	pro2 += y[i];
      }
    }
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    chi2 += 4.0*pow(fabs(1.0 - ratstd/ratpro),2.0);
  }

  if(Spcode >= 10.0 && Spcode < 20.5) {
    /* A-type star luminosity criteria -- Hydrogen lines */
    chi2 = hydrogen_luminosity(x,y,n);
    /*  } else if(Spcode <= 20.5) {
    chi2 = 0.0;
    chi2 = 0.005*hydrogen_luminosity(x,y,n); ???
    pro1 = pro2 = std1 = std2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4166.0 && x[i] <= 4180.0) {
	std1 += y[i];
	pro1 += Y[i];
      }
      if(x[i] >= 4146.0 && x[i] <= 4162.0) {
	std2 += y[i];
	pro2 += Y[i];
      }
    }
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    chi2 += pow(fabs(1.0 - ratstd/ratpro),2.0);

    pro1 = pro2 = std1 = std2 = 0.0;
    for(i=0;i<n;i++) {
      if((x[i] >= 4410.0 && x[i] <= 4419.0) || 
	 (x[i] >= 4438.0 && x[i] <= 4448.0)) {
	std1 += y[i];
	pro1 += Y[i];
      }
      if(x[i] >= 4419.0 && x[i] <= 4438.0) {
	std2 += y[i];
	pro2 += Y[i];
      }
    }
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    chi2 += pow(fabs(1.0 - ratstd/ratpro),2.0); */
  } else if(Spcode < 27.5) {
    /* Early F-star luminosity criteria -- FeII TiII lines */
    pro1 = pro2 = std1 = std2 = 0.0;
    chi2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4166.0 && x[i] <= 4180.0) {
	std1 += y[i];
	pro1 += Y[i];
      }
      if(x[i] >= 4146.0 && x[i] <= 4162.0) {
	std2 += y[i];
	pro2 += Y[i];
      }
    }
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    chi2 += pow(fabs(1.0 - ratstd/ratpro),2.0);

    pro1 = pro2 = std1 = std2 = 0.0;
    for(i=0;i<n;i++) {
      if((x[i] >= 4410.0 && x[i] <= 4419.0) || 
	 (x[i] >= 4438.0 && x[i] <= 4448.0)) {
	std1 += y[i];
	pro1 += Y[i];
      }
      if(x[i] >= 4419.0 && x[i] <= 4438.0) {
	std2 += y[i];
	pro2 += Y[i];
      }
    }
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    chi2 += pow(fabs(1.0 - ratstd/ratpro),2.0);
  } else if(Spcode < 31.5) {
  /* Late F and early G luminosity criteria -- FeII, TiII, SrII lines */
    pro1 = pro2 = std1 = std2 = 0.0;
    chi2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4166.0 && x[i] <= 4180.0) {
	std1 += y[i];
	pro1 += Y[i];
      }
      if(x[i] >= 4146.0 && x[i] <= 4162.0) {
	std2 += y[i];
	pro2 += Y[i];
      }
    }
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    chi2 += pow(fabs(1.0 - ratstd/ratpro),2.0);

    pro1 = pro2 = std1 = std2 = 0.0;
    for(i=0;i<n;i++) {
      if((x[i] >= 4410.0 && x[i] <= 4419.0) || 
	 (x[i] >= 4438.0 && x[i] <= 4448.0)) {
	std1 += y[i];
	pro1 += Y[i];
      }
      if(x[i] >= 4419.0 && x[i] <= 4438.0) {
	std2 += y[i];
	pro2 += Y[i];
      }
    }
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    chi2 += pow(fabs(1.0 - ratstd/ratpro),2.0);
  
    /* Don't include the Sr II lines in the luminosity determination if
       the star is a Barium star */

   if(Ba == 0) {
     pro1 = pro2 = std1 = std2 = 0.0;
     for(i=0;i<n;i++) {
        if(x[i] >= 4041.0 && x[i] <= 4058.0) {
  	  std2 += y[i];
	  pro2 += Y[i];
        }
        if(x[i] >= 4073.0 && x[i] <= 4082.0) {
	  std1 += y[i];
	  pro1 += Y[i];
        }
      }
      ratstd = std1/std2;
      ratpro = pro1/pro2;
      chi2 += pow(fabs(1.0 - ratstd/ratpro),2.0);

      pro1 = pro2 = std1 = std2 = 0.0;
      for(i=0;i<n;i++) {
        if(x[i] >= 4212.0 && x[i] <= 4218.0) {
	  std1 += y[i];
	  pro1 += Y[i];
        }
        if(x[i] >= 4220.0 && x[i] <= 4230.0) {
	  std2 += y[i];
	  pro2 += Y[i];
        }
      }
      ratstd = std1/std2;
      ratpro = pro1/pro2;
      chi2 += pow(fabs(1.0 - ratstd/ratpro),2.0);
   }    
 
  } else if(Spcode < 39.0) {
      /* Late G and early K-type luminosity criteria -- Sr II */
   pro1 = pro2 = std1 = std2 = 0.0;
   chi2 = 0.0;
   if(Ba == 0) {
   for(i=0;i<n;i++) {
      if(x[i] >= 4041.0 && x[i] <= 4058.0) {
	std2 += y[i];
	pro2 += Y[i];
      }
      if(x[i] >= 4073.0 && x[i] <= 4082.0) {
	std1 += y[i];
	pro1 += Y[i];
      }
    }
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    chi2 += pow(fabs(1.0 - ratstd/ratpro),2.0);

    pro1 = pro2 = std1 = std2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4212.0 && x[i] <= 4218.0) {
	std1 += y[i];
	pro1 += Y[i];
      }
      if(x[i] >= 4220.0 && x[i] <= 4230.0) {
	std2 += y[i];
	pro2 += Y[i];
      }
    }
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    chi2 += pow(fabs(1.0 - ratstd/ratpro),2.0);
   }


    /* 4215 CN band */
    pro1 = pro2 = std1 = std2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4140.0 && x[i] <= 4210.0) {
	std1 += y[i];
	pro1 += Y[i];
      }
      if((x[i] >= 4219.0 && x[i] <= 4264.0) || (x[i] >= 4043 && x[i] <= 4088.0)) {
	std2 += y[i];
	pro2 += Y[i];
      }
    }
    ratstd = 9*std1/(7*std2);
    ratpro = 9*pro1/(7*pro2);
    chi2 += 3.0*pow(fabs(1.0 - ratstd/ratpro),2.0);

    if(Spcode >= 32.0 && whigh > 5278.0) {
    /* 5250/5269 */
      pro1 = pro2 = std1 = std2 = 0.0;
      for(i=0;i<n;i++) {
        if(x[i] >= 5239.0 && x[i] <= 5256.0) {
	  std1 += y[i];
	  pro1 += Y[i];
        }
        if(x[i] >= 5256.0 && x[i] <= 5278.0) {
	  std2 += y[i];
	  pro2 += Y[i];
        }
      }
      ratstd = std1/std2;
      ratpro = pro1/pro2;
      if(ratpro != 0.0) chi2 += 5.0*(1.0 - ratstd/ratpro)*(1.0 - ratstd/ratpro);
    }
  }
  else if(Spcode >= 39.0) {
    chi2 = 0.0;
    /* The 5050 region */
   pro1 = pro2 = std1 = std2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4980.0 && x[i] <= 5030.0) {
	std1 += y[i];
	pro1 += Y[i];
      }
      if(x[i] >= 5100.0 && x[i] <= 5150.0) {
	std2 += y[i];
	pro2 += Y[i];
      }
    }
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    if(ratpro != 0.0) chi2 += 0.0*(1.0 - ratstd/ratpro)*(1.0 - ratstd/ratpro);

    /* The MgH band */
   pro1 = pro2 = std1 = std2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 4720.0 && x[i] <= 4745.0) {
	std1 += y[i];
	pro1 += Y[i];
      }
      if(x[i] >= 4760.0 && x[i] <= 4785.0) {
	std2 += y[i];
	pro2 += Y[i];
      }
    }
    ratstd = std2/std1;
    ratpro = pro2/pro1;
    if(ratpro != 0.0) chi2 += (1.0 - ratstd/ratpro)*(1.0 - ratstd/ratpro);
    /* 5250/5269 */
   pro1 = pro2 = std1 = std2 = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 5239.0 && x[i] <= 5256.0) {
	std1 += y[i];
	pro1 += Y[i];
      }
      if(x[i] >= 5256.0 && x[i] <= 5278.0) {
	std2 += y[i];
	pro2 += Y[i];
      }
    }
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    if(ratpro != 0.0) chi2 += 5.0*(1.0 - ratstd/ratpro)*(1.0 - ratstd/ratpro);

    /* MgI band */ 
    stdc1 = stdc2 = proc1 = proc2 = stdMgI = proMgI = 0.0;
    for(i=0;i<n;i++) {
      if(x[i] >= 5004.0 && x[i] <= 5058.0) {
	stdc1 += y[i];
	proc1 += Y[i];
      }
      if(x[i] >= 5158.0 && x[i] <= 5212.0) {
	stdMgI += y[i];
	proMgI += Y[i];
      }
      if(x[i] >= 5312.0 && x[i] <= 5366.0) {
	stdc2 += y[i];
	proc2 += Y[i];
      }
    }
    istdMgI = stdMgI/(stdc1+stdc2);
    iproMgI = proMgI/(proc1+proc2);
    if(iproMgI != 0.0) chi2 += 2.0*(1.0-istdMgI/iproMgI)*(1.0-istdMgI/iproMgI);
  }

  free_vector(x,0,N);
  free_vector(y,0,N);

  return(chi2);
}

float hydrogen_luminosity(float *x, float *y, int n)
{
  extern float *X,*Y;
  extern int N;
  extern int sf;
  float yr,Yr,scont,pcont;
  int i;
  float chi;
  float chi2 = 0.0;
  float ssumc1,ssumc2,psumc1,psumc2;
  float sindx,pindx;

  /* H-delta line */

  chi2 = 0.0;
  ssumc1 = ssumc2 = psumc1 = psumc2 = 0.0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4012.0 && x[i] <= 4052.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4152.0 && x[i] <= 4192.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  ssumc1 /= 40.0;
  psumc1 /= 40.0;
  ssumc2 /= 40.0;
  psumc2 /= 40.0;
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
  for(i=0;i<n;i++) {
    if(x[i] >= 4250.0 && x[i] <= 4290.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4390.0 && x[i] <= 4430.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  ssumc1 /= 40.0;
  psumc1 /= 40.0;
  ssumc2 /= 40.0;
  psumc2 /= 40.0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4300.0 && x[i] <= 4380.0) {
      scont = ssumc1 + (ssumc2-ssumc1)*(x[i] - 4270.0)/(4410.0-4270.0);
      pcont = psumc1 + (psumc2-psumc1)*(x[i] - 4270.0)/(4410.0-4270.0);
      yr = y[i]/scont;
      Yr = Y[i]/pcont;
      chi2 += (yr-Yr)*(yr-Yr);
    }
  }

  /* H Beta */
  ssumc1 = ssumc2 = psumc1 = psumc2 = 0.0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4791.0 && x[i] <= 4831.0) {
      ssumc1 += y[i];
      psumc1 += Y[i];
    }
    if(x[i] >= 4911.0 && x[i] <= 4951.0) {
      ssumc2 += y[i];
      psumc2 += Y[i];
    }
  }
  ssumc1 /= 40.0;
  psumc1 /= 40.0;
  ssumc2 /= 40.0;
  psumc2 /= 40.0;
  for(i=0;i<n;i++) {
    if(x[i] >= 4841.0 && x[i] <= 4901.0) {
      scont = ssumc1 + (ssumc2-ssumc1)*(x[i] - 4811.0)/(4931.0-4811.0);
      pcont = psumc1 + (psumc2-psumc1)*(x[i] - 4811.0)/(4931.0-4811.0);
      yr = y[i]/scont;
      Yr = Y[i]/pcont;
      chi2 += (yr-Yr)*(yr-Yr);
    }
  }

  return(chi2);
}

float ratioCaIFeII()
{
  extern float *X,*Y;
  extern int N;
  float psumc1,psumc2,cont,ratio;
  float count = 0;
  float count1 = 0;
  float count2 = 0;
  int i;

  /* Ca I 4226 ratioed with Fe II 4233 */
  psumc1 = psumc2 = cont = 0.0;
  for(i=0;i<N;i++) {
    if((X[i] >= 4219.0 && X[i] <= 4222.0) || (X[i] > 4242.0 && X[i] <= 4245.0)) {
      cont += Y[i];
      count += 1.0;
    }
  }
  cont /= count;
  for(i=0;i<N;i++) {
    if(X[i] >= 4223.0 && X[i] <= 4230.0) {
      psumc1 += cont - Y[i];
      count1 += 1.0;
    }
    if(X[i] >= 4231.0 && X[i] <= 4235.0) {
      psumc2 += cont - Y[i];
      count2 += 1.0;
    }
  }
  psumc1 /= count1;
  psumc2 /= count2;
  ratio = psumc1/psumc2;

  return(ratio);
}
