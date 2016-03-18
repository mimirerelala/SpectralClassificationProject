// Applies various criteria for the classification of K and M-type stars
#include <stdio.h>
#include <math.h>
#include "util.h"
void sp2class(float spcode, float lumcode, float metcode, 
	      float *x, float *y, int *n);

float sptKM(float spcode)
{
  extern float *X,*Y;
  extern float Lumcode;
  extern int N;
  float *x,*y;
  float chi,sumy,sumY,avey,aveY,sigy,sigY;
  float chi2 = 0.0;
  float ch2 = 0.0;
  float pro1,pro2,std1,std2,ratpro,ratstd;
  float proc1,proc2,proc3,stdc1,stdc2,stdc3,proCaI,proG,stdCaI,stdG;
  float stdMgI,proMgI,istdMgI,iproMgI;
  float stdMgH,proMgH,istdMgH,iproMgH,stdCaOH,proCaOH,istdCaOH,iproCaOH;
  float iproG,iproCaI,istdG,istdCaI;
  int i,nq,n;

  x = vector(0,N);
  y = vector(0,N);

  sp2class(spcode,Lumcode,0.0,x,y,&n);

  pro1=pro2=std1=std2=ratpro=ratstd=0.0;
  proc1=proc2=proc3=stdc1=stdc2=stdc3=proCaI=proG=stdCaI=stdG=0.0;
  chi2 = 0.0;
  for(i=0;i<n;i++) {
    if(x[i] >= 3900.0 && x[i] <= 4250.0) ch2 += (Y[i]-y[i])*(Y[i]-y[i]);
  }
  ch2 /= 350.0;

  /* Early K-type star temperature classification */
  if(spcode <= 39.0) {
  /* Ca I, G-band indices */
    for(i=0;i<n;i++) {
      if(x[i] >= 4202.5 && x[i] <= 4220.5) {
	stdc1 += y[i];
	proc1 += Y[i];
      }
      if(x[i] >= 4225.5 && x[i] <= 4227.5) {
	stdCaI += y[i];
	proCaI += Y[i];
      }
      if(x[i] >= 4233.0 && x[i] <= 4248.0) {
	stdc2 += y[i];
	proc2 += Y[i];
      }
      if(x[i] >= 4297.0 && x[i] <= 4314.0) {
	stdG += y[i];
	proG += Y[i];
      }
      if(x[i] >= 4355.0 && x[i] <= 4378.0) {
	stdc3 += y[i];
	proc3 += Y[i];
      }
    }
    stdc1 /= 18.0;
    proc1 /= 18.0;
    stdCaI /= 2.0;
    proCaI /= 2.0;
    stdc2 /= 15.0;
    proc2 /= 15.0;
    stdG /= 17.0;
    proG /= 17.0;
    stdc3 /= 15.0;
    proc3 /= 15.0;
    istdG = stdG/(0.484*stdc2 + 0.516*stdc3);
    iproG = proG/(0.484*proc2 + 0.516*proc3);
    istdCaI = stdCaI/(stdc1+stdc2);
    iproCaI = proCaI/(proc1+proc2);
    chi2 += pow(fabs(1.0-istdG/iproG),2.0);
    chi2 += pow(fabs(1.0-istdCaI/iproCaI),2.0);

    /* Mg I index */


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
    if(iproMgI != 0.0) chi2 += 2.0*pow(fabs(1.0-istdMgI/iproMgI),2.0);

    /* MgH band */

    stdc1 = stdc2 = proc1 = proc2 = stdMgH = proMgH = 0.0;

    for(i=0;i<n;i++) {
      if(x[i] >= 4722.0 && x[i] <= 4750.0) {
	stdc1 += y[i];
	proc1 += Y[i];
      }
      if(x[i] >= 4750.0 && x[i] <= 4791.0) {
	stdMgH += y[i];
	proMgH += Y[i];
      }
      if(x[i] >= 4791.0 && x[i] <= 4816.0) {
	stdc2 += y[i];
	proc2 += Y[i];
      }
    }
    istdMgH = stdMgH/(stdc1+stdc2);
    iproMgH = proMgH/(proc1+proc2);
    if(iproMgH != 0.0) chi2 += 2.0*pow(fabs(1.0-istdMgH/iproMgH),2.0);

    chi2 /= 6.0;

  } else { 

    /* Mg I index */

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
    if(iproMgI != 0.0) chi2 += 2.0*pow(fabs(1.0-istdMgI/iproMgI),2.0);

    /* MgH band */

    stdc1 = stdc2 = proc1 = proc2 = stdMgH = proMgH = 0.0;

    for(i=0;i<n;i++) {
      if(x[i] >= 4722.0 && x[i] <= 4750.0) {
	stdc1 += y[i];
	proc1 += Y[i];
      }
      if(x[i] >= 4750.0 && x[i] <= 4791.0) {
	stdMgH += y[i];
	proMgH += Y[i];
      }
      if(x[i] >= 4791.0 && x[i] <= 4816.0) {
	stdc2 += y[i];
	proc2 += Y[i];
      }
    }
    istdMgH = stdMgH/(stdc1+stdc2);
    iproMgH = proMgH/(proc1+proc2);
    if(iproMgH != 0.0) chi2 += 2.0*pow(fabs(1.0-istdMgH/iproMgH),2.0);

    /* CaOH band */

    stdc1 = proc1 = stdCaOH = proCaOH = 0.0;

    for(i=0;i<n;i++) {
      if(x[i] >= 5499.0 && x[i] <= 5524.0) {
	stdc1 += y[i];
	proc1 += Y[i];
      }
      if(x[i] >= 5532.5 && x[i] <= 5548.0) {
	stdCaOH += y[i];
	proCaOH += Y[i];
      }
    }
    istdCaOH = stdCaOH/stdc1;
    iproCaOH = proCaOH/proc1;
    if(iproCaOH != 0.0) chi2 += 2.0*pow(fabs(1.0-istdCaOH/iproCaOH),2.0); 

    /* TiO band head measures */

  pro1 = pro2 = std1 = std2 = 0.0;

  for(i=0;i<n;i++) {
    if(x[i] >= 4731.0 && x[i] <= 4751.0) {
      std1 += y[i];
      pro1 += Y[i];
    }
    if(x[i] >= 4766.0 && x[i] <= 4786.0) {
      std2 += y[i];
      pro2 += Y[i];
    }
  }
  if(std2 != 0.0 && pro2 != 0.0) {
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    chi2 += pow(fabs(1.0-ratstd/ratpro),2.0);
  }

  for(i=0;i<n;i++) {
    if(x[i] >= 4929.0 && x[i] <= 4949.0) {
      std1 += y[i];
      pro1 += Y[i];
    }
    if(x[i] >= 4959.0 && x[i] <= 4979.0) {
      std2 += y[i];
      pro2 += Y[i];
    }
  }
  if(std2 != 0.0 && pro2 != 0.0) {
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    chi2 += pow(fabs(1.0-ratstd/ratpro),2.0);
  }

  for(i=0;i<n;i++) {
    if(x[i] >= 5140.0 && x[i] <= 5160.0) {
      std1 += y[i];
      pro1 += Y[i];
    }
    if(x[i] >= 5172.0 && x[i] <= 5192.0) {
      std2 += y[i];
      pro2 += Y[i];
    }
  }
  if(std2 != 0.0 && pro2 != 0.0) {
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    chi2 += pow(fabs(1.0-ratstd/ratpro),2.0);
  }

  for(i=0;i<n;i++) {
    if(x[i] >= 5424.0 && x[i] <= 5444.0) {
      std1 += y[i];
      pro1 += Y[i];
    }
    if(x[i] >= 5456.0 && x[i] <= 5476.0) {
      std2 += y[i];
      pro2 += Y[i];
    }
  }
  if(std2 != 0.0 && pro2 != 0.0) {
    ratstd = std1/std2;
    ratpro = pro1/pro2;
    chi2 += pow(fabs(1.0-ratstd/ratpro),2.0);
  }
  chi2 /= 8.0;
  }
  free_vector(x,0,N); 
  free_vector(y,0,N);

  return(chi2+ch2);
}
