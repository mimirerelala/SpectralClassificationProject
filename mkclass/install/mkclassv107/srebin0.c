#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "util.h"
float *x2,*y2;
void getspectrum();
void rebin();

main(int argc, char *argv[])
{
  float *xo,*yo;
  int l1,l2;
  int i,k,N,ni;
  double start,end,dw;
  char infile[80],outfile[80];
  FILE *in,*out;

  x2 = vector(0,500000);
  y2 = vector(0,500000);
  xo = vector(0,500000);
  yo = vector(0,500000);

  if(argc != 6) {
    printf("\nEnter name of input spectrum > ");
    ggets(infile);
    if((in = fopen(infile,"r")) == NULL) {
      printf("\nCannot find input file\n");
      exit(1);
    }
    fclose(in);
    printf("\nEnter name of output rebinned spectrum > ");
    ggets(outfile);
    out = fopen(outfile,"w");

    printf("\nEnter lambda(start),lambda(end),delta_lambda > ");
    ni = scanf("%lf,%lf,%lf",&start,&end,&dw);
  } else {
    strcpy(infile,argv[1]);
    strcpy(outfile,argv[2]);
    out = fopen(outfile,"w");
    start = atof(argv[3]);
    end = atof(argv[4]);
    dw = atof(argv[5]);
  }

  k = 0;
  getspectrum(infile,xo,yo,&k);
  rebin(xo,yo,k,x2,y2,&l1,&l2,start,end,dw);
  N = (end-start)/dw;  

  for(i=0;i<N;i++) {
    /*    if(y2[i] == 0.0) continue; */
    fprintf(out,"%lf %lg\n",x2[i],y2[i]);
  }

  free_vector(x2,0,500000);
  free_vector(y2,0,500000);
  free_vector(xo,0,500000);
  free_vector(yo,0,500000);
  fclose(out);
}

void rebin(x,y,k,x2,y2,l1,l2,start,end,dw)
float *x,*y;
float *x2,*y2;
double start,end,dw;
int k,*l1,*l2;
{
    double wave;
    int l,j,i;
    char tmp[10];

    wave = start;
    *l2 = 0;
    while(wave <= end) {
      x2[*l2] = wave;
      wave += dw;
      (*l2)++;
    }
    wave = start;

    *l1 = 0;
    if(start > x[0]) *l1 = 0;
    else {
      while(wave <= x[0]) {
	wave += dw;
	(*l1)++;
      }
    }

    while(end > x[k-1]) {
      end -= dw;
      (*l2)--;
    }

    l = 1;
    wave = start + (*l1)*dw;
    for(i=*l1;i<=*l2;i++) {
      for(j=l;j<k;j++) {
	if(start > x[j-1] && start >= x[j]) continue;
	if(wave <= x[j] && wave > x[j-1]) {
	   l = j;
	   break;
	}
      }
      y2[i] = y[l-1] + (y[l]-y[l-1])*(wave-x[l-1])/(x[l]-x[l-1]);
      wave += dw;
    }
}

