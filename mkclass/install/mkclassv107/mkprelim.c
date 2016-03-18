#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "util.h"
long size = 2*600000;
double velshift(char *infile, char *template, double start, double end);
void rebin();
void getspectrum();
void gaussj(float **a, int n, float **b, int m);
int approx();

int main(int argc, char *argv[])
{
  char infile[80],template[80],outfile[80],keywords[80];
  double start,end,dw,vel;
  double *wavein,*waveout,*flux,*fluxout;
  double flxmax;
  int i,k,l1=0,l2=0,n;
  int flag_shift = 0;
  int flag_norm = 0;
  double c = 2.99792458e+05;
  FILE *out,*nor;

  if(argc != 8) {
    printf("\nUsage: mkprelim infile template outfile start end dwave keywords\n");
    exit(1);
  }

  wavein = (double*)calloc(size,sizeof(double));
  waveout = (double*)calloc(size,sizeof(double));
  flux = (double*)calloc(size,sizeof(double));
  fluxout = (double*)calloc(size,sizeof(double));

  strcpy(infile,argv[1]);
  strcpy(template,argv[2]);
  strcpy(outfile,argv[3]);
  start = atof(argv[4]);
  end = atof(argv[5]);
  dw = atof(argv[6]);
  strcpy(keywords,argv[7]);

  if(strstr(keywords,"shift") != NULL) flag_shift=1;
  if(strstr(keywords,"norm") != NULL) flag_norm=1;

  n = floor((end-start)/dw);

  out = fopen(outfile,"w");
  if(out == NULL) {
    printf("\nCannot open output file %s\n",outfile);
    exit(1);
  }

  if(flag_norm == 1) {
    getspectrum(infile,wavein,flux,&k);
    flxmax = 0;
    for(i=0;i<k;i++) {
      if(wavein[i] > 4490.0 && wavein[i] < 4520.0) {
	if(flux[i] > flxmax) flxmax = flux[i];
      }
    }
    for(i=0;i<k;i++) flux[i] /= flxmax;
    nor = fopen("normal.dat","w");
    if(nor == NULL) {
       printf("Cannot open normal.dat in mkprelim\n");
       exit(1);
    }
    for(i=0;i<k;i++) fprintf(nor,"%f %f\n",wavein[i],flux[i]);
    fclose(nor);
    strcpy(infile,"normal.dat");
  }

  vel = 0.0;
  if(flag_shift == 1) vel = velshift(infile,template,start+100,end-100);

  k = 0;
  getspectrum(infile,wavein,flux,&k);

  for(i=0;i<k;i++) wavein[i] = wavein[i]/(1.0 + vel/c);

  rebin(wavein,flux,k,waveout,fluxout,&l1,&l2,start,end,dw);

  for(i=0;i<n;i++) fprintf(out,"%7.2f %f\n",waveout[i],fluxout[i]);
  fclose(out); 

  return(0);
}

double velshift(char *infile, char *template, double start, double finish)
{
  double *wtmplt,*ftmplt;
  double *wobs,*fobs;
  double *wtrbn,*ftrbn;
  double *worbn,*forbn;
  double *wtrbn2,*ftrbn2;
  double *wtrbn3,*ftrbn3;
  double chimin = 1.0e+32;
  double chi2,chim1,chi0,chip1,ch,v;
  double vm1,v0,vp1,q;
  double vmin,dwav=0.0;;
  double begin,end,vel;
  double vstep = 2.0,V,Vend,wcenter;
  double wrest;
  FILE *tmplt,*obs,*out,*tst,*shft;
  int i,j,k,l,m,n,l1=0,l2=0,M,N,ni;
  double c = 2.99792458e+05;
  double wavcen,radvel;
  int flag = 1;
  float **A;
  float **a;
  float **y;

  wtmplt = (double*)calloc(size,sizeof(double));
  ftmplt = (double*)calloc(size,sizeof(double));
  wobs = (double*)calloc(size,sizeof(double));
  fobs = (double*)calloc(size,sizeof(double));
  wtrbn = (double*)calloc(size,sizeof(double));
  ftrbn = (double*)calloc(size,sizeof(double));
  worbn = (double*)calloc(size,sizeof(double));
  forbn = (double*)calloc(size,sizeof(double));
  wtrbn2 = (double*)calloc(size,sizeof(double));
  ftrbn2 = (double*)calloc(size,sizeof(double));
  wtrbn3 = (double*)calloc(size,sizeof(double));
  ftrbn3 = (double*)calloc(size,sizeof(double));
  A = matrix(1,3,1,3);
  a = matrix (1,3,1,1);
  y = matrix(1,3,1,1);

  out = fopen("out.out","w");
  tst = fopen("test.out","w");

  k = 0;
  getspectrum(template,wtmplt,ftmplt,&k);
  l = 0;
  getspectrum(infile,wobs,fobs,&l);

  rebin(wobs,fobs,l,worbn,forbn,&l1,&l2,start,finish,0.01);
  rebin(wtmplt,ftmplt,k,wtrbn,ftrbn,&l1,&l2,start,finish,0.01);
  wcenter = (finish+start)/2.0;
  N = (finish-start)/0.01;

  m = 0;
  n = N-1;
  for(i=0;i<N;i++) {
    if(approx(wtrbn[i],start+7.0,0.005) == 1) m = i;
    if(approx(wtrbn[i],finish-7.0,0.005) == 1) n = i;
  }

  V = -12.0*c/wcenter;
  Vend = 12.0*c/wcenter;
  while(V <= Vend) {
    for(i=0;i<N;i++) {
      dwav = wtrbn[i]*V/c;
      wtrbn2[i] = wtrbn[i] + dwav;
      ftrbn2[i] = ftrbn[i];
    }
    for(i=0;i<N;i++) wtrbn3[i] = ftrbn3[i] = 0.0;
    rebin(wtrbn2,ftrbn2,N,wtrbn3,ftrbn3,&l1,&l2,start,finish,0.01);
    /*    if(flag == 1) {
      for(i=0;i<N;i++) fprintf(tst,"%f %f\n",wtrbn3[i],ftrbn3[i]);
      flag = 0;
      } */
    chi2 = 0.0;
    M = 0;
    for(i=m;i<=n;i++) {
      if(ftrbn3[i] == 0.0) continue;
      chi2 += (ftrbn3[i] - forbn[i])*(ftrbn3[i] - forbn[i]);
      M += 1;
    }
    chi2 /= (double)M;
    fprintf(out,"%lf %le\n",V,chi2);
    if(chimin > chi2) {
      chimin = chi2;
      vmin = V;
    }
    /* printf("V = %f  chi2 = %f\n",V,chi2); */
    V += vstep;
  }
  fclose(out);
  out = fopen("out.out","r");
  chim1 = chi0 = chip1 = vm1 = v0 = vp1 = 0.0;
  while(fscanf(out,"%lf %le",&v,&ch) != EOF) {
    chim1 = chi0;
    vm1 = v0;
    chi0 = ch;
    v0 = v;
    if(fabs(v0 - vmin) < vstep/4.0) {
      ni = fscanf(out,"%lf %lf",&vp1,&chip1);
      break;
    }
  }
  A[1][1] = vm1*vm1;
  A[1][2] = vm1;
  A[1][3] = A[2][3] = A[3][3] = 1.0;
  A[2][1] = v0*v0;
  A[2][2] = v0;
  A[3][1] = vp1*vp1;
  A[3][2] = vp1;
  y[1][1] = chim1;
  y[2][1] = chi0;
  y[3][1] = chip1;
  gaussj(A,3,y,1);
  /*  printf("\n\nPolished velocity shift = %6.2fkm/s\n",-y[2][1]/(2.0*y[1][1])); */
  vel = -y[2][1]/(2.0*y[1][1]); 

  free(wtmplt);
  free(ftmplt);
  free(wobs);
  free(fobs);
  free(wtrbn);
  free(ftrbn);
  free(worbn);
  free(forbn);
  free(wtrbn2);
  free(ftrbn2);
  free(wtrbn3);
  free(ftrbn3);
  free_matrix(A,1,3,1,3);
  free_matrix(a,1,3,1,1);
  free_matrix(y,1,3,1,1);
  return(vel);
}

void rebin(x,y,k,x2,y2,l1,l2,start,end,dw)
double *x,*y;
double *x2,*y2;
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

void getspectrum(infile,x,y,k)
     char infile[];
     double x[],y[];
     int *k;
{

  int l;
  char buffer[250],*tmp;
  FILE *in;

  if((in = fopen(infile,"r")) == NULL) {
    printf("\nCannot open %s\n",infile);
    exit(1);
  }
  l = 0;
  while(fgets(buffer,200,in) != NULL) {
    if(buffer[0] == '#') continue;
    tmp = strtok(buffer," ");
    x[l] = atof(tmp);
    tmp = strtok(NULL," ");
    y[l] = atof(tmp);
    l++;
  }
  *k = l;
  fclose(in);
  return;
}

int approx(x1,x2,prec)
double x1,x2,prec;
{
  double diff;

  diff = fabs(x1 - x2);
  if(diff > -1.0*prec && diff < prec) return(1);
  else return(0);
}
