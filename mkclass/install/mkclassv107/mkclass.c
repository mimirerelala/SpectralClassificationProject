#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "util.h"

typedef struct {
  char formatout[80];
  double chi2;
} results;
int N=10000;
int I,NI=1;
int sf = 0;
int Ba = 0;  /* flag denoting Ba peculiarity */
float spt2min(float pcode[]);
float spt2hyd(float pcode[]);
float spthydrogen(float spcode);
float hydrogen_profile(float spcode);
float hydrogen_index(float spcode);
float sptGlines(float spcode);
float sptmetal(float spcode);
float sptCaK(float spcode);
float sptKM(float spcode);
float lummin(float lumcode);
float lumratiomin(float lumcode);
int eqwHeI(float spcode, float lumcode);
float sptHeII(float spcode);
float heIpec(float spcode, float lumcode);
float hydrogen_profile_hot(float lumcode);
float hydD2();
float sptHeImet(float spcode);
float MgII(float spcode, float lumcode);
int lam_boo(float spcode, float lumcode,float *htype, float *metaltype);
int lam_boo2(float spcode, float lumcode,float *htype, float *metaltype);
float ratioCaIFeII();
int emission(float *wave, float *flx, int N);
FILE *ffopen();
float *X,*Y;
float Lumcode;
float Spcode;
char SPT[10];
char LUM[10];
char lib[40];
char MKLIB[300];
void code2spt(float spcode);
void code2lum(float lumcode);
float sptcode(char ispt[]);
float subclass(float sp);
void getspectrum();
void rebin();
void Oclass();
void Bclass();
void Aclass();
void FGclass();
void KMclass();
void roughtype(char *lib, char *name, float *isp, float *ilt);
void roughtype2(char *lib, char *name, float *isp, float *ilt);
int libconform(char *name);
void sp2class(float spcode, float lumcode, float metcode, 
	      float *x, float *y, int *n);
float match(float spt,float lum);
int peculiarity(float spcode, float lumcode, int *he, int *sr, int *si, 
		int *eu, int *cr, int *ba, int *ch, int *cn );
int barium(float spcode, float lumcode);
float carbon4737(float spcode, float lumcode);
float CN4215(float spcode, float lumcode);
float CHband2(float spcode, float lumcode);
void quality(float chi2, char *qual,float spt);
void templateDSO(float *xin,float *yin,float *xt,float *yt,int kin);
double spt2code(char spt[]);
double GbandCaI(double x_in, double y_in);
double CaKHe4471(double x_in, double y_in);
double TiOIndex(double x);
int DetectNN2(float *x, float *y, int n);
float Carbon(float *wave, float *flx, int N);
double max(double x, double y);
FILE *out;
FILE *Log;
float wlow=3800.0;
float whigh=4600.0;
float space=0.5;
int done = 0;
int iterate = 0;
int Nq;
char qual[50];
double shot,scool;
float spt,lum;
char libmatch[80];
char Note[20];
char to_out[80];
int To_out;
char name[80];
int flagname = 0;
results Iter[10];

int main(int argc, char *argv[])
{
  char buffer[150],output[80],ispt[10],irbn[40],corname[80];
  float *pcode,**xi;
  float isp,ilt,q,fret,ax,bx,cx,fa,fb,fc,spth,sptm,spta,fe;
  int he,sr,si,eu,cr,ba,ch,cn,pec,ni,len;
  float subh,subm;
  float CHI2;
  int i,j,l,n,nq,err;
  int flagi = 0;
  int flaglib = 0;
  int flagR = 0;
  int flag_template = 0;
  int iter = 0;
  FILE *in,*mk,*cor;
  char key1[10],key2[10],key3[10],key4[10],Lib[40],splow[6],sphigh[6],key5[10];
  char prelim[80],LOG[80],template[5];
  float *xt,*yt;
  char *p,LIB[300];
  int nt,e;
  int normal;
  int flagextra = 0;
  int nerr = 0;
  float ratio;
  int l1=0,l2=0,k=0;
  float *x0,*y0;

  he = sr = si = eu = cr = ba = ch = cn = 0;

  X = vector(0,N);
  Y = vector(0,N);
  x0 = vector(0,N);
  y0 = vector(0,N);
  xt = vector(0,N);
  yt = vector(0,N);
  pcode = vector(1,2);
  xi = matrix(1,2,1,2);

  for(i=1;i<=2;i++) {
    for(j=1;j<=2;j++) {
      if(i == j) xi[i][j] = 1.0;
      else xi[i][j] = 0.0;
    }
  }

  if(argc < 6 || argc > 7) {
    printf("Usage:  mkclass spectrum lib output log R [NI]\n\n"); 
    printf("\nInput spectrum to be classified > ");
    ggets(name);
    printf("\nEnter standards library (libr18) > ");
    ggets(lib);
    printf("\nEnter output file > ");
    ggets(output);
    printf("\nEnter details log file > ");
    ggets(LOG);
    printf("\nEnter Rough type algorithm flag (1 or 2) > ");
    ni = scanf("%d",&flagR);
    printf("\nEnter number of iterations > ");
    ni = scanf("%d",&NI);
  } else {
    strcpy(name,argv[1]);
    strcpy(lib,argv[2]);
    strcpy(output,argv[3]);
    strcpy(LOG,argv[4]);
    flagR = atoi(argv[5]);
    if(argc == 7) {
      NI = atoi(argv[6]);
      flagi = 1;
    }
  }
  Log = ffopen(LOG,"a");
  out = ffopen(output,"a");
  fprintf(Log,"MKCLASS v1.07: %s %s\n",lib,name);
  fprintf(Log,"Roughtype = %d\n",flagR);

  // Get environment variable MKLIB

  p = getenv("MKLIB");
  if(!p) {
    fprintf(Log,"MKLIB not set: using default MKLIB path: /usr/local/mkclass\n");
    strcpy(MKLIB,"/usr/local/mkclass");
  } else {
    strcpy(MKLIB,p);
    fprintf(Log,"MKLIB read from environment as %s\n",MKLIB);
  }

  // Strip trailing / from MKLIB

  l = strlen(MKLIB);
  if(MKLIB[l-1] == '/') MKLIB[l-1] = '\0';

  // Sets the names for the various output files.
  len = strcspn(name,".");
  strncpy(corname,name,strcspn(name,"."));
  corname[len] = '\0';
  strcpy(libmatch,corname);
  strcat(corname,".cor");
  strcat(libmatch,".mat");

  printf("MKCLASS v1.07 %s\n",name);

  in = fopen(name,"r");
  if(in == NULL) {
    printf("\nCannot find input spectrum file\n");
    exit(1);
  }
  fclose(in);

  pcode[1] = 16.0;
  pcode[2] = 4.5;

  strcpy(irbn,"t160l50p00.rbn");
  if(flagi == 1) {
    pcode[1] = isp = sptcode(ispt);
    sprintf(irbn,"t%2.0f0l50p00.rbn",isp);
  }

  // Get details of spectral library from mkclass.lib
  flaglib = 0;
  strcpy(LIB,MKLIB);
  strcat(LIB,"/mkclass.lib");
  mk = ffopen(LIB,"r");
  while(fscanf(mk,"%s %s",key1,Lib) != EOF) {
    ni = fscanf(mk,"%s %f %f %f",key2,&wlow,&whigh,&space);
    ni = fscanf(mk,"%s %s %s",key3,splow,sphigh);
    ni = fscanf(mk,"%s %s",key4,prelim);
    ni = fscanf(mk,"%s %s",key5,template);
    ni = fscanf(mk,"%s",buffer);
    ni = strcmp(Lib,lib);
    if(strcmp(Lib,lib) == 0) {
      flaglib = 1;
      break;
    }
  }

  // Check to see if input spectrum conforms to library specifications
  nerr = libconform(name);
  if(nerr == 1 || nerr == 2) {
    fprintf(Log,"Wavelength range does not conform to library specifications\n");
    fprintf(out,"%s Cannot classify -- see log\n",name);
    exit(1);
  }

  fprintf(Log,"prelim = %s\n",prelim);
  shot = spt2code(splow);
  scool = spt2code(sphigh);
  if(strcmp(template,"yes") == 0) flag_template = 1;
  else flag_template = 0;

  if(flaglib == 0) {
    printf("\nCannot find standards library\n");
    exit(1);
  }


  /* Now get a rough initial type.  If you are classifying a rectified
spectrum, use roughtype 1, if a flux-calibrated spectrum, roughtype 2 is
better, but experiment to see which works for your spectra */


  if(flagR == 1) roughtype(lib,name,&isp,&ilt);
  if(flagR == 2) roughtype2(lib,name,&isp,&ilt);
  fprintf(Log,"t%03.0fl%1.0f0p00.rbn\n",10*isp,ilt);
  sprintf(irbn,"t%03.0fl%1.0f0p00.rbn",10*isp,ilt);
  fprintf(Log,"Initial type = %s\n",irbn);

  getspectrum(name,x0,y0,&k);
  rebin(x0,y0,k,X,Y,&l1,&l2,wlow,whigh,space);
  Nq = (whigh-wlow)/space;

  /* After the rough type, we really need to check to see if star is normal,
We first do a preliminary check to see if the star is an emission-line star.
Then we pass the star to the normality-checking routine DetectNN2.  At the 
moment, we are doing this only if flagR is 2 */
  e = emission(X,Y,Nq);
  if(isp >= 34 && e == 0) {
    // Check to see if star is carbon star
    ratio = Carbon(X,Y,Nq);
    fprintf(Log,"Carbon ratio = %f\n",ratio);
    if(ratio >= 3.0) {
      fprintf(Log,"Classified as a Carbon Star\n");
      fprintf(out,"%s  | Carbon star\t\t|         | \\\\\n",name);
      fclose(out);
      fprintf(Log,"================================================\n\n");
      fclose(Log);
      free_vector(X,0,N);
      free_vector(Y,0,N);
      free_vector(xt,0,N);
      free_vector(yt,0,N);
      exit(1);
    } else {
      normal = 0;
      strcpy(Note," \\\\ ");
    }
  } else {
    normal = DetectNN2(X,Y,Nq);
    strcpy(Note," \\\\ ");
    if(normal == 14) {
      fprintf(out,"%s  | Unclassifiable             |         | \\\\ \n",name);
      fclose(out);
      exit(1);
    }
    if(normal == 13) strcpy(Note," ? \\\\ ");
    if(normal != 0 && normal != 7 && normal != 13) {
      fprintf(Log,"Star is not normal: code = %d\n",normal);
      if(normal == 1) fprintf(out,"%s  | DA\t\t\t\t|         | \\\\\n",name);
      if(normal == 2) fprintf(out,"%s  | DB\t\t\t\t|         | \\\\\n",name);
      if(normal == 3) fprintf(out,"%s  | DO\t\t\t\t|         | \\\\\n",name);
      if(normal == 4) fprintf(out,"%s  | DZ\t\t\t\t|         | \\\\\n",name);
      if(normal == 5) fprintf(out,"%s  | DQ\t\t\t\t|         | \\\\\n",name);
      if(normal == 8) fprintf(out,"%s  | Carbon star\t\t|         | \\\\\n",name);
      if(normal == 9) fprintf(out,"%s  | emission-line?\t\t|         | \\\\\n",name);
      if(normal == 10) fprintf(out,"%s  | WN\t\t\t\t|         | \\\\\n",name);
      if(normal == 11) fprintf(out,"%s  | Helium nova\t\t|         | \\\\\n",name);
      if(normal == 12) fprintf(out,"%s  | WC\t\t\t\t|         | \\\\\n",name);
      if(normal == 6) fprintf(out,"%s  | ??\t\t\t\t|         | \\\\\n",name);
      fclose(out);
      fprintf(Log,"================================================\n\n");
      fclose(Log);
      free_vector(X,0,N);
      free_vector(Y,0,N);
      free_vector(xt,0,N);
      free_vector(yt,0,N);
      exit(1);
    } else fprintf(Log,"Star appears normal\n");
  } 

  free_vector(x0,0,N);
  free_vector(y0,0,N);
  free_vector(X,0,N);
  free_vector(Y,0,N);

  X = vector(0,N);
  Y = vector(0,N);

  // Pass the spectrum through mkprelim
  sprintf(buffer,"mkprelim %s %s/%s/%s temp.out %7.1f %7.1f %3.1f %s",
	  name,MKLIB,lib,irbn,wlow,whigh,space,prelim);
  err = system(buffer);

  // Read in the radial-velocity corrected spectrum produced by mkprelim */
  in = ffopen("temp.out","r");

  i = 0;
  while(fscanf(in,"%f %f",&X[i],&Y[i]) != EOF) i++;
  Nq = n = i;
  Nq -= 1;
  fclose(in);

  spt = isp;
  lum = ilt;  


  /* Now correct energy distribution to initial type, provided the library
specifies a template should be used.  Rectified libraries should not apply
templates.  It is probably advantageous to use templates for flux libraries,
as this will help to correct the SED before the use of powell below. */

  // Note: changed spt limit from 39.0 to scool
  if(spt <= scool && flag_template == 1) 
    fprintf(Log,"Spectral template applied\n");


  if(spt <= scool && flag_template == 1) {
    sp2class(spt,lum,0.0,xt,yt,&nt);
    templateDSO(X,Y,xt,yt,n);
  }

  // Now refine initial type using least squares with powell

  pcode[1] = isp;
  pcode[2] = ilt;

  powell(pcode,xi,2,0.01,&iter,&fret,spt2min);  

  spt = Spcode = pcode[1];
  lum = Lumcode = pcode[2];

  if(spt >= scool) {
    fprintf(Log,"Initial spectral type %f is later than library limit.  Adjusting to within library range\n",spt);
    spt = scool - 3.0;
    fprintf(Log,"Initial spectral type adjusted to %f\n",spt);
  }
  if(lum > 5.0) lum = Lumcode = 5.0;


  code2spt(spt);
  code2lum(lum);
  fprintf(Log,"spt = %f  lum = %f\n",spt,lum);

  fprintf(out,"%s  | ",name);
  flagname = 1;
  fprintf(Log,"Initial Spectral type estimate = %s %s\n",SPT,LUM);
  done = 0;
  I = 1;
  flagextra = 0;
  while(I <= NI) {
    iterate = 0;
    // If quality is poor or fair, add one iteration
    if(I == NI && (strstr(qual,"fair") != NULL || strstr(qual,"poor") != NULL) &&
       flagextra == 0) {
      fprintf(Log,"Going for an extra iteration because of low quality\n");
      NI++;
      flagextra = 1;
    }
    /* If iterating, use the matched spectrum as a template to correct the
       energy distribution */
    lum = Lumcode;
    spt = Spcode;
    if(lum > 5.2) lum = 5.2;
    // Note: changed spt limit from 39.0 to scool
    if(I >= 2 && spt <= scool && flag_template == 1) {
      fprintf(Log,"Spectral template spt = %4.1f  lum = %3.1f\n",spt,lum);
      sp2class(spt,lum,0.0,xt,yt,&nt);
      templateDSO(X,Y,xt,yt,n);
    }
    /* Print out the flux-corrected spectrum into name.cor if spt <= 39, or
       the uncorrected if spt > 39 */
    if(I >=2 && flag_template == 1) {
      cor = fopen(corname,"w");
      for(i=0;i<n;i++) fprintf(cor,"%f %f\n",X[i],Y[i]);
      fclose(cor);
    }
 
    if(spt < 5.0 && done == 0) Oclass();
    if(spt >= 5.0 && spt < 15.0 && done == 0) Bclass();
    if(spt >= 15.0 && spt < 23.0 && done == 0) Aclass();
    if(spt >= 23.0 && spt <= 34.0 && done == 0) FGclass();
    if(spt > 34.0 && done == 0) KMclass();
    I++;
  }
  for(i=1;i<=NI;i++) fprintf(Log,"%d:  %s %3.1e\n",i,Iter[i].formatout,Iter[i].chi2);

  fprintf(Log,"================================================\n\n");

  fclose(out);
  free_vector(X,0,N);
  free_vector(Y,0,N);
  free_vector(xt,0,N);
  free_vector(yt,0,N);
  return(0);
}


int libconform(char *name)
{
  float *x0,*y0;
  int k;
  int rtn = 0;

  x0 = (float *) calloc(5000,sizeof(float));
  y0 = (float *) calloc(5000,sizeof(float));
  getspectrum(name,x0,y0,&k);

  /* Currently, this function only checks to see if the wavelength limits
     conform. */
  if(x0[0] > wlow + 100.0) rtn = 1;
  if(x0[k-1] < whigh - 100.0) rtn = 2;

  free(x0);
  free(y0);

  return(rtn);
}

// This routine supplies a rough initial spectral type for rectified spectra.
void roughtype(char *lib, char *name, float *isp, float *ilt)
{
  float *x0,*y0,*x,*y,xt,yt;
  FILE *in9;
  int l1=0,l2=0,i,j,k,l,ni;
  char tmp[80];
  char *SD[10] = {"t030l50p00.rbn","t070l50p00.rbn","t120l50p00.rbn","t190l50p00.rbn","t230l50p00.rbn","t260l50p00.rbn","t320l50p00.rbn","t360l50p00.rbn","t400l50p00.rbn","t425l50p00.rbn"};
  char *SG[10] = {"t030l30p00.rbn","t070l30p00.rbn","t120l30p00.rbn","t190l30p00.rbn","t230l30p00.rbn","t260l30p00.rbn","t320l30p00.rbn","t360l30p00.rbn","t400l30p00.rbn","t425l30p00.rbn"};
  char *SS[10] = {"t030l10p00.rbn","t070l10p00.rbn","t120l10p00.rbn","t190l10p00.rbn","t230l10p00.rbn","t260l10p00.rbn","t320l10p00.rbn","t360l10p00.rbn","t400l10p00.rbn","t425l10p00.rbn"};

  float ISP[10] = {3.0,7.0,12.0,19.0,23.0,26.0,32.0,36.0,40.0,42.5};
  float ILT[3] = {5.0,3.0,1.0};
  float chi[3][10],chimin;
  int errstat,icool,ihot;

  x0 = (float *) calloc(10000,sizeof(float));
  y0 = (float *) calloc(10000,sizeof(float));
  x = (float *) calloc(N,sizeof(float));
  y = (float *) calloc(N,sizeof(float));

  ihot = icool = 0;
  for(i=0;i<9;i++) {
    if(scool > ISP[i] && scool <= ISP[i+1]) icool = i+1;
    if(shot >= ISP[i] && shot < ISP[i+1]) ihot = i+1;
  }

  getspectrum(name,x0,y0,&k);
  rebin(x0,y0,k,x,y,&l1,&l2,wlow,whigh,space);


  for(l=0;l<3;l++) {
    for(i=ihot;i<=icool;i++) {
      chi[l][i] = 0.0;
      if(l==0) sprintf(tmp,"%s/%s/%s",MKLIB,lib,SD[i]);
      else if(l==1) sprintf(tmp,"%s/%s/%s",MKLIB,lib,SG[i]); 
      else sprintf(tmp,"%s/%s/%s",MKLIB,lib,SS[i]); 
      in9 = ffopen(tmp,"r");
      j = 0;
      while(fscanf(in9,"%f %f",&xt,&yt) != EOF) {
        if(xt < wlow+100.0 || xt > whigh-100.0 || yt == 0.0) {
	  j++;
          continue;
	}
        chi[l][i] += (yt-y[j])*(yt-y[j]);
        j++;
      }
      errstat = fclose(in9);
    }
  }
  free(x0);
  free(y0);
  free(x);
  free(y);

  j = 0;
  chimin = 1.0e+30;
  for(l=0;l<3;l++) {
    for(i=ihot;i<=icool;i++) {
      if(chimin > chi[l][i]) {
        chimin = chi[l][i];
        j = i;
	k = l;
      }
    }
  }
  *isp = ISP[j];
  *ilt = ILT[k];
  return;
}

// Module to classify O-type stars.  Still under development
void Oclass()
{
  float ax,bx,cx,fa,fb,fc;
  float q,sptheII;
  float CHI2 = 0.0;
  int J;

  fprintf(Log,"Classifying as an O-type star\n");
  
  Lumcode = lum;

  ax = spt-1;
  bx = spt;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptHeII);
  q = brent(ax,bx,cx,sptHeII,0.001,&sptheII);
  fprintf(Log,"HeII type = %f\n",sptheII);

  Spcode = spt = sptheII;

  if(Spcode >= 6.5) {
    Bclass();
    return;
  }

  code2spt(spt);
  code2lum(lum);

  CHI2 = match(spt,lum);
  quality(CHI2,qual,spt);
  Iter[I].chi2 = CHI2;
  if(I<=NI) {
    sprintf(to_out,"%s %s",SPT,LUM);
    To_out = strlen(to_out);
    if(To_out <= 9) sprintf(Iter[I].formatout,"%s %s \t\t\t%s",SPT,LUM,qual);
    else if(To_out > 9 && To_out <= 17) sprintf(Iter[I].formatout,"%s %s \t\t%s",SPT,LUM,qual);
    else sprintf(Iter[I].formatout,"%s %s \t%s",SPT,LUM,qual);
  }
  if(I == NI) {
    J = findbest(Iter,NI);
    fprintf(out,"%s\n",Iter[J].formatout);
    fprintf(Log,"\nBest iteration: I = %d\n",J);
    done = 1;
  }
  fprintf(Log,"%d: %s %s  %f\n",I,SPT,LUM,CHI2);

  return;
}


// Module to classify B-type stars.  Still requires some development
void Bclass()
{
  extern int sf;
  float *qcode,**xi;
  float ax,bx,cx,fa,fb,fc;
  float q,spth,sptm,sptK,sptheII;
  int he,sr,si,eu,cr,ba,ch,cn,pec;
  float lum;
  int i,j,iter,J;
  float fret,diff,mg;
  int flagHe = 0;
  char PEC[80];
  float CHI2 = 0.0;
  float D2,htype,metaltype;
  int flaglb = 0;
  char kSPT[10],hSPT[10],mSPT[10];

  fprintf(Log,"Classifying as a B-type star\n");
  iterate++;

  /* First, check to make certain that this is not really an A-type star
     by looking at the Ca K-line */

  ax = spt-1;
  bx = spt;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptCaK);
  q = brent(ax,bx,cx,sptCaK,0.001,&sptK);
  fprintf(Log,"K-line type = %f\n",sptK);

  if(sptK >= 15.0) {
    spt = sptK;
    Aclass();
    return;
  }

  // Lets make certain this is not a DA star with very wide hydrogen lines
  D2 = hydD2();
  fprintf(Log,"D2 = %4.1f\n",D2);
  if(D2 >= 30.0) {
    fprintf(out," DA\t\t\t\t|         |\\\\\n");
    fprintf(Log,"D2 = %4.1f  DA\n",D2);
    I = NI;
    return;
  }

  /* Let us now check that this is not an O-type star by examining an
     He II line */

  ax = spt-1;
  bx = spt;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptHeII);
  q = brent(ax,bx,cx,sptHeII,0.001,&sptheII);
  fprintf(Log,"HeII type = %f\n",sptheII);

  if(sptheII <= 6.5) {
    spt = sptheII;
    Oclass();
    return;
  }

  /* We now determine an approximate temperature type based on Helium I
     strengths and metallic-line ratios (more appropriate for early-B type 
     stars */

  ax = spt-1;
  bx = spt;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptHeImet);
  q = brent(ax,bx,cx,sptHeImet,0.001,&spt);

  fprintf(Log,"Helium/metal spectral type = %f\n",spt);

  Spcode = spt;

  /* Now improve the luminosity type */

  ax = lum-0.5;
  bx = lum;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,lumratiomin);
  q = brent(ax,bx,cx,lumratiomin,0.001,&lum);
  if(lum > 5.2) lum = 5.2;
  fprintf(Log,"Luminosity type = %f\n",lum);

  /* Polish the temperature type */ 

  Lumcode = lum;

  ax = spt-1;
  bx = spt;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptHeImet);
  q = brent(ax,bx,cx,sptHeImet,0.001,&spt);

  fprintf(Log,"Polished helium/metal spectral type = %f\n",spt);

  /* Now check for helium peculiarities */

  diff = heIpec(spt,lum);

  fprintf(Log,"Helium strength = %f\n",diff);

  if(diff <= -0.1) strcpy(PEC,"Helium weak");
  if(diff >= 0.1) strcpy(PEC,"Helium strong");

  /* Check for any other obvious peculiarities */

  pec = peculiarity(spt,lum,&he,&sr,&si,&eu,&cr,&ba,&ch,&cn);

  /* Check to see if the star is a Lambda Boo star by looking at the
       Mg II 4481 line */
  if(spt >= 13.5 && pec == 0) {
    mg = MgII(spt,lum);
    fprintf(Log,"mg = %f\n",mg);
    if(mg <= -0.005) fprintf(Log,"This looks like it could be a Lambda Boo star\n");
    if(mg <= -0.005) flaglb = lam_boo(spt,lum,&htype,&metaltype);
    if(flaglb == 1) {
      fprintf(Log,"Yes, it is a Lambda Boo star\n");
      code2spt(htype);
      strcpy(hSPT,SPT);
      code2spt(metaltype);
      strcpy(mSPT,SPT);
      fprintf(Log,"%s m%s V Lam Boo\n",hSPT,mSPT);
      sprintf(to_out,"%s m%s V Lam Boo  ",hSPT,mSPT);
      sprintf(Iter[I].formatout,"%s m%s V Lam Boo  ",hSPT,mSPT);
      Iter[I].chi2 = 1.0;
      if(I == NI) {
	To_out = strlen(to_out);
	if(To_out <= 9) fprintf(out,"%s m%s V Lam Boo  \t\t\t|          | %s\n",hSPT,mSPT,Note);
	else if(To_out > 9 && To_out <= 17) fprintf(out,"%s m%s V Lam Boo  \t\t|          | %s\n",hSPT,mSPT,Note);
	else fprintf(out,"%s m%s V Lam Boo  \t|          | %s\n",hSPT,mSPT,Note);
      }
      return;
    } else if(mg <= -0.005)  
      fprintf(Log,"No, it is not a Lambda Boo\n");
  }

  pec = peculiarity(spt,lum,&he,&sr,&si,&eu,&cr,&ba,&ch,&cn);

  strcpy(PEC," ");
  if(pec == 1) {
     strcpy(PEC," ");
     if(sr == 1) strcat(PEC,"(Sr)");
     if(sr == 2) strcat(PEC,"Sr");
     if(si == 1) strcat(PEC,"Si");
     if(eu == 1) strcat(PEC,"Eu");
     if(cr == 1) strcat(PEC,"Cr");
  }

  code2spt(spt);
  code2lum(lum);

  CHI2 = match(spt,lum);
  quality(CHI2,qual,spt);
  Iter[I].chi2 = CHI2;
  if(I<=NI) {
    sprintf(to_out,"%s %s %s",SPT,LUM,PEC);
    To_out = strlen(to_out);
    if(To_out <= 9) sprintf(Iter[I].formatout,"%s %s %s \t\t\t%s",SPT,LUM,PEC,qual);
    else if(To_out > 9 && To_out <= 17) sprintf(Iter[I].formatout,"%s %s %s \t\t%s",SPT,LUM,PEC,qual);
    else sprintf(Iter[I].formatout,"%s %s %s \t%s",SPT,LUM,PEC,qual);
  }
  fprintf(Log,"%d: %s %s %s  %f\n",I,SPT,LUM,PEC,CHI2);
  if(I==NI) {
    // Repress hypergiants, as the code is not competent to classify them
    if(lum < 0.0) fprintf(out,"Unclassifiable             |         | \\\\ \n");
    else {   
      J = findbest(Iter,NI); 
      fprintf(out,"%s\n",Iter[J].formatout);
      fprintf(Log,"\nBest iteration: I = %d\n",J);
      done = 1;
    }
  }

  return;
}


//Module to classify A-type stars
void Aclass()
{

  float ax,bx,cx,fa,fb,fc;
  float q,spth,sptm,sptK;
  float chi2,fret;
  float *qcode,**xi,lumh;
  int he,sr,si,eu,cr,ba,ch,cn,pec;
  int Amflag = 0;
  int flag_early = 0;
  int i,j,iter;
  char kSPT[10],hSPT[10],mSPT[10],tmpout[80];
  extern int sf;
  char PEC[20];
  float CHI2 = 0.0;
  float D2,mg;
  int flaglb = 0;
  float htype,metaltype;
  int flagvsini = 0;
  int J;

  sf = 2;
  fprintf(Log,"\nClassifying this star as an A-type star\n");
  iterate++;

  he = sr = si = eu = cr = ba = ch = cn = 0;

  // Lets make certain this is not a DA star with very broad hydrogen lines
  D2 = hydD2();
  fprintf(Log,"D2 = %4.1f\n",D2);
  if(D2 >= 30.0) {
    fprintf(out,"DA \t\t\t|          | \\\\ \n");
    fprintf(Log,"D2 = %4.1f  DA\n",D2);
    I = NI;
    return;
  }

  /* We first begin by improving the luminosity type assuming the
preliminary spectral type.  This can help in the case of stars in which
the hydrogen-line profiles have been poorly rectified */

  ax = lum-0.5;
  bx = lum;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,lumratiomin);
  q = brent(ax,bx,cx,lumratiomin,0.001,&lum);
  if(lum > 5.2) lum = 5.2;
  fprintf(Log,"lum = %f\n",lum);

  /* We also check at this point for obvious peculiarities */

  pec = peculiarity(spt,lum,&he,&sr,&si,&eu,&cr,&ba,&ch,&cn);


  // Early A-type stars
  To_out = 0;
  if(spt <= 19.5) {
    flag_early = 1;
    // Determine the K-line type
    ax = spt-1;
    bx = spt;
    mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptCaK);
    q = brent(ax,bx,cx,sptCaK,0.001,&sptK);
    fprintf(Log,"K-line type = %f\n",sptK);

    ax = sptK-1;
    bx = sptK;
    mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptmetal);
    q = brent(ax,bx,cx,sptmetal,0.001,&sptm);
    fprintf(Log,"Metallic-line type = %f\n",sptm);

    /* if the resolution is poor, the metallic-line type can be faulty, so
       if sptm < sptK, set sptm = sptK */
    if(sptm < sptK-2) {
       sptm = sptK;
       fprintf(Log,"Correcting metallic-line type to K-line type\n");
    }

    /* Check to see if the star is a metal-weak A-type star by looking at the
       Mg II 4481 line */
    mg = MgII(spt,lum);
    fprintf(Log,"mg = %f\n",mg);
    if(mg <= -0.005 && pec == 0) fprintf(Log,"This looks like it could be a metal-weak star\n");
    flaglb = 0;
    if(mg <= -0.005 && pec == 0) flaglb = lam_boo2(spt,lum,&htype,&metaltype);
    if(flaglb == 1) {
      fprintf(Log,"It appears to be a Lambda Boo star\n");
      code2spt(htype);
      strcpy(hSPT,SPT);
      code2spt(metaltype);
      strcpy(mSPT,SPT);
      fprintf(Log,"%s m%s V Lam Boo\n",hSPT,mSPT);
      sprintf(to_out,"%s m%s V Lam Boo",hSPT,mSPT);
      sprintf(Iter[I].formatout,"%s m%s V Lam Boo",hSPT,mSPT);
      Iter[I].chi2 = 1.0; //fake chi2 for peculiar star
      /* If the final solution is lam boo, print that out without selecting
	 the "best" solution on the basis of chi2 */
      if(I == NI) {
	To_out = strlen(to_out);
        if(To_out <= 9) fprintf(out,"%s m%s V Lam Boo \t\t\t|         |\\\\\n",hSPT,mSPT);
	else if(To_out > 9 && To_out <= 17) fprintf(out,"%s m%s V Lam Boo \t\t|         |\\\\\n",hSPT,mSPT);
	else fprintf(out,"%s m%s V Lam Boo \t|         |\\\\\n",hSPT,mSPT);
      }
      return;
    } else if(mg <= -0.005 && flaglb == 0)  
      fprintf(Log,"No, it is not metal-weak; it may be a luminous A-type star\n");
    if(flaglb == 2) {
      fprintf(Log,"It appears to be a metal-weak star\n");
      code2spt(htype);
      strcpy(hSPT,SPT);
      code2spt(metaltype);
      strcpy(mSPT,SPT);
      code2lum(Lumcode);
      fprintf(Log,"%s m%s %s metal-weak\n",hSPT,mSPT,LUM);
      sprintf(to_out,"%s m%s %s metal-weak",hSPT,mSPT,LUM);
      sprintf(Iter[I].formatout,"%s m%s %s metal-weak",hSPT,mSPT,LUM);
      Iter[I].chi2 = 1.0; //fake chi2 for peculiar star
      /* If the final solution is metal-weak, print that out without selecting
	 the "best" solution on the basis of chi2 */
      if(I == NI) {
	To_out = strlen(to_out);
	if(To_out <= 9) fprintf(out,"%s m%s %s metal-weak \t\t\t|         |\\\\\n",hSPT,mSPT,LUM);
	else if(To_out > 9 && To_out <= 17) fprintf(out,"%s m%s %s metal-weak \t\t|         |\\\\\n",hSPT,mSPT,LUM);
	else fprintf(out,"%s m%s %s metal-weak \t|         |\\\\\n",hSPT,mSPT,LUM);
      }
      return;
    } else if(mg <= -0.005 && flaglb == 0)  
      fprintf(Log,"No, it is not metal-weak\n");
    if(flaglb == 3) {
      fprintf(Log,"This star is probably rapidly rotating\n");
      flagvsini = 1;
    }

    Amflag = 0;
    if((sptm - sptK) > 2.0) Amflag = 1;
    spth = (sptm + sptK)/2.0;  /* because it is difficult to determine
                                  spectral type from the hydrogen lines
                                  alone if star is early A-type */
  } else {
  /* Late A-type star: Determine H-line type, assuming preliminary 
     luminosity type */

    /* Check to see if the star is a Lambda Boo star by looking at the
       Mg II 4481 line */
    mg = MgII(spt,lum);
    fprintf(Log,"mg = %f\n",mg);
    if(mg <= -0.005 && pec == 0) fprintf(Log,"This looks like it could be a metal-weak star\n");
    if(mg <= -0.005 && pec == 0) {
      flaglb = lam_boo2(spt,lum,&htype,&metaltype);
      fprintf(Log,"flaglb = %d\n",flaglb);
    }
    if(flaglb == 1) {
      fprintf(Log,"It appears to be a Lambda Boo star\n");
      code2spt(htype);
      strcpy(hSPT,SPT);
      code2spt(metaltype);
      strcpy(mSPT,SPT);
      fprintf(Log,"%s m%s V Lam Boo\n",hSPT,mSPT);
      sprintf(to_out,"%s m%s V Lam Boo",hSPT,mSPT);
      Iter[I].chi2 = 1.0; //fake chi2 for peculiar star
      sprintf(Iter[I].formatout,"%s m%s V Lam Boo",hSPT,mSPT);
      if(I == NI) {
	To_out = strlen(to_out);
	if(To_out <= 9) fprintf(out,"%s m%s V Lam Boo \t\t\t|         |\\\\\n",hSPT,mSPT);
	else if(To_out > 9 && To_out <= 17) fprintf(out,"%s m%s V Lam Boo \t\t|         |\\\\\n",hSPT,mSPT);
	else fprintf(out,"%s m%s V Lam Boo \t|         |\\\\\n",hSPT,mSPT);
      }
      return;
    } else if(mg <= -0.005)  
      fprintf(Log,"No, it is not a Lambda Boo\n");
    if(flaglb == 2) {
      fprintf(Log,"It appears to be a metal-weak star\n");
      code2spt(htype);
      strcpy(hSPT,SPT);
      code2spt(metaltype);
      strcpy(mSPT,SPT);
      code2lum(Lumcode);
      fprintf(Log,"%s m%s %s metal-weak\n",hSPT,mSPT,LUM);
      sprintf(to_out,"%s m%s %s metal-weak",hSPT,mSPT,LUM);
      Iter[I].chi2 = 1.0;//fake chi2 for peculiar star
      sprintf(Iter[I].formatout,"%s m%s %s metal-weak",hSPT,mSPT,LUM); 
      if(I == NI) {
	To_out = strlen(to_out);
	if(To_out <= 9) fprintf(out,"%s m%s %s metal-weak\t\t\t|         | %s\n",hSPT,mSPT,LUM,Note);
	else if(To_out > 9 && To_out <= 17) fprintf(out,"%s m%s %s metal-weak\t\t|         | %s\n",hSPT,mSPT,LUM,Note);
	else fprintf(out,"%s m%s %s metal-weak\t|         | %s\n",hSPT,mSPT,LUM,Note);
      }
      return;
    } else if(mg <= -0.005)  
      fprintf(Log,"No, it is not a metal-weak star\n");
    if(flaglb == 3) {
      flagvsini = 1;
      fprintf(Log,"This is probably a rapidly rotating star\n");
    }
    if(flaglb == 0) fprintf(Log,"This is not a metal-weak star\n");


    flag_early = 0;
    ax = spt-1;
    bx = spt;
    mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,hydrogen_profile);
    q = brent(ax,bx,cx,hydrogen_profile,0.001,&spth);
    fprintf(Log,"Hydrogen-line type = %f\n",spth);

    // Determine the K-line type
    ax = spt-1;
    bx = spt;
    mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptCaK);
    q = brent(ax,bx,cx,sptCaK,0.001,&sptK);
    fprintf(Log,"K-line type = %f\n",sptK);

    // Pass it back to FGclass if it is F-type but not an Am star

    if(sptK > 22.8 && spth > 22.8 && iterate < 5) {
      spt = spth;
      FGclass();
      return;
    }

    // Determine metallic-line type, assuming preliminary luminosity type

    ax = spt-1;
    bx = spt;
    mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptmetal);
    q = brent(ax,bx,cx,sptmetal,0.001,&sptm);
    fprintf(Log,"Metallic-line type = %f\n",sptm);

/* if the resolution is poor, the metallic-line type can be faulty, so
       if sptm < sptK, set sptm = sptK */
    if(sptm < sptK-2) {
      sptm = sptK;
      fprintf(Log,"Correcting metallic-line type to K-line type\n");
    }

  /* Determine K-line type, assuming preliminary luminosity type */
    ax = spt-1;
    bx = spt;
    mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptCaK);
    q = brent(ax,bx,cx,sptCaK,0.001,&sptK); 
    fprintf(Log,"K-line type = %f\n",sptK);
    fprintf(Log,"Metallic-line type = %f\n",sptm);
    fprintf(Log,"Hydrogen-line type = %f\n",spth);
    Amflag = 0;
    if((sptm-sptK) >= 2.0) Amflag = 1;
    else if((sptm+sptK)/2 < spth-2) Amflag = 2;
    fprintf(Log,"Amflag = %d\n",Amflag);
  }


  if(Amflag == 0) {
    /* Take a mean temperature type.  For early A-type stars, overweight the metallic and K-line types */
    if(spt <= 19.5) Spcode = spt = (sptm + 3.0*sptK)/4.0;
    else Spcode = spt = (spth + sptm + sptK)/3.0;
    fprintf(Log,"Merged spt = %f\n",spt);

    if(spt < 15.0 && iterate < 5) {
      Bclass();
      return;
    }

    // Now improve the luminosity type
     ax = lum-0.5;
     bx = lum;
     mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,lumratiomin);
     q = brent(ax,bx,cx,lumratiomin,0.001,&lum);
     if(lum > 5.2) lum = 5.2;
     fprintf(Log,"lum = %f\n",lum);



     pec = peculiarity(spt,lum,&he,&sr,&si,&eu,&cr,&ba,&ch,&cn);
     if(pec == 1) {
       strcpy(PEC," ");
       if(sr == 1) strcat(PEC,"(Sr)");
       if(sr == 2) strcat(PEC,"Sr");
       if(si == 1) strcat(PEC,"Si");
       if(eu == 1) strcat(PEC,"Eu");
       if(cr == 1) strcat(PEC,"Cr");
     }

     Lumcode = lum;

     // Introduce chi2 selection code starting here

     code2spt(spt);
     code2lum(lum);
     CHI2 = match(spt,lum);
     quality(CHI2,qual,spt);
     Iter[I].chi2 = CHI2;
     if(pec == 0) {
       if(I<=NI) {
         if(flagvsini == 1) strcat(LUM,"n");
         sprintf(to_out,"%s %s",SPT,LUM);
         To_out = strlen(to_out);
	 if(To_out <= 9) sprintf(Iter[I].formatout,"%s %s \t\t\t%s",SPT,LUM,qual);
         else if(To_out > 9 && To_out <= 17) sprintf(Iter[I].formatout,"%s %s \t\t%s",SPT,LUM,qual);
	 else sprintf(Iter[I].formatout,"%s %s \t%s",SPT,LUM,qual);
       }
       if(I==NI) {
	 // Suppress hypergiants as code is not competent
	 if(lum < 0.0) fprintf(out,"Unclassifiable             |         | \\\\ \n");
	 else {
	   J = findbest(Iter,NI);
	   fprintf(out,"%s\n",Iter[J].formatout);
           fprintf(Log,"\nBest iteration: I = %d\n",J);
	   done = 1;
	 }
       }
       fprintf(Log,"%d: %s %s %7.1e\n\n",I,SPT,LUM,CHI2);
       } else {
	 if(I<=NI) {
           if(flagvsini == 1) strcat(LUM,"n");
           sprintf(to_out,"%s %s %s",SPT,LUM,PEC);
	   To_out = strlen(to_out);
           if(To_out <= 9) sprintf(Iter[I].formatout,"%s %s %s \t\t\t%s",SPT,LUM,PEC,qual);
	   else if(To_out > 9 && To_out <= 17) sprintf(Iter[I].formatout,"%s %s %s \t\t%s",SPT,LUM,PEC,qual);
	   else sprintf(Iter[I].formatout,"%s %s %s \t%s",SPT,LUM,PEC,qual);
	 }
       if(I==NI) {
	 if(lum < 0.0) fprintf(out,"Unclassifiable             |         | \\\\ \n");
	 else {
	   J = findbest(Iter,NI);
	   fprintf(out,"%s\n",Iter[J].formatout);
           fprintf(Log,"\nBest iteration: I = %d\n",J);
	   done = 1;
	 }
       }
       fprintf(Log,"%d: %s %s %s %7.1e\n\n",I,SPT,LUM,PEC,CHI2);
       }
  } else if(Amflag == 1) {
    code2spt(sptK);
    strcpy(kSPT,SPT);
    code2spt(spth);
    strcpy(hSPT,SPT);
    code2spt(sptm);
    strcpy(mSPT,SPT);
    spt = spth;
    strcpy(PEC," ");
    pec = peculiarity(spt,lum,&he,&sr,&si,&eu,&cr,&ba,&ch,&cn);
    if(pec == 1) {
       strcpy(PEC," ");
       if(sr == 1) strcat(PEC,"Sr");
       if(si == 1) strcat(PEC,"Si");
       if(eu == 1) strcat(PEC,"Eu");
       if(cr == 1) strcat(PEC,"Cr");
    }
    fprintf(Log,"%d: k%sh%sm%s %s\n",I,kSPT,hSPT,mSPT,PEC);
    Iter[I].chi2 = 1.0;
    sprintf(Iter[I].formatout,"k%sh%sm%s %s",kSPT,hSPT,mSPT,PEC);
    if(I==NI) {
      sprintf(to_out,"k%sh%sm%s %s",kSPT,hSPT,mSPT,PEC);
      To_out = strlen(to_out);
      if(To_out <= 9) fprintf(out,"k%sh%sm%s %s\t\t\t|         |%s\n",kSPT,hSPT,mSPT,PEC,Note);
      else if(To_out > 9 && To_out <= 17) fprintf(out,"k%sh%sm%s %s\t\t|         |%s\n",kSPT,hSPT,mSPT,PEC,Note);
      else fprintf(out,"k%sh%sm%s %s\t|         |%s\n",kSPT,hSPT,mSPT,PEC,Note); 
      done = 1;
    }
  } else if(Amflag == 2) {
    code2spt(sptm);
    strcpy(mSPT,SPT);
    code2spt(spth);
    strcpy(hSPT,SPT);
    spt = spth;
    fprintf(Log,"%d: %s m%s  metal weak\n",I,hSPT,mSPT);
    Iter[I].chi2 = 1.0;
    sprintf(Iter[I].formatout,"%s m%s metal weak",hSPT,mSPT);
    if(I==NI) {
      sprintf(to_out,"%s m%s metal weak",hSPT,mSPT);
      To_out = strlen(to_out);
      if(To_out <= 9) fprintf(out,"%s m%s metal weak\t\t\t|         | %s\n",hSPT,mSPT,Note);
      else if(To_out > 9 && To_out <= 17) fprintf(out,"%s m%s metal weak\t\t|         | %s\n",hSPT,mSPT,Note);
      else fprintf(out,"%s m%s metal weak\t|         | %s\n",hSPT,mSPT,Note);
      done = 1;
    }
  } else {
    if(I==NI) {
      fprintf(out,"Undetermined peculiar A star\t|         |\\\\\n");
      done = 1;
    }
    fprintf(Log,"%d: Undetermined peculiar A-type star\n",I);
  }

  return;
}

// Module for classifying F and G-type stars
void FGclass()
{

  float ax,bx,cx,fa,fb,fc,lumorig,sptorig,lumBa;
  float q,spth,sptm,sptK,spthm,subh,subm,fe,spta;
  int Amflag = 0;
  int pec = 0;
  int flagmetal = 0;
  int flagrich = 0;
  int flaglb = 0;
  int he,sr,si,eu,ba,cr,ch,cn;
  float htype,metaltype,mg;
  extern int sf;
  char PEC[20],hSPT[20],mSPT[20],tmpout[80];
  float CHI2 = 0.0;
  float C2,CN,CH;
  int flagfe = 0;
  int J;

  strcpy(PEC," ");

  fprintf(Log,"\nClassifying this star as an F/G star\n");
  iterate++;

  he = sr = si = eu = cr = ba = ch = cn = 0;

  sf = 3;

  lumorig = lum;
  sptorig = spt;

  // Determine the hydrogen-line type
  ax = spt-1;
  bx = spt;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,hydrogen_index);
  q = brent(ax,bx,cx,hydrogen_index,0.001,&spth);
  code2spt(spth);
  fprintf(Log,"Hydrogen-line type = %s %f\n",SPT,spth);
 
  // If it is an A-type star, pass it back to the A module
  if(spth <= 22.8 && iterate < 5) {
    spt = spth;
    Aclass();
    return;
  }

  // If it is a K-type star, pass it to the KM module
  if(spth >= 34.0 && iterate < 5) {
    if(spth >= 36.0) spt = 36.0;
    else spt = spth;
    KMclass();
    return;
  }

  // Determine the hydrogen to metal spectral type
  spthm = spth;
  if(spth >= 30.0) { 
    ax = spth-1;
    bx = spth;
    mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptGlines);
    q = brent(ax,bx,cx,sptGlines,0.001,&spthm);
    code2spt(spthm);
    fprintf(Log,"Hydrogen to metallic line type = %s %f\n",SPT,spthm);
  }
 
  ax = spth-1;
  bx = spth;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptmetal);
  q = brent(ax,bx,cx,sptmetal,0.001,&sptm);
  code2spt(sptm);
  fprintf(Log,"Metallic-line type = %s %f\n\n",SPT,sptm);

  /* If hydrogen-line type and metallic line type differ too widely, set
     a peculiarity flag for future reference */
  if(fabs(sptm - spth) > 2.0 || fabs(spthm - spth) > 3.0) flagmetal = 1;
  if(sptm - spth > 1.0) flagrich = 1;
  if(sptm - spth > 2.0) flagrich = 2;
  if(flagrich != 0) fprintf(Log,"Star metal rich?\n");  

  fprintf(Log,"flagmetal = %d\n",flagmetal);

  /* If things look normal, try to determine a luminosity type */
  if(flagmetal == 0) {
    spt = (3.0*spth + 2.0*spthm + sptm)/6.0;
    fprintf(Log,"Normal star: spt = %f\n",spt);
    Spcode = spt;
    Ba = 0; /* So that all luminosity criteria are considered */
    ax = lum - 0.5;
    bx = lum;
    mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,lumratiomin);
    q = brent(ax,bx,cx,lumratiomin,0.001,&lum);
    fprintf(Log,"lum = %f\n",lum);
    if(lum > 5.0) lum = 5.0;
    code2spt(spt);
    code2lum(lum);
  }


  // Now, begin check for peculiarities

  // If Hydrogen-line type is earlier than F1, check to see if Lambda Boo
  if(spth <= 23.5 && flagmetal == 1) {
    mg = MgII(spt,lum);
    fprintf(Log,"mg = %f\n",mg);
    if(mg <= -0.005) fprintf(Log,"This looks like it could be a metal-weak star\n");
    if(mg <= -0.005) flaglb = lam_boo2(spt,lum,&htype,&metaltype);
    if(flaglb == 1) {
      fprintf(Log,"It appears to be a Lambda Boo star\n");
      code2spt(htype);
      strcpy(hSPT,SPT);
      code2spt(metaltype);
      strcpy(mSPT,SPT);
      fprintf(Log,"%s m%s V Lam Boo\n",hSPT,mSPT);
      if(I == NI) fprintf(out,"%s m%s V Lam Boo\n",hSPT,mSPT);
      return;
    }
  }

  // If Hydrogen-line type is earlier than F5, check for Am characteristics

  if(spth <= 27.0) {
    ax = spth-1;
    bx = spth;
    mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptCaK);
    q = brent(ax,bx,cx,sptCaK,0.001,&sptK);

    fprintf(Log,"sptK = %f sptm = %f\n",sptK,sptm);

    /* If K-line type is earlier than metallic-line type, and the metallic type
       is equal to or later than the hydrogen-line type, it is probably an Am 
       star, so pass to Aclass function */
    if((sptm-sptK) >= 2.0 && (sptm-spth) >= -0.5 && iterate < 5) {
      spt = spth;
      Aclass();
      return;
    }
    // Also check for Ap-type peculiarities

    pec = peculiarity(spth,lum,&he,&sr,&si,&eu,&cr,&ba,&ch,&cn);
    if(pec == 1) {
       strcpy(PEC," ");
       if(sr == 1) strcat(PEC,"(Sr)");
       if(sr == 2) strcat(PEC,"Sr");
       if(si == 1) strcat(PEC,"Si");
       // if(eu == 1) strcat(PEC,"Eu");
       if(cr == 1) strcat(PEC,"Cr");
    }
  }

  /* Determine an approximate luminosity class, leaving out Sr II lines.
     Then, check for a barium peculiarity.  */

  Spcode = spth;
  Ba = 1; /* So that Sr II is left out of luminosity classification */
  lumBa = lum;
  ax = lumBa - 0.5;
  bx = lumBa;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,lumratiomin);
  q = brent(ax,bx,cx,lumratiomin,0.001,&lumBa);
  fprintf(Log,"lumBa = %f\n",lumBa);
  if(lumBa > 5.0) lumBa = 5.0;

  /* The barium function checks for Sr peculiarity if earlier than G5,
     for a Ba peculiarity if G5 or later */
  Ba = barium(spt,lumBa);  /* Check for Ba (Sr) peculiarity */

  // Check for carbon peculiarities
  fprintf(Log,"C2 ratio = %f\n",carbon4737(spt,lum));
  fprintf(Log,"CN ratio = %f\n",CN4215(spt,lum));
  C2 = CN = CH = 0.0;
  if(lum < 4.5  && spt < 39.0 ) C2 = carbon4737(spt,lum);
  if(lum < 4.5 && spt < 39.0) CN = CN4215(spt,lum);
  if(spt >= 32 && spt < 39.0 && lum < 4.5) {  
    CH = CHband2(spt,lum);
    fprintf(Log,"CH ratio = %f\n",CH);
  }
  if(CN > 0.025 && CN < 0.05) { 
    strcat(PEC,"");
    pec = 1;
  }
  if(CN >= 0.05 && CN < 0.075) {
    strcat(PEC," CN1");
    pec = 1;
  }
  if(CN >= 0.075) {
    strcat(PEC," CN2");
    pec = 1;
  }
  if(CN < -0.025 && CN > -0.05) {
    strcat(PEC,"");
    pec = 1;
  }
  if(CN <= -0.05) {
    strcat(PEC," CN-1");
    pec = 1;
  }
  if(spt >= 32) {
    if(CH > 0.075 && CH < 0.150) {
       strcat(PEC,"");
       pec = 1;
    }
    if(CH >= 0.150 && CH < 0.200) {
       strcat(PEC," CH1");
       pec = 1;
    }
    if(CH > 0.200) {
       strcat(PEC," CH2");
       pec = 1;
    }
    if(CH <= -0.075 && CH > -0.150) {
      strcat(PEC,"");
      pec = 1;
    }
    if(CH <= -0.150 && CH > -0.200) {
      strcat(PEC," CH-1");
      pec = 1;
    }
    if(CH <= -0.200) {
      strcat(PEC," CH-2");
      pec = 1;
    }
  }

  if(Ba == 1 && spt < 32) {
    pec = 1;
    lum = lumBa;
    strcat(PEC," Sr");
  }
  if(Ba == 1 && spt >= 32) {
    pec = 1;
    lum = lumBa;
    strcat(PEC," Ba");
  }

  code2spt(spt);
  code2lum(lum);

  Spcode = spt;
  Lumcode = lum;

  // If Ba is not peculiar, then accept the earlier luminosity class

  subh = subclass(spth);
  subm = subclass(sptm);

  fe = 0.0;
  if(subh - subm > 2.0) fe = -0.13*(subh-subm) - 0.26;
  if(subm - subh > 2.0) fe = 0.25*(subm - subh);
  if(subm - subh > 1.0 && subm - subh < 2.0) flagrich = 1; 

  To_out = 0;
  if(fabs((subh - subm) < 2.0 || fabs(fe) < 0.50)) {
    CHI2 = match(spt,lum);
    quality(CHI2,qual,spt);
    if(I<=NI) {
      sprintf(Iter[I].formatout,"%s %s",SPT,LUM);
      Iter[I].chi2 = CHI2;
      sprintf(to_out,"%s %s",SPT,LUM);
      To_out += strlen(to_out);
    }
    fprintf(Log,"%d: %s %s",I,SPT,LUM);
    if(flagrich == 1) fprintf(Log," ((metal-rich))");
    if(pec == 0) {
      if(I<=NI) {
	sprintf(tmpout," \t\t\t%s",qual);
        strcat(Iter[I].formatout,tmpout);
      }
      if(I==NI) {
	J = findbest(Iter,NI);
	fprintf(out,"%s\n",Iter[I].formatout);
        fprintf(Log,"\nBest iteration: I = %d\n",J);
	done = 1;
      }
      fprintf(Log," %7.1e\n",CHI2);
    }
    else {
      if(I<=NI) {
	To_out += strlen(PEC);
	if(To_out <= 9) {
          sprintf(tmpout,"%s \t\t\t%s",PEC,qual);
	  strcat(Iter[I].formatout,tmpout);
	}
	else if(To_out > 9 && To_out <= 17) {
          sprintf(tmpout,"%s \t\t%s",PEC,qual);
          strcat(Iter[I].formatout,tmpout);
	}
        else {
          sprintf(tmpout,"%s \t%s",PEC,qual);
	  strcat(Iter[I].formatout,tmpout);
	}
      }
      if(I==NI) {
	J = findbest(Iter,NI);
	fprintf(out,"%s\n",Iter[J].formatout);
      fprintf(Log,"\nBest iteration: I = %d\n",J);
	done = 1;
      }
      fprintf(Log,"%s %7.1e\n",PEC,CHI2);
    }
  } else {
    // metal-weak or metal-rich stars
    code2spt(spth);
    Spcode = spt = spth;
    fprintf(Log,"Metal weak or rich: spt = %f\n",spt);
    ax = 2.5;
    bx = 3.0;
    mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,lumratiomin);
    q = brent(ax,bx,cx,lumratiomin,0.001,&lum);
    if(lum > 5.0) {
      lum = Lumcode = 5.0;
      ax = spth-1;
      bx = spth;
      mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptmetal);
      q = brent(ax,bx,cx,sptmetal,0.001,&sptm); 
      code2spt(sptm);
      subh = subclass(spth);
      subm = subclass(sptm);
      fe = -0.13*(subh-subm) - 0.26;
    }
    Lumcode = lum;
    code2spt(spth);
    code2lum(lum);
    CHI2 = match(spth,lum);
    Iter[I].chi2 = CHI2;
    quality(CHI2,qual,spth);

    flagfe = 0;
    if(fabs(fe) > 0.50) {
      if(I<=NI) {
	sprintf(Iter[I].formatout,"%s %s Fe%+4.1f",SPT,LUM,fe);
        sprintf(to_out,"%s %s Fe%+4.1f",SPT,LUM,fe);
	To_out = strlen(to_out);
      }
      if(I==NI) {
	J = findbest(Iter,NI);
	done = 1;
      }
      fprintf(Log,"%d:  %s %s Fe%+4.1f",I,SPT,LUM,fe);
      flagfe = 1;
    }
    else {
      if(I<=NI) {
	sprintf(Iter[I].formatout,"%s %s ",SPT,LUM);
      }
      if(I==NI) {
	J = findbest(Iter,NI);
	sprintf(to_out,"%s %s ",SPT,LUM);
	To_out = strlen(to_out);
	done = 1;
      }
      fprintf(Log,"%d:   %s %s ",I,SPT,LUM);
    }
    if(pec == 0 && flagfe == 0) {
      if(I<=NI) {
	if(To_out <= 9) sprintf(tmpout,"\t\t\t%s\n",qual);
	else if(To_out > 9 && To_out <= 17) sprintf(tmpout,"\t\t%s\n",qual);
        else sprintf(tmpout,"\t%s\n",qual);
	strcat(Iter[I].formatout,tmpout);
      }
      if(I==NI) {
	J = findbest(Iter,NI);
	fprintf(out,"%s\n",Iter[J].formatout);
        fprintf(Log,"\nBest iteration: I = %d\n",J);
	done = 1;
      }
    }
    else {
      if(I<=NI) {
	To_out += strlen(PEC);
	if(To_out <=9) sprintf(tmpout,"%s \t\t\t%s",PEC,qual);
	else if(To_out > 9 && To_out <= 17) sprintf(tmpout,"%s \t\t%s",PEC,qual);
        else sprintf(tmpout,"%s \t%s",PEC,qual);
	strcat(Iter[I].formatout,tmpout);
      }
      if(I==NI) {
	J = findbest(Iter,NI);
	fprintf(out,"%s\n",Iter[J].formatout);
        fprintf(Log,"\nBest iteration: I = %d\n",J);
	done = 1;
      }
      fprintf(Log,"%s %7.1e\n",PEC,CHI2);
    }
  }
}

// Module for K and M type stars
void KMclass()
{
  extern int sf;
  extern float Lumcode;
  extern float Spcode;
  float ax,bx,cx,fa,fb,fc;
  float q,sptkm,lumkm,spth,sptm;
  int pec = 0;
  int he,sr,si,eu,ba,cr,ch,cn,Ba;
  char PEC[20];
  float CHI2 = 0.0;
  float C2,CN,CH;
  float lumBa;
  int i,J;
  int flagmetal = 0;

  he = sr = si = eu = cr = ba = ch = cn = Ba = 0;
  strcpy(PEC,"");

  // Iterate on temperature and luminosity types twice

  fprintf(Log,"\nClassifying this star as a K/M star\n");
  iterate++;

  /* First look at the hydrogen-line spectral type; hydrogen lines 
don't constitute a good criterion for K/M stars, but can exclude a
G spectral type.  Note: hydrogen index will give a latest hydrogen-line
type of K2. */

  ax = spt-1;
  bx = spt;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,hydrogen_index);
  q = brent(ax,bx,cx,hydrogen_index,0.001,&spth);
  code2spt(spth);
  fprintf(Log,"%d: Hydrogen-line type = %s\n",I,SPT);

  /* Now the metallic-line type.  Again, not great for K/M stars,
     but this can help detect an early K, metal-weak star */

  ax = spth-1;
  bx = spth;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptmetal);
  q = brent(ax,bx,cx,sptmetal,0.001,&sptm);
  code2spt(sptm);
  fprintf(Log,"%d: Metallic-line type = %s\n\n",I,SPT);

  if(spth < 33.0) {
    spt = spth;
    FGclass();
    return;
  }

  /* If the following is satisfied, the star is a metal-weak early K-type
star.  Caution, application of this criterion beyond K5 will pick up the
natural decline in blue-violet line strengths in late K and M stars because
of increased violet opacity */
  if(spth >= 34.0 && spth <= 36.0 && sptm < 32.5) flagmetal = 1;
  fprintf(Log,"flagmetal = %d\n",flagmetal);
  // Check to see if there is a Barium peculiarity
  Ba = barium(spt,lum);

  // If star is normal, iterate on temperature and luminosity classifications
  To_out = 0;
  if(flagmetal == 0) {
    sptkm = spt;
    Spcode = sptkm;
    lumkm = lum;
    Lumcode = lumkm;
    for(i=0;i<1;i++) {
      ax = sptkm-1;
      bx = sptkm;
      mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptKM);
      q = brent(ax,bx,cx,sptKM,0.001,&sptkm);
      if(sptkm > scool) {
        sptkm = scool;
	fprintf(Log,"Caution, spectral type may be cooler than library limit\n");
      }
      Spcode = spt = sptkm;
      code2spt(spt);
      fprintf(Log,"%d: sptkm = %f\n",I,spt);
      ax = lumkm;
      bx = lumkm + 0.5;
      mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,lumratiomin);
      q = brent(ax,bx,cx,lumratiomin,0.001,&lumkm);
      if(lumkm > 5.2) lumkm = Lumcode = 5.2;
      /* because of deficiency of standard library (no supergiants after
	 M2), if lum < 3.0 for stars later than M2, set lum = 3.0 */
      if(spt > 42.5 && lumkm < 3.0) lumkm = 3.0;
      /* no dwarf standards later than library limit, so adjust everything to giant */
      if(spt >= scool && lumkm > 3.0) lumkm = 3.0;
      Lumcode = lum = lumkm;
      code2lum(lum);
      // Check on carbon peculiarities
      fprintf(Log,"C2 ratio = %f\n",carbon4737(spt,lum));
      fprintf(Log,"CN ratio = %f\n",CN4215(spt,lum));
      C2 = CN = CH = 0.0;
      if(spt < 39.0 && lum < 4.5) {
        C2 = carbon4737(spt,lum);
        CN = CN4215(spt,lum);
        CH = CHband2(spt,lum);
      }
      fprintf(Log,"CH ratio = %f\n",CH);
      if(CN >= 0.025 && CN < 0.05) strcat(PEC,"");
      if(CN >= 0.05 && CN < 0.075) strcat(PEC," CN1");
      if(CN >= 0.075) strcat(PEC,"CN2");
      if(CN < -0.025 && CN > -0.05) strcat(PEC,"");
      if(CN <= -0.05) strcpy(PEC," CN-1");
      if(CH >= 0.075 && CH < 0.150) strcat(PEC,"");
      if(CH >= 0.150 && CH < 0.200) strcat(PEC," CH1");
      if(CH >= 0.200) strcat(PEC," CH2");
      if(CH <= -0.075 && CH > -0.150) strcat(PEC,"");
      if(CH <= -0.150 && CH > -0.200) strcat(PEC," CH-1");
      if(CH < -0.200) strcat(PEC," CH-2");
      Ba = barium(spt,lum);
      // Suppress Ba peculiarity for M stars.
      if(Ba == 1 && spt < 40.0) strcat(PEC," Ba");
      fprintf(Log,"%d:  lumkm = %f\n",I,lumkm);
    }

    if(sptkm < 33.0 && iterate < 5) {
      spt = sptkm;
      Spcode = spt;
      lum = lumkm;
      Lumcode = lum;
      FGclass();
      return;
    }
    // Slight tendency to classify K and M dwarfs as IV-V
    if(I==NI) if(lum > 4.50) lum = 5.0;

    Spcode = spt;
    Lumcode = lum;
    code2lum(lum);
    code2spt(spt);
    CHI2 = match(spt,lum);
    quality(CHI2,qual,spt);
    sprintf(to_out,"%s %s %s",SPT,LUM,PEC);
    To_out = strlen(to_out);
    if(I<=NI) {
      if(To_out <= 9) sprintf(Iter[I].formatout,"%s %s %s \t\t\t%s",SPT,LUM,PEC,qual);
      else if(To_out > 9 && To_out <= 17) sprintf(Iter[I].formatout,"%s %s %s \t\t%s",SPT,LUM,PEC,qual);
      else sprintf(Iter[I].formatout,"%s %s %s \t%s",SPT,LUM,PEC,qual);
      Iter[I].chi2 = CHI2;
    }
    if(I == NI) {
      J = findbest(Iter,NI);
      fprintf(out,"%s\n",Iter[J].formatout);
      fprintf(Log,"Best iteration: I = %d\n",J);
      done = 1;
    }
    fprintf(Log,"%d:  %s %s %s %7.1e\n",I,SPT,LUM,PEC,CHI2);
  } else {
    // Here we deal with metal-weak early K-type stars
    code2spt(spth);
    ax = lum - 0.5;
    bx = lum;
    mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,lumratiomin);
    q = brent(ax,bx,cx,lumratiomin,0.001,&lumkm);
    if(lumkm > 5.2) lumkm = Lumcode = 5.2;
    code2lum(lumkm);
    Lumcode = lum = lumkm;
    CHI2 = match(spt,lum);
    quality(CHI2,qual,spt);
    sprintf(to_out,"%s %s metal-weak",SPT,LUM);
    To_out = strlen(to_out);
    if(I<=NI) {
      if(To_out <= 9) sprintf(Iter[I].formatout,"%s %s metal-weak \t\t\t%s",SPT,LUM,qual);
      else if(To_out > 9 && To_out <= 17) sprintf(Iter[I].formatout,"%s %s metal-weak \t\t%s",SPT,LUM,qual);
      else sprintf(Iter[I].formatout,"%s %s metal-weak \t%s",SPT,LUM,qual);
      Iter[I].chi2 = CHI2;
    }
    if(I==NI) {
      J = findbest(Iter,NI);
      fprintf(out,"%s\n",Iter[J].formatout);
      fprintf(Log,"Best iteration: I = %d\n",J);
      done = 1;
    }
    fprintf(Log,"%d:  %s %s metal-weak %7.1e\n",I,SPT,LUM,CHI2);
  }

  fprintf(Log,"KM return\n");
  return;
}

// Returns a rough running number equivalent of a temperature type
float sptcode(char spt[])
{
  float code = -10.0;

  if(strncmp(spt,"O",1) == 0) code = 3.0;
  if(strncmp(spt,"B",1) == 0) code = 12.0;
  if(strncmp(spt,"A",1) == 0) code = 20.0;
  if(strncmp(spt,"F",1) == 0) code = 26.0;
  if(strncmp(spt,"G",1) == 0) code = 32.0;
  if(strncmp(spt,"K",1) == 0) code = 36.0;
  if(strncmp(spt,"M",1) == 0) code = 43.0;

  return(code);
} 

// I forget what this one does !
float subclass(float sp)
{
  float spcode[21] = {20.0,21.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,
		    30.0,31.0,32.0,33.0,34.0,35.0,36.0,37.0,38.0,39.0,
		    40.0,41.0};
  float sub[21] = {0.0,2.0,5.0,7.0,8.0,10.0,11.0,13.0,14.0,15.0,17.0,20.0,
		 23.0,25.0,26.0,27.0,28.0,29.0,30.0,32.0,35.5};
  int i,k;
  float subclass;

  for(i=0;i<20;i++) {
    if(sp >= spcode[i] && sp < spcode[i+1]) {
      k = i;
      break;
    }
  }

  subclass = sub[k] + (sub[k+1]-sub[k])*(sp-spcode[k])/(spcode[k+1]-spcode[k]);
  return(subclass);
}

// Rebins the spectrum
void rebin(x,y,k,x2,y2,l1,l2,start,end,dw)
float *x,*y;
float *x2,*y2;
float start,end,dw;
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

// Reads in the ASCII spectrum ignoring headers/tailers denoted with '#'
void getspectrum(infile,x,y,k)
     char infile[];
     float x[],y[];
     int *k;
{

  int l;
  char buffer[200],*tmp;
  FILE *in;

  in = fopen(infile,"r");
  l = 0;
  while(fgets(buffer,190,in) != NULL) {
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

// Pulls out the library match to the program spectrum
float match(float Spt,float Lum)
{
  float *x,*y;
  FILE *mch;
  int i,n,k;
  float chi2;

  x = vector(0,N);
  y = vector(0,N);

  sp2class(Spt,Lum,0.0,x,y,&n);

  mch = fopen(libmatch,"w");

  for(i=0;i<n;i++) fprintf(mch,"%f %f\n",x[i],y[i]);

  // Compute total chi2

  chi2 = 0.0;
  k=0;
  /* If spectral type is later than M0.5, use only wave > 4400.0 in the
     chi2 computation.  Also corrected error in Boolean expression (v1.06) */
  if(Spt < 41.0) {
    for(i=0;i<n;i++) {
      if(x[i] >= wlow+100.0 && x[i] <= whigh-100) { 
        chi2 += (y[i]-Y[i])*(y[i]-Y[i]);
        k++;
      }
    }
  } else {
    for(i=0;i<n;i++) {
      if((x[i] >= wlow+100.0 && x[i] >= 4400.0) && x[i] <= whigh-100.0) {
	chi2 += (y[i]-Y[i])*(y[i]-Y[i]);
        k++;
      }
    }
  } 

  fclose(mch);
  free_vector(x,0,N);
  free_vector(y,0,N);

  return(chi2/(float)k);
}

// Assigns quality tag
void quality(float chi2, char *qual, float spt)
{

  if(spt >= 39.0) {
    if(chi2 < 1.0e-4)      strcpy(qual,"|  excel  |");
    else if(chi2 < 5.0e-3) strcpy(qual,"|  vgood  |");
    else if(chi2 < 5.0e-2) strcpy(qual,"|  good   |");
    else if(chi2 < 1.0e-1) strcpy(qual,"|  fair   |");
    else strcpy(qual,"|  poor   |");
  } else {
    if(chi2 < 1.0e-4)      strcpy(qual,"|  excel  |");
    else if(chi2 < 1.0e-3) strcpy(qual,"|  vgood  |");
    else if(chi2 < 1.0e-2) strcpy(qual,"|  good   |");
    else if(chi2 < 5.0e-2) strcpy(qual,"|  fair   |");
    else strcpy(qual,"|  poor   |");
  }
  strcat(qual,Note);
  return;
}

/* This routine determines the initial rough type using various indices and
   ratios. */
void roughtype2(char *lib, char *name, float *isp, float *ilt)
{
  float *x0,*y0,*x,*flx,xt,yt;
  int i,j,n,k,l1,l2;
  double C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13;
  double CaK,HeI4026,Hdelta,CaI,Gband,HeI4471,TiO;
  double sp1,sp2,sp3,spt,close;
  double t[31] = {9.0,10.0,12.0,13.0,14.0,15.0,16.0,17.0,19.0,20.0,21.0,23.0,
                 24.0,25.0,26.0,27.0,28.0,29.0,30.0,31.0,32.0,33.0,34.0,36.0,
                 37.0,39.0,40.0,40.7,42.5,44.0,45.5};

  x0 = (float *) calloc(10000,sizeof(float));
  y0 = (float *) calloc(10000,sizeof(float));
  x = (float *) calloc(N,sizeof(float));
  flx = (float *) calloc(N,sizeof(float));

  k = 0;
  getspectrum(name,x0,y0,&k);
  rebin(x0,y0,k,x,flx,&l1,&l2,wlow,whigh,space);

  n = (whigh - wlow)/space;

  TiO = 0.0;
  C1=C2=C3=C4=C5=C6=C7=C8=C9=C10=C11=C12=C13 = 0.0;
  for(i=0;i<n;i++) {
	  if(x[i] >= 3918.0 && x[i] <= 3925.0) C1 += flx[i];
	  if(x[i] >= 3927.0 && x[i] <= 3937.0) CaK += flx[i];
	  if(x[i] >= 4022.0 && x[i] <= 4052.0) C2 += flx[i];
	  if(x[i] >= 4062.0 && x[i] <= 4142.0) Hdelta += flx[i];
	  if(x[i] >= 4152.0 && x[i] <= 4182.0) C3 += flx[i];
	  if(x[i] >= 4014.0 && x[i] <= 4020.0) C4 += flx[i];
	  if(x[i] >= 4020.0 && x[i] <= 4032.0) HeI4026 += flx[i];
          if(x[i] >= 4032.0 && x[i] <= 4038.0) C5 += flx[i];
          if(x[i] >= 4210.0 && x[i] <= 4221.0) C6 += flx[i];
          if(x[i] >= 4221.0 && x[i] <= 4232.0) CaI += flx[i];
          if(x[i] >= 4232.0 && x[i] <= 4243.0) C7 += flx[i];
          if(x[i] >= 4233.0 && x[i] <= 4248.0) C8 += flx[i];
          if(x[i] >= 4297.0 && x[i] <= 4314.0) Gband += flx[i];
          if(x[i] >= 4355.0 && x[i] <= 4378.0) C9 += flx[i];
          if(x[i] >= 4452.0 && x[i] <= 4462.0) C10 += flx[i];
          if(x[i] >= 4462.0 && x[i] <= 4480.0) HeI4471 += flx[i];
          if(x[i] >= 4480.0 && x[i] <= 4490.0) C11 += flx[i];
	  if(scool >= 40.0) {
	    if(x[i] >= 4918.0 && x[i] <= 4948.0) C12 += flx[i];
	    if(x[i] >= 4958.0 && x[i] <= 4988.0) C13 += flx[i];
	  }
  }
  C1 /= 7.0;
  CaK /= 10.0;
  C2 /= 30.0;
  Hdelta /= 80.0;
  C3 /= 30.0;
  C4 /= 6.0;
  HeI4026 /= 12.0;
  C5 /= 6.0;
  C6 /= 11.0;
  CaI /= 11.0;
  C7 /= 11.0;
  C8 /= 15.0;
  Gband /= 17.0;
  C9 /= 23.0;
  C10 /= 10.0;
  HeI4471 /= 12.0;
  C11 /= 10.0;
  CaK /= C1;
  Hdelta /= (C2 + C3);
  HeI4026 /= (C4 + C5);
  CaI /= (C6 + C7);
  Gband /= (0.484*C8 + .516*C9);
  HeI4471 /= (C10 + C11);
  if(scool >= 40.0) TiO = C12/C13;

  fprintf(Log,"Gband = %f CaI = %f HeI4471 = %f CaK = %f TiO = %f\n",Gband,CaI,HeI4471,CaK,TiO);


  sp1 = GbandCaI(CaI,Gband);
  sp2 = CaKHe4471(CaK,HeI4471);
  sp3 = scool;
  if(scool >= 40.0) sp3 = TiOIndex(TiO);
  fprintf(Log,"Roughtype2: sp1 = %f  sp2 = %f  sp3 = %f\n",sp1,sp2,sp3);

  spt = 20.0;
  if(sp1 < 0) spt = sp3;
  if(sp1 >= 0 && sp1 <= 20) spt = sp2;
  if(sp2 >= 20 || sp1 >= 20) spt = sp1;
  // if(sp1 <= 20 && sp2 <= 20 && sp3 > 39.0) spt = (max(sp1,sp2)+sp3)/2.0;
  if(sp1 <= 0 && sp2 <=0) spt = sp3;
  if(scool >= 40.0 && sp3 >= 40.0 && sp1 >= 38.0) spt = sp3;
  if(TiO > 1.10) spt = sp3; /* ensures M-type dwarfs are not missed */
  if(spt > scool) spt = scool - 2.0;
  if(spt < shot) spt = shot + 2.0;
  fprintf(Log,"Roughtype2: final spt = %f\n",spt);

  close = 50.0;
  for(i=0;i<31;i++) {
    if(fabs(spt-t[i]) <= close) {
      close = fabs(spt-t[i]);
      j = i;
    }
  }
  spt = t[j];
 
  *isp = spt;
  *ilt = 5.0;  
  return;
}

// Gband and Ca I index used in Roughtype 2
double GbandCaI(double x_in, double y_in)
{
double temp;
temp = 0.0;
// coefficients
double a = -1.3053346170487083E+02;
double b = 3.0348986510742679E+02;
double c = 3.2394730249323806E+02;
double d = -3.1688518509163600E+02;
double f = -1.8121782978289525E+02;
double g = -1.4438268090700953E+02;
temp = a;
temp += b * x_in;
temp += c * y_in;
temp += d * x_in*x_in;
temp += f * y_in*y_in;
temp += g * x_in * y_in;
return temp;
}

// CaK, He 4471 index used in Roughtype 2
double CaKHe4471(double x_in, double y_in)
{
double temp;
temp = 0.0;
// coefficients
double a = -1.6177868547180907E+02;
double b = 2.0657351636101680E+02;
double c = 1.3934317702166081E+02;
double d = -3.8066771401657036E+01;
double f = 1.4860705532285428E+02;
double g = -2.3867986009510940E+02;
temp = a;
temp += b * x_in;
temp += c * y_in;
temp += d * x_in*x_in;
temp += f * y_in*y_in;
temp += g * x_in * y_in;
return temp;
}

// TiO Index used in Roughtype 2
double TiOIndex(double x)
{
  double TiO;

  // TiO = 26.8899 + 16.8833*x -3.18332*x*x;
  TiO = 26.8899 + 16.8833*x;
  if(TiO > 45.5) TiO = 45.5;
  return(TiO);
}

/* This function tries to decide if a star is a lambda boo, but is superceded
   with the lam_boo2 function below */
int lam_boo(float spcode, float lumcode, float *htype, float *metaltype)
{
  float *x,*y;
  extern float *X,*Y;
  extern int N;
  float ax,bx,cx,fa,fb,fc;
  float spth,sptm,q;
  int lb = 0;

  /* Let us assume the star is a dwarf, and then find the spectral type
     that best matches the hydrogen lines */

  Lumcode = 5.0;

  fprintf(Log,"Classifying this star as a Lambda Boo star\n");

  /* Start H-line type out as F0, so that the routine has to work back
to the hydrogen maximum, rather than starting near the hydrogen maximum
and spuriously moving into the B-type stars */

  spth = 23.0;

  ax = spth - 1;
  bx = spth;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,hydrogen_profile);
  q = brent(ax,bx,cx,hydrogen_profile,0.001,&spth);
  fprintf(Log,"Hydrogen-line type = %f\n",spth);
  *htype = spth; 

  /* Now the metallic-line type */

  ax = spt-1;
  bx = spt;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptmetal);
  q = brent(ax,bx,cx,sptmetal,0.001,&sptm);
  fprintf(Log,"Metallic-line type = %f\n",sptm);
  *metaltype = sptm;

  /* If the hydrogen line type is much later than F0, it probably is not a
     Lambda Boo star, but a metal-weak FBS, although this point is debateable */
  if(spth-sptm > 2.0 && spth >= 17.5 && spth <= 23.5) lb = 1;
  if(spth < 17.5) lb = 1;  /* Note: if the hydrogen-line type is too early, 
                              then determining the metallic-line type is 
                              too difficult.  So, we are relying here for
                              the classification solely on the weakness 
                              of the Mg II line.  We might add in a criterion
                              based on the K-line?? */

  return(lb);
}

/* This function tries to decide whether an A-type star is indeed metal
weak, and then attempts to distinguish between a Lambda Boo star and another
type of metal-weak A-type star, such as a horizontal-branch star.  In
order to not over interpret the spectrum, metal-weak non-Lam Boo stars
are simply labelled "metal-weak".  This code also takes care not to label
luminous A-type stars, which can have weak 4481 lines as Lambda Boo. */

int lam_boo2(float spcode, float lumcode, float *htype, float *metaltype)
{
  float *x,*y;
  extern float *X,*Y;
  extern int N;
  float ax,bx,cx,fa,fb,fc;
  float spth,sptm,sptK,q,lum,ratio;
  int lb = 0;

  /* The following is to exclude stars which are luminous from being
     incorrectly classified as Lambda Boo stars */
  if(lumcode < 3.0) return(0);

  /* We assume a preliminary luminosity type of IV-V, which is a compromise
between the dwarf status of Lambda Boo stars and IV or III type of HB stars.
We use that preliminary luminosity type to find the spectral type that best 
matches the hydrogen lines, and then we iterate further. */

  Lumcode = 4.5;

  fprintf(Log,"Classifying this star as a metal-weak A-type star. spt = %f\n",spt);

  /* Start H-line type out as F0, so that the routine has to work back
to the hydrogen maximum, rather than starting near the hydrogen maximum
and spuriously moving into the B-type stars */

  spth = 23.0;

  ax = spth - 1;
  bx = spth;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,hydrogen_profile);
  q = brent(ax,bx,cx,hydrogen_profile,0.001,&spth);
  fprintf(Log,"Hydrogen-line type = %f\n",spth);
  *htype = Spcode = spth; 

  /* If the hydrogen-line type is earlier than F0, we iterate on the luminosity
     type and the hydrogen-line type */
  if(spth < 23.0) {
    lum = Lumcode;
    ax = lum-0.5;
    bx = lum;
    mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,lumratiomin);
    q = brent(ax,bx,cx,lumratiomin,0.001,&lum);
    if(lum > 5.2) lum = 5.2;
    Lumcode = lum;
    fprintf(Log,"Luminosity type = %f\n",lum);

    ax = spth - 1;
    bx = spth;
    mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,hydrogen_profile);
    q = brent(ax,bx,cx,hydrogen_profile,0.001,&spth);
    fprintf(Log,"Hydrogen-line type = %f\n",spth);
    *htype = Spcode = spth;
  }

  /* Now the metallic-line type */

  // first, the K-line type.
  ax = 16.0;
  bx = 17.0;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptCaK);
  q = brent(ax,bx,cx,sptCaK,0.001,&sptK);
  fprintf(Log,"K-line type = %f\n",sptK);

  // then the metallic-line type

  ax = 16.0;
  bx = 17.0;
  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,sptmetal);
  q = brent(ax,bx,cx,sptmetal,0.001,&sptm);
  fprintf(Log,"Metallic-line type = %f\n",sptm);

  /* If the dispersion is low, the metallic-line type can be faulty, 
     especially in metal-weak stars.  So, if sptm < sptK, then correct 
     sptm to sptK */

  if(sptm < sptK - 0.5) {
    sptm = sptK;
    fprintf(Log,"Correcting metallic-line type to K-line type\n");
  }

  *metaltype = sptm;

  /* If the hydrogen line type is much later than F0, it probably is not a
     Lambda Boo star, but a metal-weak FBS, although this point is debateable.
     We limit LB classifications to stars with H - M spectral type differences
     greater than 2.5, to keep high vsini stars from being spuriously classified
     as LBs */
  if(spth-sptm >= 2.5 && spth >= 17.5 && spth <= 23.5) lb = 1;
  if(spth-sptm < 2.5) lb = 3; //More likely these stars are high vsini instead of mild LB.
  if(spth < 17.5) lb = 1;  /* Note: if the hydrogen-line type is too early, 
                              then determining the metallic-line type is 
                              too difficult.  So, we are relying here for
                              the classification solely on the weakness 
                              of the Mg II line.  We might add in a criterion
                              based on the K-line?? */

  /* We use the luminosity type to differentiate between Lambda Boo and
     FHB stars.  For stars earlier than A5, we use the actual luminosity
type, determined above.  For stars later than A5, however, we use the ratio
between Ca I 4226 and Fe II 4233. */

  if(spth < 20.0) {
    if(lum < 4.3 && lb == 1) lb = 2;
  } else {
    ratio = ratioCaIFeII();
    fprintf(Log,"CaIFeII ratio = %f\n",ratio);
    if(ratio < 1.3) lb = 2;
  }
  /* Check once again the metal - hydrogen difference */
  if(fabs(spth-sptm) < 2.5) lb = 3;

  /* If the star is a Lambda Boo star, assume the luminosity class is V
     and redetermine the hydrogen-line type */

  return(lb);
}

double max(double x, double y)
{
  if(x >= y) return(x);
  else return(y);
}

int findbest(results *Kinter,int ni)
{
  int i,j;
  double chimin = 1.0e+30;

  j = ni;
  for(i=1;i<=ni;i++) {
    if(Kinter[i].chi2 <= chimin) {
      chimin = Kinter[i].chi2;
      j = i;
    }
  }
  return(j);
}
