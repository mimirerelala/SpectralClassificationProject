// Translates running number into spectral and luminosity types
#include <stdio.h>
#include <string.h>
int between(float x,float low, float high);

/* O3 0.0, O4 1.0, O5 2.0, O6 3.0, O7 4.0, O8 5.0, O9 6.0, B0 7.0,
   B1 8.0 */

void code2spt(float spcode)
{
  extern char SPT[10];
  strcpy(SPT,"??");
  if(between(spcode,-0.5,0.5) == 1) strcpy(SPT,"O3");
  if(between(spcode,0.5,1.5) == 1) strcpy(SPT,"O4");
  if(between(spcode,1.5,2.5) == 1) strcpy(SPT,"O5");
  if(between(spcode,2.5,3.5) == 1) strcpy(SPT,"O6");
  if(between(spcode,3.5,4.5) == 1) strcpy(SPT,"O7");
  if(between(spcode,4.5,5.5) == 1) strcpy(SPT,"O8");
  if(between(spcode,5.5,6.5) == 1) strcpy(SPT,"O9");
  if(between(spcode,6.5,7.5) == 1) strcpy(SPT,"B0");
  if(between(spcode,7.5,8.5) == 1) strcpy(SPT,"B1");
  if(between(spcode,8.5,9.5) == 1) strcpy(SPT,"B2");
  if(between(spcode,9.5,10.5) == 1) strcpy(SPT,"B3");
  if(between(spcode,10.5,11.5) == 1) strcpy(SPT,"B4");
  if(between(spcode,11.5,12.3) == 1) strcpy(SPT,"B5");
  if(between(spcode,12.3,12.7) == 1) strcpy(SPT,"B6");
  if(between(spcode,12.7,13.5) == 1) strcpy(SPT,"B7");
  if(between(spcode,13.5,14.5) == 1) strcpy(SPT,"B8");
  if(between(spcode,14.5,15.3) == 1) strcpy(SPT,"B9");
  if(between(spcode,15.3,15.7) == 1) strcpy(SPT,"B9.5");
  if(between(spcode,15.7,16.5) == 1) strcpy(SPT,"A0");
  if(between(spcode,16.5,17.5) == 1) strcpy(SPT,"A1");
  if(between(spcode,17.5,18.5) == 1) strcpy(SPT,"A2");
  if(between(spcode,18.5,19.3) == 1) strcpy(SPT,"A3");
  if(between(spcode,19.3,19.7) == 1) strcpy(SPT,"A4");
  if(between(spcode,19.7,20.3) == 1) strcpy(SPT,"A5");
  if(between(spcode,20.3,20.7) == 1) strcpy(SPT,"A6");
  if(between(spcode,20.7,21.5) == 1) strcpy(SPT,"A7");
  if(between(spcode,21.5,22.1) == 1) strcpy(SPT,"A8");
  if(between(spcode,22.1,22.75) == 1) strcpy(SPT,"A9");
  if(between(spcode,22.75,23.25) == 1) strcpy(SPT,"F0");
  if(between(spcode,23.25,23.75) == 1) strcpy(SPT,"F1");
  if(between(spcode,23.75,24.5) == 1) strcpy(SPT,"F2");
  if(between(spcode,24.5,25.25) == 1) strcpy(SPT,"F3");
  if(between(spcode,25.25,25.75) == 1) strcpy(SPT,"F4");
  if(between(spcode,25.75,26.5) == 1) strcpy(SPT,"F5");
  if(between(spcode,26.5,27.25) == 1) strcpy(SPT,"F6");
  if(between(spcode,27.25,27.75) == 1) strcpy(SPT,"F7");
  if(between(spcode,27.75,28.6) == 1) strcpy(SPT,"F8");
  if(between(spcode,28.6,29.5) == 1) strcpy(SPT,"F9");
  if(between(spcode,29.5,30.25) == 1) strcpy(SPT,"G0");
  if(between(spcode,30.25,30.75) == 1) strcpy(SPT,"G1");
  if(between(spcode,30.75,31.25) == 1) strcpy(SPT,"G2");
  if(between(spcode,31.25,31.50) == 1) strcpy(SPT,"G3");
  if(between(spcode,31.50,31.75) == 1) strcpy(SPT,"G4");
  if(between(spcode,31.75,32.25) == 1) strcpy(SPT,"G5");
  if(between(spcode,32.25,32.55) == 1) strcpy(SPT,"G6");
  if(between(spcode,32.55,32.85) == 1) strcpy(SPT,"G7");
  if(between(spcode,32.85,33.25) == 1) strcpy(SPT,"G8");
  if(between(spcode,33.25,33.75) == 1) strcpy(SPT,"G9");
  if(between(spcode,33.75,34.5) == 1) strcpy(SPT,"K0");
  if(between(spcode,34.5,35.5) == 1) strcpy(SPT,"K1");
  if(between(spcode,35.5,36.5) == 1) strcpy(SPT,"K2");
  if(between(spcode,36.5,37.5) == 1) strcpy(SPT,"K3");
  if(between(spcode,37.5,38.5) == 1) strcpy(SPT,"K4");
  if(between(spcode,38.5,39.4) == 1) strcpy(SPT,"K5");
  if(between(spcode,39.4,39.6) == 1) strcpy(SPT,"K6");
  if(between(spcode,39.6,40.5) == 1) strcpy(SPT,"K7");
  if(between(spcode,40.5,40.8) == 1) strcpy(SPT,"M0");
  if(between(spcode,40.8,41.2) == 1) strcpy(SPT,"M0.5");
  if(between(spcode,41.2,41.8) == 1) strcpy(SPT,"M1");
  if(between(spcode,41.8,42.2) == 1) strcpy(SPT,"M1.5");
  if(between(spcode,42.2,42.8) == 1) strcpy(SPT,"M2");
  if(between(spcode,42.8,43.2) == 1) strcpy(SPT,"M2.5");
  if(between(spcode,43.2,43.8) == 1) strcpy(SPT,"M3");
  if(between(spcode,43.8,44.2) == 1) strcpy(SPT,"M3.5");
  if(between(spcode,44.2,44.8) == 1) strcpy(SPT,"M4");
  if(between(spcode,44.8,45.25) == 1) strcpy(SPT,"M4.5");
  if(between(spcode,45.25,46.0) == 1) strcpy(SPT,"M5");
  if(between(spcode,46.0,47.0) == 1) strcpy(SPT,"M6");
  if(between(spcode,47.0,48.0) == 1) strcpy(SPT,"M7");
  if(between(spcode,48.0,49.0) == 1) strcpy(SPT,"M8");
  if(between(spcode,49.0,50.0) == 1) strcpy(SPT,"M9");

}

int between(float x,float low, float high)
{
  if(x >= low && x < high) return(1);
  else return(0);
}

void code2lum(float lumcode)
{
  extern char LUM[10];
  if(lumcode > 6.0) strcpy(LUM,"?");
  if(lumcode < -1.5) strcpy(LUM,"?");
  if(between(lumcode,5.25,6.0) == 1) strcpy(LUM,"V-");
  if(between(lumcode,4.75,5.25) == 1) strcpy(LUM,"V");
  if(between(lumcode,4.25,4.75) == 1) strcpy(LUM,"IV-V");
  if(between(lumcode,3.75,4.25) == 1) strcpy(LUM,"IV");
  if(between(lumcode,3.25,3.75) == 1) strcpy(LUM,"III-IV");
  if(between(lumcode,2.75,3.25) == 1) strcpy(LUM,"III");
  if(between(lumcode,2.25,2.75) == 1) strcpy(LUM,"II-III");
  if(between(lumcode,1.75,2.25) == 1) strcpy(LUM,"II");
  if(between(lumcode,1.25,1.75) == 1) strcpy(LUM,"Ib-II");
  if(between(lumcode,0.75,1.25) == 1) strcpy(LUM,"Ib");
  if(between(lumcode,0.25,0.75) == 1) strcpy(LUM,"Iab");
  if(between(lumcode,-0.5,0.25) == 1) strcpy(LUM,"Ia");
  if(between(lumcode,-1.5,-0.5) == 1) strcpy(LUM,"0");
  return;
}
