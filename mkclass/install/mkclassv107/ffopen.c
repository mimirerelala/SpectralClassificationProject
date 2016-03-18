// A fancy substitute for fopen
#include <stdio.h>
#include <stdlib.h>

FILE *ffopen(Name,mode)
char *Name,*mode;
{
  FILE *fp;
  extern FILE *Log;
  extern FILE *out;
  extern char name[80];
  extern int flagname;

  if((fp = fopen(Name,mode)) == NULL) {
    // printf("Cannot open file %s\n",Name);
    // printf("Now exiting program\n");
    if(flagname == 0) fprintf(out,"%s  | ?     \t\t\t|         | \\\\\n",name);
    else fprintf(out,"?     \t\t\t|         | \\\\\n");
    fprintf(Log,"Cannot open file %s\n",Name);
    exit(1);
  }
  return(fp);
}
