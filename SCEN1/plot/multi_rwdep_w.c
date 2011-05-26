/******************************************************************
 * This program reads multiple output files created by FLEXPART   *
 * and creates a single output including all time steps           *
 * outputs c-binary                                               *
 ******************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "flexpart_helpers.h" /* dimensions of flexpart output */

#define STRSIZ 80
#define NSKIP 2
#define NREAD (NX*NY)


int
  main(void)
  {
    char out_name[STRSIZ]= "./wetdep.dat" ;
    char in_list[STRSIZ] = "./input_list";
    char fnam[STRSIZ];
    FILE *inlstp, *finp,
         *foutp;
    int num, i, status;

    inlstp = fopen(in_list,"r");
    foutp = fopen(out_name,"wb");
    for (status=fscanf(inlstp, "%s", fnam);
    status != EOF;
    status = fscanf(inlstp, "%s", fnam)) 
    {
      printf("Reading  %s\n", fnam); 
      finp  = fopen(fnam,"r");
      skip_lines2(finp,NSKIP); 
      read_write(finp, foutp, NREAD);
      fclose(finp);
    } 
    fclose(inlstp);
  }
