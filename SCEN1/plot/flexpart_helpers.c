#include <stdio.h>
#include <stdlib.h>
#include "flexpart_helpers.h" /* dimensions of flexpart output */

void read_write(FILE *in, FILE *out, int nlines ){
int i;
float data, err;
    for ( i=0; i <= (nlines-1); ++i) 
       {
        fscanf(in, "%f%f", &data, &err);
        /* printf("%f\t%f\n", data, err); */
        fwrite(&data, sizeof(float), 1, out);
       }
}

void skip_lines(FILE *fpointer, int nlines ) {
  int i;
  char str[100];
    for (i=0; i <= nlines-1; i++)
    {
     fgets(str,1800,fpointer);
    }
}
void skip_lines2(FILE *fpointer, int nlines ) {
  int i;
  char *line = NULL;
  size_t len = 0;
    for (i=0; i <= nlines-1; i++)
    {
     getline(&line, &len, fpointer);
    }
}
