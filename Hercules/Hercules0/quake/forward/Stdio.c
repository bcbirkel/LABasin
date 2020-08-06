/*
 *
 *  Stdio.c: some stdio.h functions with automatic error and exit
 *
 *  Copyright (c) 2006 Leonardo Ramirez-Guzman
 *  All rights reserved.  May not be used, modified, or copied
 *  without permission.
 *
 *  Contact:
 *  Leonardo Ramirez-Guzman
 *  Civil and Environmental Engineering
 *  Carnegie Mellon University
 *  5000 Forbes Avenue
 *  Pittsburgh, PA 15213
 *  lramirez@andrew.cmu.edu
 *
 */

#include "Stdio.h"

/**
 *   Fopen: same as fopen
 *
 *  return: file pointer of the opened file
 *
 */
FILE* Fopen( char filetoopen[256], char type[5] ){
    
    FILE *filepointer;

    filepointer = fopen(filetoopen,type);
    if ( filepointer == NULL) {
	fprintf(stderr, "Error opening %s", filetoopen);
	exit(1);
    }
    
    return filepointer;
}


/*
 *
 *  Fscanf: reduced version of fscanf, it only retrieves one value
 *          at a time.
 *
 *
 */

int Fscanf( FILE *stream, const char *format, void *result ){

    int ret;
        
    if (strcmp("%d", format) == 0)
	ret = fscanf(stream, "%d", (int *)result);
    else if (strcmp("%f", format) == 0)
	ret = fscanf(stream, "%f", (float *)result);
    else if (strcmp("%lf", format) == 0)
	ret = fscanf(stream, "%lf", (double *)result);
    else{
	fprintf(stderr, "Fscanf: %s  not available format", format);
	exit(1);
    }
	          

    if( result == NULL || ret ==0 ){
	fprintf(stderr, "Fscanf: unable to read value");
	exit(1);
    }

    return ret;

}



int FscanfE( FILE *stream, const char *format, void *result,int a,int b ){

    int ret;
        
    if (strcmp("%d", format) == 0)
	ret = fscanf(stream, "%d", (int *)result);
    else if (strcmp("%f", format) == 0)
	ret = fscanf(stream, "%f", (float *)result);
    else if (strcmp("%lf", format) == 0)
	ret = fscanf(stream, "%lf", (double *)result);
    else{
	fprintf(stderr, "Fscanf: %s  not available format", format);
	exit(1);
    }
	          

    if( result == NULL || ret ==0 ){
	fprintf(stderr, "Fscanf: unable to read value, integers %d %d",a,b);
	exit(1);
    }

    return ret;

}


/**
 * parsetext: Parse a text file and return the value of a match string.
 *
 *
 *
 */
void Parse_Text(FILE *fp, const char *querystring, const char type,
              void* result){
    int32_t res = 0, found = 0;

    /* Start from the beginning */
    rewind(fp);

    /* Look for the string until found */
    while (!found) {
        char line[LINESIZE];
        char delimiters[] = " =\n\r\t";
        char *name, *value;

        /* Read in one line */
        if (fgets(line, LINESIZE, fp) == NULL)
            break;
        
        name = strtok(line, delimiters);
	if ((name != NULL) && (strcmp(name, querystring) == 0)) {
            found = 1;
            value = strtok(NULL, delimiters);
            
            switch (type) {
            case 'i':
                res = sscanf(value, "%d", (int *)result);
                break;
	    case 'I':
                res = sscanf(value, "%d", (int32_t *)result);
                break;
            case 'f':
                res = sscanf(value, "%f", (float *)result);
                break;
            case 'd':
                res = sscanf(value, "%lf", (double *)result);
                break;
            case 's':
                res = 1;
                strcpy((char *)result, value);
                break;
            default:
                fprintf(stderr, "parsetext: unknown type %c\n", type);
                exit(1);
            }
        }
        
    }

    if (found == 0){
	fprintf(stderr,"\n Cannot find %s in input file", querystring);
	exit(1);
    }
    if (res == -1){
	fprintf(stderr, "parsetext: unknown type %c\n", type);
	exit(1);
    }
    
    return ;
}


/**
 * parsedarray: Parse array of doubles .
 *
 * - return 0 if OK, -1 on error
 *
 */
int Parse_Double_Array( FILE *fp, const char *querystring,int sizecol, int sizerow,
			double **array){

  rewind (fp);
  int found=0;
  int iSize, jSize;
  
  while (!found) {
    
    char line[LINESIZE];
    char *name;
    char delimiters[] = " =\n";

    /* Read in one line */
    if (fgets(line, LINESIZE, fp) == NULL)
       break;
    
    name = strtok(line, delimiters);
    
    if ( (name != NULL) && (strcmp(name, querystring) == 0)){
       found = 1;
       for ( iSize = 0; iSize < sizerow; iSize++ )
	   for ( jSize = 0; jSize < sizecol; jSize++)
	       Fscanf ( fp,"%lf", &(array[iSize][jSize] ));
	
	return 0;
    }
  }
  return -1;
}


