/*
 *
 *  Stdio.h:
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>
#include <errno.h>
#include <stdarg.h>

#define LINESIZE        512

FILE* Fopen( char filetoopen[256], char type[5] );
int Fscanf( FILE *stream, const char *format, void *result );
int FscanfE( FILE *stream, const char *format, void *result,int a,int b );
void Parse_Text(FILE *fp, const char *querystring, const char type,
		       void* result);
int Parse_Double_Array( FILE *fp, const char *querystring,int sizecol, int sizerow,double **array);

