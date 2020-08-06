/******************************************************************************
 * Code: geodataq_version2
 *
 * Spatial Velocity Model Database Manager. 
 *
 * Copyright (c) 2006 (c) 2011 Leonardo Ramirez-Guzman
 * ALL RIGHTS RESERVED. MAY NOT BE USED, MODIFIED OR COPIED WITHOUT 
 * PERMISSION.  
 *
 * Contact:
 * Leonardo Ramirez-Guzman
 * lramirezguzman@usgs.gov or leoramirezg@gmail.com
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "geodataq.h"

/*
 * print_usage()
 */
void  print_usage(){
  fputs ( "Usage: geodataq typeofquery velorveldepthorbdrortopo db in rout \n "
	  "              typeofquery:  (0) points or (1) lines or planes\n"
	  " velorveldepthorbdrortopo:   (0) velocity & density (asl) \n"
	  "                             (1) velocity & density (depth from surface) \n"
	  "                             (2) velocity & density (squeeze to msl != 1) \n"
	  "                             (3) elevation topo\n"
	  "                             (4) elevation bedrock\n"
	  "                             (5) depth to bedrock\n"
	  "                       db:   path to a database\n"
	  "                       in:   input file, query specifications\n"
	  "                     rout:   output path (directory), \n"
	  "                             output file named output.out\n",stderr);
  exit(1);
}


/******************************************************************* ********* 
 *                             Main                                          
 *                                                                           
 *   Notes:                                                                  
 *         Boundaries of units have to be given in geodetic coordinates,     
 *         i.e. lon,lat, elevation (msl), they will be shifted to datum      
 *         to use depth from theDatum value.                                 
 *                                                                           
 *****************************************************************************/
int main (int argc, char **argv){

  int typeofquery, velorveldepthorbdrortopo,squeeze;    
  FILE *fpInputQuery;
  static char inputPlanesQueryPath[256], outputPath[256],databasePath[256];
  database_t vm; 

  /* check input args */
  if (argc != 6) print_usage();
   
  typeofquery             =atoi(argv[1]);
  velorveldepthorbdrortopo=atoi(argv[2]);
   
  strcpy(databasePath     ,argv[3]);
  strcpy(inputPlanesQueryPath,argv[4]);
  strcpy(outputPath          ,argv[5]);
    
  fpInputQuery = fopen(inputPlanesQueryPath,"r");
    
  Init_Database(&vm,databasePath,1);    

  if(typeofquery == 0){ /* point by point */
    switch(velorveldepthorbdrortopo) 
      { /*  properties can include velocity,density or containing unit */
      case 0: extract_print_point_query(argv,&vm,0);break; 
      case 1: extract_print_point_query(argv,&vm,1);break; 
      case 2: extract_print_point_query(argv,&vm,2);break; 
      }
  }
  else{
    switch(velorveldepthorbdrortopo)
     {	/*  properties can include velocity,density or containing unit */
     case 0:extract_print_properties_grid (argv,&vm,0);break;
     case 1:extract_print_properties_grid (argv,&vm,1);break;
     case 2:extract_print_properties_grid (argv,&vm,2);break;
     case 3:extract_print_bedrockelevation(argv,&vm)  ;break;
     case 4:extract_print_topoelevation   (argv,&vm)  ;break;
     case 5:extract_print_bedrockdepth    (argv,&vm)  ;break;
     }
  }    
  return 0;
} 

