/**
 * mesh_query.c - Query the properties of the mesh at a point
 *
 * Copyright (c) 2010 Leonardo Ramirez-Guzman
 * All rights reserved.  May not be used, modified, or copied 
 * without permission.
 *
 * Leonardo Ramirez-Guzman
 * leoramirezg@gmail.com or lramirezguzman@usgs.gov
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include <string.h>
#include "Stdio.h"

#include "query_mesh_property.h"

void linspace(double xa,double xb,int numpoints,double *result){

    int i;
    double dx;
  
    dx=(xb-xa)/(numpoints-1);
    
    for (i=0; i<numpoints; i++)	
	result[i]=xa+i*dx;     
    return;    
}


int main(int argc, char **argv)
{
    etree_t *mep;
    FILE *result_fp;
    FILE *inputFile,*fpInputQuery;
    out_hdr_t result_hdr;
    double Lon,Lat,depth, vs, vp,rho;
    static char  inputPath[256], outputPath[256], etreePath[256],auxString[256],etreeInfoPath[256];
    FILE *fpOutProp[3],*fpOutQpQs;
    int typeofquery,numberofstations,iStation,hov;
    double originLon,originLat,originDepth;
    int numLon,numLat, numD,iLon,iLat,iD,i,numAlong,j;
    double stepLon,stepLat,stepD,lengthLon,lengthLat,totDepth,finalLon,finalLat,finalDepth;

    if (argc != 6) {
    	fputs ( "Usage: geodata typeofquery db in rout \n "
		" typeofquery:   (0) points or (1) lines and planes\n"
		"       etree:    path etree\n"
		"     infosim:    path etree data simulation\n"
    		"          in:   input file, query specifications\n"
		"        rout:   output path (directory), output file named output.out\n",stderr);	
    }
    
    typeofquery=atoi(argv[1]);
    strcpy(etreePath,argv[2]);
    strcpy(etreeInfoPath,argv[3]);
    strcpy(inputPath,argv[4]);
    strcpy(outputPath,argv[5]);
    
    fpInputQuery = fopen(inputPath,"r");
    

    /* Open mesh databases */
    if ((mep = etree_open(etreePath, O_RDONLY, 0, 0, 0)) == NULL) {
        fprintf(stderr, "Cannot open mesh element etree %s\n", argv[1]);
        exit(1);
    }

    /* Open 4D output */
    if ((result_fp = fopen(etreeInfoPath, "r")) == NULL) {
        fprintf(stderr, "Cannot open 4D output file %s\n", argv[2]);
        exit(1);
    }

    /* Get the metadata */
    if (fread(&result_hdr, sizeof(out_hdr_t), 1, result_fp) != 1) {
        fprintf(stderr, "Cannot read 4D-out result header info from %s.\n", 
                argv[2]);
        exit(1);
    }


    if(typeofquery == 0){
	 
	sprintf(auxString, "%s/output.out", argv[5]);
	fpOutProp[0] = fopen(auxString,"w");
	fscanf( fpInputQuery,"%d",&numberofstations);
	for (iStation=0; iStation<numberofstations;iStation++){
	    fscanf(fpInputQuery,"%lf %lf %lf", &Lon,&Lat,&depth);

	    if (mesh_point(Lat, Lon, depth, mep,result_hdr, &vp,&vs,&rho) != 0) {
	
		vp=0;
		vs=0;
		rho=0;
	    }
	    
	    //fprintf(fpOutProp[0],"\n%f %f %f %f %f %f", Lat,Lon,depth,vp,vs,rho);
	    	    fprintf(fpOutProp[0],"\n%f ",vs);

	}	    		
	
	fclose(fpOutProp[0]);
	
    }else{	

	/* AUXILIARY FILES TO OUTPUT PROPERTIES */
	sprintf(auxString, "%s/vp.out", argv[5]);
	fpOutProp[0] = fopen(auxString,"w");
	sprintf(auxString, "%s/vs.out", argv[5]);
	fpOutProp[1] = fopen(auxString,"w");
	sprintf(auxString, "%s/rho.out", argv[5]);
	fpOutProp[2] = fopen(auxString,"w");
	sprintf(auxString, "%s/qpqs.out", argv[5]);
	fpOutQpQs  = fopen(auxString,"w");
		
	/* DEFINITION OF THE BOX WHERE WE ARE QUERYING THE DATA */
	/* IT CAN BE A HORIZONTAL PLANE ( 0 ) OR A VERTICAL PLANE (1)*/	
	Parse_Text( fpInputQuery,"hororvert"  , 'i',&hov);
	Parse_Text( fpInputQuery,"originLon"  , 'd',&originLon);	 
	Parse_Text( fpInputQuery,"originLat"  , 'd',&originLat);	 
	Parse_Text( fpInputQuery,"originDepth", 'd',&originDepth);
	
	if(hov == 0 ){
	    Parse_Text( fpInputQuery,"lengthLon", 'd',&lengthLon);
	    Parse_Text( fpInputQuery,"lengthLat", 'd',&lengthLat);
	    Parse_Text( fpInputQuery,"numLon", 'i',&numLon);
	    Parse_Text( fpInputQuery,"numLat", 'i',&numLat);
	    if(numLon > 1) stepLon = lengthLon /(numLon-1);
	    else           stepLon = lengthLon;
	    if(numLat > 1) stepLat = lengthLat /(numLat-1);
	    else           stepLat = lengthLat;
	    numD = 1;
	    stepD   = 1;
	    
	    /* OUTPUT PLANES */
	    /* WRITE HEADER  */
	    for (i=0;i<3;i++)
		fprintf(fpOutProp[i],"%d %e %e %e %d %d %d %e %e %e",hov,
			originLon, originLat,originDepth,
			numLon,    numLat,       numD,
			lengthLon, lengthLat,   totDepth);	    
	    for (iLon=0;iLon < numLon; iLon++){	
		Lon=originLon+stepLon*iLon;		
		for (iLat=0;iLat < numLat; iLat++){	    
		    Lat= originLat+stepLat*iLat;		    
		    for (iD=0; iD < numD ; iD++){	      
			depth=originDepth+stepD*iD;
			if (mesh_point(Lat, Lon, depth, mep,result_hdr, &vp,&vs,&rho) != 0) {			    
			    vp=0;
			    vs=0;
			    rho=0;
			}					
			fprintf(fpOutProp[0]," %e ", vp );
			fprintf(fpOutProp[1]," %e ", vs );
			fprintf(fpOutProp[2]," %e ", rho);		
		    }	
		}
	    }    	
	}
	
	if(hov == 1 ){
	    Parse_Text( fpInputQuery,"finalLon"  , 'd',&finalLon);
	    Parse_Text( fpInputQuery,"finalLat"  , 'd',&finalLat);
	    Parse_Text( fpInputQuery,"finalDepth", 'd',&finalDepth);
	    Parse_Text( fpInputQuery,"numAlong", 'i',&numAlong);
	    Parse_Text( fpInputQuery,"numD"    , 'i',&numD);	
	    double *longitudes, *latitudes, *depths;
	    longitudes=malloc(sizeof(double)*numAlong);
	    latitudes =malloc(sizeof(double)*numAlong);
	    depths    =malloc(sizeof(double)*numD);
	    linspace(originLon  ,finalLon,numAlong,longitudes);
	    linspace(originLat  ,finalLat,numAlong,latitudes );
	    linspace(originDepth,finalDepth,  numD,depths    );
	    /* PRINT HEADER */
	    for (i=0;i<3;i++)
		fprintf(fpOutProp[i],"%d %e %e %e %e %e %e %d %d \n",hov,
			originLon, originLat,originDepth,finalLon,finalLat,finalDepth,
			numAlong,       numD);
	    for (i=0;i<numAlong;i++)
		for (j=0;j<numD;j++){		
		    
		    if (mesh_point( latitudes[i],longitudes[i],depths[j], mep,result_hdr, &vp,&vs,&rho) != 0) {			    
			vp=-1;
			vs=-1;
			rho=-1;
		    }
		    
		    fprintf(fpOutProp[0]," %e ", vp );
		    fprintf(fpOutProp[1]," %e ", vs );
		    fprintf(fpOutProp[2]," %e ", rho);	
		}
	}
    }
    

    etree_close(mep);
   
    return 0;
}



