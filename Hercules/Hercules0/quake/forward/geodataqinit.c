#include <stdio.h>
#include <string.h>

#include "Stdio.h"
#include "geodataq.h"

int theIsBinary=0;

/***************************************************************************** 
 *                                                                           *
 *                                                                           *
 *  Init_Surface: INIT A SURFACE DATA STRUCTURE, SIMILAR TO A RASTER FILE    *
 *                                                                           *
 *                                                                           *
 /****************************************************************************/
double  Init_Surface(FILE *fp, surface_t *currentsurface){
    
    int i,iX,iY,iCount, nx, ny; 
    double *auxiliar, maxVal=-1e10;
    char surfacefile[1024];
    
    fscanf(fp," %d %d ", &nx, &ny);
    currentsurface->nxgrid = nx;
    currentsurface->nygrid = ny;    
    currentsurface-> xgrid = (double* ) malloc (sizeof(double)* nx);
    currentsurface-> ygrid = (double* ) malloc (sizeof(double)* ny);
    
    if(  currentsurface-> xgrid == NULL ||  currentsurface-> ygrid == NULL) {
	fprintf(stdout,"Error allocating memory for surfaces");
	fflush(stdout);
	exit(1);
    }
    
    for ( i =0; i< nx; i++)FscanfE(fp,"%lf",&(currentsurface->xgrid[i]),i,0);
    for ( i =0; i< ny; i++)FscanfE(fp,"%lf",&(currentsurface->ygrid[i]),i,0);
    
    /* ALLOCATE AND LOAD THE ZGRID */
    currentsurface->zgrid = (double ** ) malloc( sizeof(double *) * nx);  
    for ( iX = 0; iX < nx; iX++){
	currentsurface->zgrid[ iX ] =  (double *) malloc( sizeof(double) * ny);
	if(currentsurface->zgrid[ iX ] == NULL ) {
	    fprintf(stdout,"Error allocating memory for surfaces");
	    fflush(stdout);
	    exit(1);	    
	}
    }
    
    iCount=0;
    for ( iX = 0; iX < nx; iX++){
	for( iY = 0; iY < ny; iY++){	  
	    FscanfE(fp,"%lf",&( currentsurface->zgrid[iX][iY]),iX,iY);
	    iCount+=1; 
	    if(maxVal<currentsurface->zgrid[iX][iY])
	       maxVal=currentsurface->zgrid[iX][iY];
	}
    }
    
    return maxVal;

}

double  ShiftDatum_Surface( surface_t *currentsurface,double datum){
    
    int i,iX,iY,iCount, nx, ny; 
    double maxVal;

    nx=currentsurface->nxgrid;
    ny=currentsurface->nygrid;    
    iCount=0;
    for ( iX = 0; iX < nx; iX++){
	for( iY = 0; iY < ny; iY++){	  
	    currentsurface->zgrid[iX][iY]=datum-currentsurface->zgrid[iX][iY];
	    iCount+=1; 
	    if(maxVal<currentsurface->zgrid[iX][iY])
	       maxVal=currentsurface->zgrid[iX][iY];
	}
    }
    
    return maxVal;

}



/***************************************************************************** 
 *                                                                           *
 *     Init_Vel_Profiles:  Init velocity distribution datastructure          *
 *                                                                           *
 ****************************************************************************/
int Init_Vel_Profiles(FILE *fp,geologicunit_t *current){
    
    int numProfiles,numPoints, iProfile,iPoint;
    
    fscanf(fp," %d ", &numProfiles);
    current->numprofiles= numProfiles;    
    current->velprofiles= malloc(sizeof(velprofile_t)* numProfiles) ;
    
    for (iProfile=0; iProfile < numProfiles; iProfile++){
	fscanf(fp," %d ",&numPoints);
	current->velprofiles[iProfile].numpoints=numPoints;
	
	current->velprofiles[iProfile].depth = malloc(sizeof(float)*numPoints);
	current->velprofiles[iProfile].vp    = malloc(sizeof(float)*numPoints);
	current->velprofiles[iProfile].vs    = malloc(sizeof(float)*numPoints);
	current->velprofiles[iProfile].rho   = malloc(sizeof(float)*numPoints);
	
	for(iPoint=0; iPoint <numPoints;iPoint++)
	    fscanf(fp,"%f",&(current->velprofiles[iProfile].depth[iPoint]) );
	
	for(iPoint=0; iPoint < numPoints; iPoint++)
	    fscanf(fp,"%f",&(current->velprofiles[iProfile].vp[iPoint]) );
	
	for(iPoint=0; iPoint < numPoints; iPoint++)
	    fscanf(fp,"%f",&(current->velprofiles[iProfile].vs[iPoint]) );
	
	for(iPoint=0; iPoint < numPoints; iPoint++)
	    fscanf(fp,"%f",&(current->velprofiles[iProfile].rho[iPoint]) );
    }
    return numProfiles;
}

/***************************************************************************** 
 *                                                                           *
 *     Init_Vel_Profiles_Binary (speeds up the reading ):                    *
 *                                                                           *
 ****************************************************************************/
int Init_Vel_Profiles_Binary(FILE *fp,geologicunit_t *current){
    
    int32_t numProfiles,numPoints, iProfile,iPoint;
    
    fread(&numProfiles,sizeof(int),1,fp);
    current->numprofiles= numProfiles;    
    current->velprofiles= malloc(sizeof(velprofile_t)* numProfiles) ;
    
    for (iProfile=0; iProfile < numProfiles; iProfile++){
	fread(&numPoints,sizeof(int),1,fp);
	current->velprofiles[iProfile].numpoints=numPoints;
	
	current->velprofiles[iProfile].depth = malloc(sizeof(float)*numPoints);
	current->velprofiles[iProfile].vp    = malloc(sizeof(float)*numPoints);
	current->velprofiles[iProfile].vs    = malloc(sizeof(float)*numPoints);
	current->velprofiles[iProfile].rho   = malloc(sizeof(float)*numPoints);
	
	fread(&(current->velprofiles[iProfile].depth[0]),sizeof(float),numPoints,fp);
	fread(&(current->velprofiles[iProfile].vp[0])   ,sizeof(float),numPoints,fp);		   
	fread(&(current->velprofiles[iProfile].vs[0])   ,sizeof(float),numPoints,fp);
	fread(&(current->velprofiles[iProfile].rho[0])  ,sizeof(float),numPoints,fp);
    }
    
    return numProfiles;
}

/***************************************************************************** 
 *                                                                           *
 *    Init_Velocity_Distribution:INIT THE DISTRIBUTION OF VELOCITIES IN THE  *
 *                               GEOLOGIC UNIT, THE DISTRIBUTION IS DEFINED  *
 *                               IN TERMS OF A SURFACE WICH INDICATES THE    *
 *                               NUMBER OF PROFILE. THE PROFILES DESCRIBE THE* 
 *                               VARIATION WITH DEPTH OF THE VELOCITIES AND  *
 *                               DENSITIES, THEY CAN OR CANNOT BE NORMALIZED *
 *                               RELATIVE TO THE DEPTH OF THE GEOLOGIC UNIT  * 
 *    fpcontrol - file pointer to control.in                                 *
 *                                                                           *
 ****************************************************************************/
int Init_Velocity_Distribution(FILE *fpcontrol,geologicunit_t *current,int id){
    
    char fileVelProfiles[256],fileVelSurf[256];
    char parseaux[256];
    FILE *fpCurrent;
    double maxValProfile;
    int numProfiles;
    double maxVal;
    
    if(current->profiletype==0){ /*VEL DISTIB.SURF AND PROFILES 0 PROFILE*/
	
	sprintf(parseaux,"object_%d_znormalizationtype",id);
	Parse_Text(fpcontrol, parseaux,'i',
		   &(current->znormalizationtype));

	
	sprintf(fileVelProfiles,"%s/vel_profiles_%s.fun",(current->path),
		(current->geolunitname));

	fpCurrent = Fopen(fileVelProfiles,"r");       
	if(theIsBinary==1)
	numProfiles=Init_Vel_Profiles_Binary(fpCurrent, current);
	else
	numProfiles=Init_Vel_Profiles(fpCurrent, current);	
	fclose(fpCurrent);
	
       	sprintf(fileVelSurf,"%s/vel_profiles_%s.surf",
	current->path,current->geolunitname);
	fpCurrent = Fopen(fileVelSurf,"r");
	maxValProfile=Init_Surface(fpCurrent,&(current->veldistribsurface));
	fclose(fpCurrent);
	
	if(((int)maxValProfile) > (numProfiles-1)){
	    fprintf(stdout,"\n Error Init_Velocity_Distribution surf %d numProf %d maxVal%d FileName=%s",
	    id,numProfiles,(int)maxValProfile,(current->geolunitname)  );
	    exit(1);	    
	}
    }
    

    if(current->profiletype==1 ||current->profiletype==2 ){
	
	/* Common to current->profiletype==1 & 2*/
	sprintf(parseaux,"object_%d_znormalizationtype",id);
	Parse_Text(fpcontrol,parseaux,'i',&(current->znormalizationtype));
	
	sprintf(parseaux,"object_%d_velprofile_function",id);
	Parse_Text(fpcontrol,parseaux,'i',&(current->velprofilefunction));
	
	/* To constrain any polynomial approximation you use the same variables */	
	sprintf(parseaux,"object_%d_min_vs",id);
	Parse_Text(fpcontrol,parseaux, 'd',&(current->minvs));
	sprintf(parseaux,"object_%d_min_vp",id);
	Parse_Text(fpcontrol,parseaux, 'd',&(current->minvp));
	sprintf(parseaux,"object_%d_max_vs",id);
	Parse_Text(fpcontrol,parseaux, 'd',&(current->maxvs));
	sprintf(parseaux,"object_%d_max_vp",id);
	Parse_Text(fpcontrol,parseaux, 'd',&(current->maxvp));
	sprintf(parseaux,"object_%d_min_rho",id);
	Parse_Text(fpcontrol,parseaux, 'd',&(current->minrho));
	sprintf(parseaux,"object_%d_max_rho",id);
	Parse_Text(fpcontrol,parseaux, 'd',&(current->maxrho));

	/* weathering scheme */
	sprintf(parseaux,"object_%d_weatheringfactor",id);
	Parse_Text(fpcontrol,parseaux, 'd',&(current->weatheringfactor));
	sprintf(parseaux,"object_%d_weatheringexp",id);
	Parse_Text(fpcontrol,parseaux, 'd',&(current->weatheringexp));

	switch (current->velprofilefunction){	  
	case 0: /* Polynomial */
	    
	    sprintf(parseaux,"object_%d_polynomial_order",id);
	    Parse_Text(fpcontrol,parseaux, 'i',&((current->poly).order));
	    //	    Check_Error((current->poly).order <0,"Polynomial Order < 0" );
	    /* Initialize the polynomial */
	    (current->poly).coefficient=malloc(sizeof(double)*((current->poly).order+1));
	    (current->poly).exponent   =malloc(sizeof(double)*((current->poly).order+1));
	    
	    sprintf(parseaux,"object_%d_polynomial_coefficient",id);
	    Parse_Double_Array(fpcontrol,parseaux,(current->poly).order+1,1,&(current->poly).coefficient);
	    sprintf(parseaux,"object_%d_polynomial_power",id);
	    Parse_Double_Array(fpcontrol,parseaux,(current->poly).order+1,1,&(current->poly).exponent);
	    sprintf(parseaux,"object_%d_vpvs_ratio_a",id);
	    Parse_Text(fpcontrol,parseaux, 'd',&(current->vpvsa));
	    sprintf(parseaux,"object_%d_vpvs_ratio_b",id);
	    Parse_Text(fpcontrol,parseaux, 'd',&(current->vpvsb));
	    break;
	case 1:
	    sprintf(parseaux,"object_%d_exponent",id);
	    Parse_Text(fpcontrol,parseaux, 'd',&(current->exponentprofile));
	    break;
	}
	    

	/* Now read the differences in profiletype 2*/
	if (current->profiletype==2){	  
	    /* extrapolate from vs30 using a polynomial or power profile */
	    /* It only accepts linear variations of vp to vs ratio and it will  
	       assumed that it the relationship covers the entire depth of the 
	       unit */		
	    sprintf(parseaux,"object_%d_factorvs30",id);
	    Parse_Text(fpcontrol,parseaux, 'd',&(current->factorvs30));		
	    sprintf(fileVelSurf,"%s/surf_%s.surfvs30",
		  current->path,current->geolunitname);
	    fpCurrent = Fopen(fileVelSurf,"r");
	    maxVal=Init_Surface(fpCurrent,&(current->vs30surface));
	    fclose(fpCurrent);
	}	
    }        
    return 1;
}


/***************************************************************************
 *                                                                         *
 *                                                                         *
 * Init_Database: GENERAL INITIALIZATION FUNCTION.                         *
 *                                                                         *
 *     TYPES OF OBJECTS:                                                   *
 *             GEOLOGIC UNITS (0) --- THEY CAN BE DEFINED IN TERMS OF      *
 *                                    SURFACES                             *
 *                                                                         *
 *     visualcheckingoffon   -  0-off 1-on will display the unit it is     *
 *                              reading (off in parallel)                  *
 ***************************************************************************/
void Init_Database(database_t *current, char *databasePath,int visualcheckingoffon){
    
  char parseaux[256], parseauxfile[256], fileSurf[256], fileVelSurf[256];
  char controlin[256];
  FILE *fp, *fpCurrent;
  int iDb, howmanySurfaces=0, howmanyMeshes=0;
  
  sprintf(controlin,"%s/control.in",databasePath);
  fp=Fopen(controlin,"r");
  
  Parse_Text(fp, "x_min", 'd',&(current->xmin));
  Parse_Text(fp, "x_max", 'd',&(current->xmax));
  Parse_Text(fp, "y_min", 'd',&(current->ymin));
  Parse_Text(fp, "y_max", 'd',&(current->ymax));
  Parse_Text(fp, "depth_min", 'd',&(current->depthmin));
  Parse_Text(fp, "depth_max", 'd',&(current->depthmax));
  Parse_Text(fp, "datum", 'd',&(current->datum));  
  Parse_Text(fp, "number_of_objects", 'i',&(current->numberofobjects));
  Parse_Text(fp, "depth0", 'd',&(current->depth0));
  
  current->howmanygeologicunits = 0;
  current->howmanymeshes   = 0;
  
  current->priority = malloc( sizeof(int) * current->numberofobjects );
  current->subtype  = malloc( sizeof(int) * current->numberofobjects );
  
  for (iDb =0;iDb<current->numberofobjects; iDb++) {
    
    sprintf(parseaux, "object_%d_type",iDb);
    Parse_Text(fp, parseaux, 'i',&(current->subtype[iDb]));
    
    sprintf(parseaux, "object_%d_priority",iDb);
    Parse_Text(fp, parseaux, 'i',&(current->priority[iDb]));
    
    if( current->subtype[iDb] == 0 ) (current->howmanygeologicunits)+=1;
    if( current->subtype[iDb] == 1 ) (current->howmanymeshes)+=1;
  }
  
  current->meshes        = malloc(sizeof(meshdb_t) * current->howmanymeshes);
  current->geologicunits = malloc(sizeof(geologicunit_t)*
  				  current->howmanygeologicunits);

  current-> aretherebuildings=0;
  current-> aretherefoundations=0;


  /* GO THROUGH GEOLOGIC UNITS BOTTOM SURFACES AND VELOCITY DISTRIBUTION*/
  for (iDb =0;iDb<current->howmanygeologicunits; iDb++) {
    /* Visual checking */
    if(visualcheckingoffon==1)
      fprintf(stdout,"\n Reading Unit %d",iDb); fflush(stdout);
    
    sprintf(parseaux, "object_%d",iDb);
    Parse_Text(fp, parseaux, 's',parseaux);
    strcpy((current->geologicunits[iDb]).geolunitname,parseaux);

    (current->geologicunits[iDb]).path=databasePath;      
    sprintf(fileSurf,"%s/surf_%s.in",databasePath,
  	    (current->geologicunits[iDb]).geolunitname);
    fpCurrent = Fopen(fileSurf,"r");
    (current->geologicunits[iDb]).reflectorsurfaces.maxsurfaceval=
      Init_Surface(fpCurrent,&((current->geologicunits[iDb]).reflectorsurfaces));

    sprintf(parseaux,"object_%i_layer_type",iDb);
    Parse_Text(fp, parseaux,'i',
  	       &(current->geologicunits[iDb].layertype));
    
    sprintf(parseaux,"object_%i_interpolationscheme",iDb);
    Parse_Text(fp, parseaux,'i',
  	       &(current->geologicunits[iDb].interpolationscheme));

    sprintf(parseaux,"object_%i_factorsurface",iDb);
    Parse_Text(fp, parseaux,'d',
  	       &(current->geologicunits[iDb].factorsurface));


    if((current->geologicunits[iDb]).layertype==10)
      current-> aretherebuildings=1;
   
    if((current->geologicunits[iDb]).layertype==11)
      current-> aretherefoundations=1;

    /* 10 building & 11 foundation */
    if((current->geologicunits[iDb]).layertype<10){
      //      fprintf(stdout,"Correcting Unit %d",iDb); fflush(stdout);
      ShiftDatum_Surface(&((current->geologicunits[iDb]).reflectorsurfaces),current->datum);
    }
   
    sprintf(parseaux,"object_%i_isbottombedrock",iDb);
    Parse_Text(fp, parseaux,'i',
  	       &(current->geologicunits[iDb].isbottombedrock));
    /* (0) 3D DIST , (1) FUNCTION (2) VS30 surface */
    sprintf(parseaux,"object_%i_velprofile_type",iDb);
    Parse_Text(fp,parseaux,'i',&(current->geologicunits[iDb].profiletype));

    Init_Velocity_Distribution( fp, &(current->geologicunits[iDb]), iDb);
	
    sprintf(parseaux,"object_%d_isexcluded",iDb);
    Parse_Text(fp,parseaux,'i',&(current->geologicunits[iDb].isexcluded));
    fclose(fpCurrent);	 

    sprintf(parseaux,"object_%i_isbottomfreesurface",iDb);
    Parse_Text(fp, parseaux,'i',
  	       &(current->geologicunits[iDb].isbottomfreesurface));
  }  

  fprintf(stdout,"therearebuildings %d", current-> aretherebuildings);
  fprintf(stdout,"therearefoundations %d", current-> aretherefoundations);

  return;
}

