/**
 * Code: geodataq
 *
 *
 * Copyright (c) 2006 Leonardo Ramirez-Guzman
 * All rights reserved.  May not be used, modified, or copied
 * without permission.
 *
 * Contact:
 * Leonardo Ramirez-Guzman * 
 * lramirezg@iingen.unam.mx
 *
 *
 */

typedef struct nodedb_t{ 
  
    int    nodeid;
    double x,y,z;

}nodedb_t;

typedef struct polyg_t{ 
  
    int numnodes;  
    int *incidences;
    double vp, vs, rho, qs, qp;
    char name;
}polyg_t;


typedef struct meshdb_t{   
    int numnodes,numpolyg;    
    int  id;
    char name;
    nodedb_t *nodes;
    polyg_t *polygons;
    int numnodesinterpolated, *nodesinterpolated, *howmanynodes;
    double **nodesinterpcoordsx,**nodesinterpcoordsy;    
}meshdb_t;

typedef struct velprofile_t{
  int numpoints, asciiorbinary ;
    float *depth,*vs,*vp,*rho;     
    
}velprofile_t;


typedef struct surface_t{
    char *surfacename;

    /*    int profiletype;*/
  int nxgrid, nygrid;
  double *xgrid, *ygrid, **zgrid;
  double maxsurfaceval, minsurfaceval;
   
  /* double minvp,minvs,minrho,maxvp,maxvs,maxrho;*/
  /* int numprofiles;    */
  /* velprofile_t *velprofiles;*/    
}surface_t;


typedef struct polynomial_t{
    int order;
    double *coefficient;
    double *exponent;
}polynomial_t;


/* ***********************************************************************
 *
 *   struct geologicunit_t: The geologic reflectorsurfaces
 *                          Three types of distribution can be included
 *                          in the geologic unit
 *                           1) function based, see functions implemented
 *                           2) A 3D grid
 *                           3) A vs30 grid based distrubution (useful only
 *                              for geotechnical layers
 *
 ************************************************************************/
typedef struct geologicunit_t{
  
  char geolunitname[256];      /* The files describing the unit will have this name*/
  char *path;
  int profiletype, znormalizationtype;
  int numprofiles;    
  int isexcluded;
  int isbottombedrock;
  int velprofilefunction; 
  int isbottomfreesurface;
  int layertype;
  int interpolationscheme; /* 0-constant 1-linear*/
  
  double minvp,minvs,minrho,maxvp,maxvs,maxrho;
  double weatheringfactor,weatheringexp;
  double bedrockthreshold;
  double exponentprofile;
  double vpvsa,vpvsb;
  double factorvs30;
  double factorsurface;
  
  velprofile_t *velprofiles;
  polynomial_t poly;
  surface_t veldistribsurface, reflectorsurfaces; /* just the bottom for now*/
  surface_t vs30surface;

 

}geologicunit_t;


/* ***********************************************************************
 *
 *    struct database_t:
 *
 ************************************************************************/
typedef struct database_t{    
  char path[1024];
  int numberofobjects, howmanymeshes,howmanygeologicunits;
  int *priority;
  int *subtype;
  int aretherebuildings, aretherefoundations;
  meshdb_t *meshes;     /* 2D meshes used to construct background models*/
  geologicunit_t *geologicunits; /* Volumes*/
  double xmin,xmax, ymin,ymax,depthmin, depthmax; 
  double datum;  
  double depth0;
}database_t;
