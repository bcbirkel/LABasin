/*
 *
 * quakesource.c: source definition of earthquakes HERCULES version.
 *
 *  
 *          Input: physics.in, numerics.in source.in
 *
 * Optional input: slip.in, rake.in and genericfault.in.
 *
 *         Output: array of time dependant forces for each node
 *
 *
 * Definitions:
 *             
 *
 *              global coordinates - N-S   = X 
 *                                   E-W   = Y
 *                                   Depth = Z
 *
 *              domain coordinates - They could be the same as global
 *                                   if there is no azimutal rotation 
 *                                   of the rectangular prism. Other-
 *                                   wise the are the coordinates of 
 *                                   the prism, where the left edge is 
 *                                   x.
 *
 *                local coordinate - Local to the fault surface
 *                                     
 *
 *
 *
 * Notes: 
 *        Any modification comment should be made in the functions not at the
 *        begining of the code.
 *
 * Copyright (c) 2005 Leonardo Ramirez-Guzman
 * Supported by the Mexican National Council of Science and Technology 
 * Scholarship program.
 *
 * Created as part of the Quake project (tools for earthquake simulation).
 * There is no warranty whatsoever.
 * 
 * All rights reserved. May not be used, modified, or copied without permission.
 *
 * Contact:
 * Leonardo Ramirez Guzman
 * Civil and Environmental Engineering
 * Carnegie Mellon University
 * 5000 Forbes Avenue
 * Pittsburgh, PA 15213
 * lramirez@andrew.cmu.edu
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <mpi.h>

#include "octor.h"
#include "cvm.h"
#include "etree.h"
#include "psolve.h"
#include "quakesource.h"
#include "nrutila.h"
#include "util.h"


#define LINESIZE  512

#define ERROR  -100


#define    SRCNEEDED_MSG  100
#define NUMSRCNEEDED_MSG  101
#define       SRC_STATM0  102 
#define       SRC_STATFO  103 
#define       SRC_STATNO  104 
#define      SRC_STATMEM  105


/* convenience shortcut macros */
#define ABORT_PROGRAM(msg)			\
    do {					\
	if (NULL != msg) {			\
	    fputs (msg, stderr);		\
	}					\
	MPI_Abort (MPI_COMM_WORLD, ERROR);	\
	exit (1);				\
    } while (0)

#define EXITPROGRAM		ABORT_PROGRAM((void*)NULL)
#define WAIT			MPI_Barrier( MPI_COMM_WORLD )
#define ABORTEXIT               MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1)


/**
 * edata_t: Mesh element data fields 
 *
 */
typedef struct edata_t {
    float edgesize, Vp, Vs, rho;
} edata_t;

static int32_t myNumberOfForces = 0;
static int32_t myNumberOfNodesLoaded = 0;
static int myNumberOfCycles = 0;
static double myMemoryAllocated = 0;


/* Static global variables. */

/* Mesh related variables */
static octree_t *myOctree;
static mesh_t *myMesh;

static vector3D_t **myForces;

static double theSurfaceCornersLong[4], theSurfaceCornersLat[4];
static double theDomainX, theDomainY, theDomainZ;

static double theMinimumEdge;

/* MPI related variables */
static int32_t myID;
static int32_t theGroupSize;

/* General variables */   
static double theRegionOriginLatDeg, theRegionOriginLongDeg,theRegionAzimuthDeg;
static double theRegionDepthShallowM, theRegionLengthEastM;
static double theRegionLengthNorthM, theRegionDepthDeepM;
static double theDeltaT,theValidFreq; 

static int32_t        theNumberOfTimeSteps;
static int32_t        theSourceIsFiltered;
static source_type_t  theTypeOfSource;

/* Filter related */
static double theThresholdFrequency;
static int theNumberOfPoles;

/* Point and plane variables */
static int32_t theSourceFunctionType;
static double theAverageRisetimeSec, theMomentMagnitude, theMomentAmplitude;
static double theRickerTs, theRickerTp;

static int32_t theNumberOfTimeWindows;
static double *theWindowDelay;

/* Point variables */
static int theLonlatOrCartesian;
static double theHypocenterLatDeg, theHypocenterLongDeg, theHypocenterDepthM;
static double theSourceStrikeDeg, theSourceDipDeg, theSourceRakeDeg;

/* Plane variables */
static double theExtendedCellSizeAlongStrikeM, theExtendedCellSizeDownDipM;
static double theExtendedLatDeg, theExtendedLongDeg;
static double theExtendedDepthM,     theExtendedAlongStrikeDistanceM;
static double theExtendedDownDipDistanceM,   theExtendedHypocenterAlongStrikeM;
static double theExtendedHypocenterDownDipM, theExtendedAverageRuptureVelocity;
static double theExtendedStrikeDeg,      theExtendedDipDeg;        
static double theExtendedMinimumEdge;

/**
 * Time when the entire fault is no longer moving, i.e. slip velocity = 0
 * in each source point.
 */
static double theExtendedRaiseTime;

static vector3D_t theExtendedHypocenterCartesian;
static vector3D_t theExtendedHypLocCart;

static int32_t theExtendedColumns, theExtendedRows;
static double ***theSlipMatrix, ***theRakeMatrix;

static vector3D_t theExtendedCorner [4];

static int32_t numStepsNecessary;


/* Plane with kinks variables */
static double *theKinkLonArray, *theKinkLatArray;
static double *theTraceLengthAccumulated;
static double theTotalTraceLength;

static int theNumberOfKinks;
static vector3D_t *theTraceVectors;

vector3D_t  theKinksDomainCoord[6];


/* General source type variables SRFH */
static int theNumberOfPointSources;  
static int  *theSourceNt1Array;         
static double *theSourceLonArray,*theSourceLatArray,*theSourceDepthArray;
static double *theSourceAreaArray,*theSourceStrikeArray,*theSourceDipArray,
    *theSourceRakeArray,*theSourceSlipArray;
static double *theSourceTinitArray,*theSourceDtArray,**theSourceSlipFunArray ;

/* General source type variables GREEN */
static int  *theSourceNt1Array;         
static double *theSourceLonArray,*theSourceLatArray,*theSourceDepthArray;
static double *theSourceWXArray,*theSourceWYArray,*theSourceWZArray;
static double **theSourceSourceFunArray ;

/* Statistics variables and total moment ampitude*/
static double  myM0 =0;

/* I/O variables */
static int BUFFERSIZE =	    100e6;  /* change from 100 to 200 for shakeout */
int theForcesBufferSize;
static FILE *fpsourcetemp;
static char *forceinprocessor, *isforceinprocessor;
static int32_t *myForcesCycle;

/* Time statistic variables */
static double theForceGenerationTime = 0;
static double theConcatTime = 0;
static double theIOInitTime = 0; 


/*------------------------FUNCTIONS-------------------------------------------*/

#define twomb 1*1048576


/* return the number of MB available */

int compute_available_mem(){

    int64_t memoryTrial =0;
    int availableMem=0;
    int64_t *ret;
    
    memoryTrial = twomb;

    /* check if I can alloc memoryTrial */
    availableMem = 0;
    ret = malloc(memoryTrial);

    if(ret == NULL) 
	return availableMem;
    
    while (ret != NULL  ){
	
	free(ret);
	
	availableMem += 1 ;
	
	memoryTrial +=twomb;
	
	ret = malloc(memoryTrial);
	
	if( ret == NULL )
	    return availableMem;
	
    }
    
    
}
	


static int read_planewithkinks (FILE *fp);


static void compute_point_source_strike (ptsrc_t* ps);


/**
 * parsedarray: Parse array of doubles .
 *
 * - return 0 if OK, -1 on error
 *
 */
int parsedarray( FILE *fp, const char *querystring,int size, double *array){

    rewind (fp);      
    int found=0;
    int iSize;
  
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
	    for ( iSize = 0; iSize < size; iSize++ )
		if ( fscanf ( fp," %lf ", &(array[iSize] ))  == 0 ){  
		    fprintf(fp,"\nUnable to read %s",querystring  );
		    return -1;
		}              
	    return 0;       
	}
    }  
    return -1;
}


/*
 * interpolate_linear: interpolate linearly a function 
 *                     if the time is larger than the one supported by the function
 *                     the last value will be assigned
 *
 */

static double  interpolate_linear ( double time,int numsamples, double samplingtime, 
				    double *discretefunction ){
    double maxtime, m, b;

    int interval;

    maxtime = (numsamples-1)*samplingtime;

    if ( time >= maxtime ){    
	return discretefunction[numsamples-1];
    }
    else{

	/* locate the interval */
	interval = floor(time/samplingtime);
    
	/* y = mx +b */
	m =  discretefunction[interval+1] - discretefunction[interval] ;
	m =  m/samplingtime;

	b =  discretefunction[interval] - m * interval * samplingtime;

	return m*time+b;    
    }    
}


/**
 *
 * compute_source_function: compute the slip source function:
 *  
 *        Input:
 *                pointSource - pointer to a point source. All
 *                              information required has to be
 *                              loaded in it.
 *        Output: void
 *
 *        Notes: 
 *               numStepsNecessary  is the number of steps required 
 *               to reach the maximum value of the slip.
 *
 */
void
compute_source_function (ptsrc_t* pointSource)
{
    double  decay, t1, T, tao;
    int     iTime;
  
    /*  printf("\n sourcefunctiontype = %d ",  pointSource->sourceFunctionType);
	exit(1);*/

    /* for ( iTime=0; iTime < pointSource->numberOfTimeSteps; iTime++) { */  
    for (iTime=0; iTime < numStepsNecessary; iTime++) {
    
	T = pointSource->dt * iTime;
	tao= T / pointSource->T0;     
    
	if (pointSource->delayTime < T) {
	    T =  T -pointSource->delayTime;
      
	    switch (pointSource->sourceFunctionType) {
	
	    case RAMP:	
		decay = (T < pointSource->T0) ? T / pointSource->T0 : 1;	
		break;
	
	    case SINE:
		decay = ((T < pointSource->T0) 
			 ? ( T / pointSource->T0
			     - sin ( 2 * PI * T / pointSource->T0) / PI / 2)
			 : 1);
		break;
	
	    case QUADRATIC: 	
		if ( T < pointSource->T0 / 2 ) {
		    decay = 2 * pow (T / pointSource->T0, 2);
		}
		else if (T <= pointSource->T0) {
		    decay = (-2 * pow (T / pointSource->T0, 2)
			     + 4 * T / pointSource->T0 - 1);
		}
		else {
		    decay = 1;
		}
		break;
	
	    case RICKER:  /* Ricker's function */
		t1 = pow ((T - pointSource->Ts) * PI / pointSource->Tp, 2);
		decay = (t1 - 0.5) * exp (-t1);
		break;
	
	    case EXPONENTIAL:
		/* decay = 1-(1 + T/pointSource->T0)*exp(-T/pointSource->T0); */
		decay = 1 - (1 + tao) * exp(-tao);
		break;

	    case DISCRETE:	      
		decay = interpolate_linear(T, pointSource->nt1,
					   pointSource->dtfunction, 
					   pointSource->slipfundiscrete);	

		break;
	
	    default:
		fprintf (stderr, "Unknown source type %d\n",
			 pointSource->sourceFunctionType);
		ABORT_PROGRAM ("Cannot compute source function");
		return;
	    };
	}     
	else {
	    decay =0;
	}        
  
	pointSource->displacement[ iTime ] += decay * pointSource->maxSlip;    

    }
   
    return;
  
}


/**
 *
 * compute_source_function_green: compute the source function:
 *  
 *        Input:
 *                pointSource - pointer to a point source. All
 *                              information required has to be
 *                              loaded in it.
 *        Output: void
 *
 *        Notes: 
 *               numStepsNecessary  is the number of steps required 
 *               to reach the maximum value of the slip.
 *
 */
void
compute_source_function_green (ptsrc_t* pointSource)
{
    double  decay, t1, T, tao;
    int     iTime;
  
    /*  printf("\n sourcefunctiontype = %d ",  pointSource->sourceFunctionType);
	exit(1);*/

    /* for ( iTime=0; iTime < pointSource->numberOfTimeSteps; iTime++) { */  
    for (iTime=0; iTime < numStepsNecessary; iTime++) {
    
	T = pointSource->dt * iTime;
	tao= T / pointSource->T0;     
    
	if (pointSource->delayTime < T) {
	    T =  T -pointSource->delayTime;
      
	    switch (pointSource->sourceFunctionType) {
	
	    case RAMP:	
		decay = (T < pointSource->T0) ? T / pointSource->T0 : 1;	
		break;
	
	    case SINE:
		decay = ((T < pointSource->T0) 
			 ? ( T / pointSource->T0
			     - sin ( 2 * PI * T / pointSource->T0) / PI / 2)
			 : 1);
		break;
	
	    case QUADRATIC: 	
		if ( T < pointSource->T0 / 2 ) {
		    decay = 2 * pow (T / pointSource->T0, 2);
		}
		else if (T <= pointSource->T0) {
		    decay = (-2 * pow (T / pointSource->T0, 2)
			     + 4 * T / pointSource->T0 - 1);
		}
		else {
		    decay = 1;
		}
		break;
	
	    case RICKER:  /* Ricker's function */
		t1 = pow ((T - pointSource->Ts) * PI / pointSource->Tp, 2);
		decay = (t1 - 0.5) * exp (-t1);
		break;
	
	    case EXPONENTIAL:
		/* decay = 1-(1 + T/pointSource->T0)*exp(-T/pointSource->T0); */
		decay = 1 - (1 + tao) * exp(-tao);
		break;

	    case DISCRETE:	      
		decay = interpolate_linear(T, pointSource->nt1,
					   pointSource->dtfunction, 
					   pointSource->sourcefunctiondiscrete );	

		break;
	
	    default:
		fprintf (stderr, "Unknown source type %d\n",
			 pointSource->sourceFunctionType);
		ABORT_PROGRAM ("Cannot compute source function");
		return;
	    };
	}     
	else {
	    decay =0;
	}        
  
	pointSource->displacement[ iTime ] += decay;    

    }
   
    return;
  
}




/**
 * search_source: search a point in the domain of the local mesh
 * 
 *   input: coordinates 
 *  output: 0 fail 1 success
 *
 */ 
static int32_t search_point( vector3D_t point, octant_t **octant  ){
  
    tick_t  xTick, yTick, zTick;
  
    xTick = point.x[0] / myMesh->ticksize;
    yTick = point.x[1] / myMesh->ticksize;
    zTick = point.x[2] / myMesh->ticksize;
  
    *octant = octor_searchoctant(myOctree, xTick, yTick, zTick, 
				 PIXELLEVEL, AGGREGATE_SEARCH);
  
    if ( (*octant == NULL) || ((*octant)->where == REMOTE) )     
	return 0;
    
    return 1;
  
}

/*
 *
 * source_initnodalforce: Init the wighting force vector for an 
 *                        element due to a point source.
 *
 */
static void source_initnodalforce ( ptsrc_t *sourcePtr ){
  
    /* Normal vector n, tangent vector t, moment tensor v*/
    double n[3], t[3], v[3][3];
  
    /* Auxiliar array to handle shapefunctions in a loop */
    double  xi[3][8]={ {-1,  1, -1,  1, -1,  1, -1, 1} ,
		       {-1, -1,  1,  1, -1, -1,  1, 1} ,
		       {-1, -1, -1, -1,  1,  1,  1, 1} };
  
  
    /* various convienient variables */
    int j, k;
    double dx, dy, dz;
    double s = sourcePtr->strike / 180.0 * PI;
    double d = sourcePtr->dip / 180.0 * PI;
    double r = sourcePtr->rake / 180.0 * PI;
    double x = sourcePtr->x;
    double y = sourcePtr->y;
    double z = sourcePtr->z;
    double h = sourcePtr->edgesize;
    double hcube = h * h * h;
  
    /* the fault normal unit vector */
    n[0] = - sin(s) * sin(d);
    n[1] =   cos(s) * sin(d);
    n[2] = - cos(d);
  
    /* the fault tangential unit vector */
    t[0] =   cos(r) * sin(PI / 2 - s) + sin(r) * sin(s) * cos(d);
    t[1] =   cos(r) * sin(s) - sin(r) * cos(s) * cos(d);
    t[2] = - sin(r) * sin(d);
  
    for ( j = 0; j < 3; j++)
	for ( k = 0; k < 3; k++)
	    v[j][k] = n[j] * t[k] + n[k] * t[j];
  

    /* calculate equivalent force on each node */
    for ( j = 0; j < 8 ; j++){
	dx= (2 * xi[0][j]) * (h + 2 * xi[1][j] * y) * (h + 2 * xi[2][j] * z) 
	    / (8 * hcube) ;
    
	dy= (2 * xi[1][j]) * (h + 2 * xi[2][j] * z) * (h + 2 * xi[0][j] * x)
	    / (8 * hcube);
      
	dz= (2 * xi[2][j]) * (h + 2 * xi[0][j] * x) * (h + 2 * xi[1][j] * y)
	    / (8 * hcube);
      
	sourcePtr->nodalForce[j][0] = 
	    v[0][0] * dx + v[0][1] * dy + v[0][2] * dz;
	sourcePtr->nodalForce[j][1] = 
	    v[1][0] * dx + v[1][1] * dy + v[1][2] * dz;
	sourcePtr->nodalForce[j][2] = 
	    v[2][0] * dx + v[2][1] * dy + v[2][2] * dz;
    }
  
    return;
  
}



/*
 *
 * source_initnodalforce_greensfunction: Init the wighting force vector for an 
 *                        element due to a point force. 
 *
 */
static void source_initnodalforce_greensfunction ( ptsrc_t *sourcePtr ){
  
    /* Auxiliar array to handle shapefunctions in a loop */
    double  xi[3][8]={ {-1,  1, -1,  1, -1,  1, -1, 1} ,
		       {-1, -1,  1,  1, -1, -1,  1, 1} ,
		       {-1, -1, -1, -1,  1,  1,  1, 1} };
  
    /* various convienient variables */
    int j;
    double fix,fiy,fiz;
    double fx = sourcePtr->weightX;
    double fy = sourcePtr->weightY;
    double fz = sourcePtr->weightZ;
    double x = sourcePtr->x;
    double y = sourcePtr->y;
    double z = sourcePtr->z;
    double h = sourcePtr->edgesize;
    double hcube = h * h * h;

    double sourceFx=0,sourceFy=0,sourceFz=0; //to make sure it is one

  
    /* calculate equivalent force on each node */
    for ( j = 0; j < 8 ; j++){
	fix= (h + 2 * xi[0][j] * x ) * (h + 2 * xi[1][j] *  y ) * ( h + 2 * xi[2][j] * z) 
	    / (8 * hcube) ;
    
	fiy=  (h + 2 * xi[0][j] * x) * (h + 2 * xi[1][j] *  y ) * ( h + 2 * xi[2][j] * z)
	    / (8 * hcube);
      
	fiz=  (h + 2 * xi[0][j] * x) * (h + 2 * xi[1][j] *  y ) * ( h + 2 * xi[2][j] * z)
	    / (8 * hcube);
      
	sourcePtr->nodalForce[j][0] = fix*fx;
	sourcePtr->nodalForce[j][1] = fiy*fy;
	sourcePtr->nodalForce[j][2] = fiz*fz;

	sourceFx+=fix*fx;
	sourceFy+=fiy*fy;
	sourceFz+=fiz*fz;
    }

    fprintf(stdout,"\n Force X= %lf , Force Y= %lf, Force Z= %lf h=%lf ",sourceFx,sourceFy,sourceFz,h);
  
    return;
  
}


/*
 *  compute_initial_time: computes the time when the rupture is initiated in
 *                        the given station. It assumes a homogeneous radial 
 *                        dependent rupture time.
 *
 *                
 */
double compute_initial_time(vector3D_t station, vector3D_t hypocenter, 
			    double rupturevelocity ){
  

    double travelTime;
  
    int32_t i;
  
    travelTime=0;
     
    for ( i = 0; i < 3;  i++ )
	travelTime += pow( station.x[ i ] -
			   hypocenter.x[ i ], 2);
     
    travelTime = sqrt( travelTime);     
    travelTime = travelTime / rupturevelocity;
  
    return travelTime;    
  
}
  
/* 
 * compute_cartesian_coords:
 * 
 *  input: Lat(Deg),Lon(Deg),Depth
 * output: vector3D_t pointInCartesian
 *
 */
static vector3D_t compute_cartesian_coords ( double lat, double lon, 
					     double depth){
  
    vector3D_t pointInCartesian;
    
    pointInCartesian.x[0] = ( lat - theRegionOriginLatDeg ) 
	* DIST1LAT;   
    pointInCartesian.x[1] = ( lon - theRegionOriginLongDeg ) 
	* DIST1LON;   
    pointInCartesian.x[2] = depth - theRegionDepthShallowM ;
  
    return pointInCartesian;   
}


/* 
 * compute_local_coords: computes the coordinates of a station in the fault 
 *                       plane related to the origin of the fault
 * 
 *  input: point
 *         totalNodes
 *         gridx, gridy
 *
 * output: vector3D_t pointInCartesian
 *
 */
vector3D_t compute_local_coords(int32_t point, int32_t totalNodes,
				double *gridx, double *gridy ){
  
    vector3D_t pointVector;
  
    pointVector.x[ 0 ] = gridx [ point % totalNodes ]; 
    pointVector.x[ 1 ] = gridy [ ( int32_t ) point / totalNodes ];
    pointVector.x[ 2 ] =  0;
  
    return pointVector;   
}

/*
 *
 *  source_compute_print_stat():
 *
 */
static int32_t source_compute_print_stat(){
  
    /* Compute the total M0 and Mw */
    if ( myID == 0 ){   
    
	int32_t fromWhom;     
	int32_t rcvNumberOfForces, rcvNumberOfNodesLoaded;
	int rcvNumberOfCycles;
     
	double rcvMemoryAllocated,theTotalM0=0,rcvM0, Mw;      
	double theTotalMemoryAllocated =0 ;
    
	int iCorner, iCoord;
        
	MPI_Status status;
	theTotalM0=myM0;
	theTotalMemoryAllocated += myMemoryAllocated;

	if ( theTypeOfSource == PLANE ) {
	    fprintf(stdout,"\n\n Extended Fault Information \n ");    
	    fprintf(stdout,"\n Fault's Corners Global Coordinates (X,Y,Z)\n");
    
	    for ( iCorner = 0; iCorner < 4; iCorner++){       
		fprintf(stdout,"\n");
      
		for( iCoord = 0; iCoord < 3; iCoord++) 	
		    fprintf(stdout," %f ",theExtendedCorner[iCorner].x[iCoord]);       
	    }     	
	}
     
	fprintf ( stdout,"\n\n Process       M0   myNumOfForces" );     
	fprintf ( stdout,"     myCycles  myNumNodesLoaded  myMemAlloc");  
	if( myNumberOfForces != 0)
	    fprintf ( stdout,"\n    %-4d   %-10e   %-10d    %-5d    %-10d   %e",myID,myM0,
		      myNumberOfForces,myNumberOfCycles, myNumberOfNodesLoaded, myMemoryAllocated);
    
	for ( fromWhom=1; fromWhom< theGroupSize ; fromWhom++ ){

	    MPI_Recv ( &rcvNumberOfCycles, 1,   MPI_INT,fromWhom , 
		       SRC_STATM0,  MPI_COMM_WORLD, &status);
	    MPI_Recv ( &rcvM0, 1,   MPI_DOUBLE,fromWhom , 
		       SRC_STATM0,  MPI_COMM_WORLD, &status);       
	    MPI_Recv ( &rcvNumberOfForces, 1, MPI_INT,fromWhom , 
		       SRC_STATFO,  MPI_COMM_WORLD, &status);       
	    MPI_Recv ( &rcvNumberOfNodesLoaded, 1, MPI_INT,fromWhom ,
		       SRC_STATNO,  MPI_COMM_WORLD, &status);       
	    MPI_Recv ( &rcvMemoryAllocated, 1, MPI_DOUBLE,fromWhom , 
		       SRC_STATMEM, MPI_COMM_WORLD, &status);         	            
      
	    theTotalM0+=rcvM0;
	    theTotalMemoryAllocated += rcvMemoryAllocated;

	    if( rcvNumberOfForces != 0)	 
		fprintf ( stdout,"\n    %-4d   %-10e   %-10d    %-5d    %-10d   %e" ,fromWhom,rcvM0,
			  rcvNumberOfForces,rcvNumberOfCycles, rcvNumberOfNodesLoaded,
			  rcvMemoryAllocated);
    
	}

	if ( theSourceIsFiltered == 1 ) 
	    fprintf( stdout,"\n  Filtered  Threshold frequency = %lf Poles =%d \n", 
		     theThresholdFrequency, theNumberOfPoles);
	else 
	    fprintf(stdout,"\n  Non Filtered \n");   

	Mw=(log10 ( theTotalM0*pow(10,7))/1.5)-10.73;

	fprintf ( stdout,"\n   M0  =  %e Nm  \n Mw = %f",theTotalM0, Mw );
	fprintf ( stdout,"\n Total Memory Allocated = %e\n", 
		  theTotalMemoryAllocated );
    }

    else {
	MPI_Ssend (&myNumberOfCycles, 1, MPI_INT, 0, SRC_STATM0, MPI_COMM_WORLD);
	MPI_Ssend (&myM0,          1, MPI_DOUBLE, 0, SRC_STATM0, MPI_COMM_WORLD);
	MPI_Ssend (&myNumberOfForces, 1, MPI_INT, 0, SRC_STATFO, MPI_COMM_WORLD);
	MPI_Ssend (&myNumberOfNodesLoaded,1,MPI_INT,0,SRC_STATNO,MPI_COMM_WORLD);
	MPI_Ssend (&myMemoryAllocated,1,MPI_DOUBLE,0,SRC_STATMEM,MPI_COMM_WORLD);
    }

    return 1;
}



/*
 *  This computes an in-place complex-to-complex FFT 
 *  x and y are the real and imaginary arrays of 2^m points.
 *  dir =  1 gives forward transform
 *  dir = -1 gives reverse transform 
 */
int FFT(short int dir,long m,double *x,double *y){
    long n,i,i1,j,k,i2,l,l1,l2;
    double c1,c2,tx,ty,t1,t2,u1,u2,z;

    /* Calculate the number of points */
    n = 1;
    for (i=0;i<m;i++) 
	n *= 2;

    /* Do the bit reversal */
    i2 = n >> 1;
    j = 0;
    for (i=0;i<n-1;i++) {
	if (i < j) {
	    tx = x[i];
	    ty = y[i];
	    x[i] = x[j];
	    y[i] = y[j];
	    x[j] = tx;
	    y[j] = ty;
	}
	k = i2;
	while (k <= j) {
	    j -= k;
	    k >>= 1;
	}
	j += k;
    }

    /* Compute the FFT */
    c1 = -1.0; 
    c2 = 0.0;
    l2 = 1;
    for (l=0;l<m;l++) {
	l1 = l2;
	l2 <<= 1;
	u1 = 1.0; 
	u2 = 0.0;
	for (j=0;j<l1;j++) {
	    for (i=j;i<n;i+=l2) {
		i1 = i + l1;
		t1 = u1 * x[i1] - u2 * y[i1];
		t2 = u1 * y[i1] + u2 * x[i1];
		x[i1] = x[i] - t1; 
		y[i1] = y[i] - t2;
		x[i] += t1;
		y[i] += t2;
	    }
	    z =  u1 * c1 - u2 * c2;
	    u2 = u1 * c2 + u2 * c1;
	    u1 = z;
	}
	c2 = sqrt((1.0 - c1) / 2.0);
	if (dir == 1) 
	    c2 = -c2;
	c1 = sqrt((1.0 + c1) / 2.0);
    }

    /* Scaling for forward transform */
    if (dir == 1) {
	for (i=0;i<n;i++) {
	    x[i] /= n;
	    y[i] /= n;
	}
    }
   
    return 1;
}


/**
 * Compute Butterworth filter value for a given frequency.
 *
 * \param f0	Filter's threshold frequency.
 * \param m	Filter's order.
 * \param f	Frequency for which the coefficient is to be computed.
 */
static double
bw_filter_coeffient (double f0, int m, double f)
{
    return sqrt ( 1.0 / (pow (f / f0, 2 * m) + 1) );
}


static int
print_filter(
	     const char*  filename,
	     unsigned int signal_length,
	     double       sampling_frequency,
	     double       threshold_frequency,
	     unsigned int m
	     )
{
    unsigned int filter_length, i;
    double	 delta_frequency;
    double	 f0 = threshold_frequency;	/* alias for threshold freq */

    FILE* fp = fopen (filename, "w");

    if (NULL == fp) {
	return -1;
    }

    /* n = ( int ) ( log ( signal_size ) / log( 2 ) ) + 1 ;  */
    /* filter_length = pow ( 2 , n ); */

    /* rely on integer arithmetic for this to work */
    filter_length = (signal_length / 2) + 1;

    fputs ("0  1\n", fp);

    delta_frequency = sampling_frequency / filter_length;

    for (i = 1; i <= filter_length; i++) {
	double f      = delta_frequency * i;	/* current frequency */
	double filter = bw_filter_coeffient (f0, m, f);
	fprintf (fp, "%lf %lf\n", f, filter);
    }

    fflush (fp);
    fclose (fp);

    return 0;
}


/**
 * Print the slip function for the "average" source function.
 */
static int
print_slip_function (const char* filename)
{
    int i;
    double* p_disp;	/* pointer to the displacement array */
    ptsrc_t src;	/* point source */

    FILE* fp = fopen (filename, "w");

    if (NULL == fp) {
	perror ("Could not open file");
	fprintf (stderr, "Filename: %s\n", filename);
	return -1;
    }

    /* Compute the "average" source function */
    src.dt		   = theDeltaT;
    src.numberOfTimeSteps  = theNumberOfTimeSteps;
    src.delayTime	   = 0;
    src.sourceFunctionType = theSourceFunctionType;
    src.T0		   = theAverageRisetimeSec;
    src.maxSlip		   = 1;

    p_disp = (double*) calloc (sizeof (double), theNumberOfTimeSteps);

    if (NULL == p_disp) { /* memory allocation failed */
	fprintf (stderr, "Could not allocate memory for point source "
		 "displacement\nTrying to allocate %d doubles, size = %zu\n",
		 theNumberOfTimeSteps, theNumberOfTimeSteps * sizeof (double));
	fclose (fp);
	return -1;
    }

    src.displacement = p_disp;
    compute_source_function (&src);

    for (i = 0; i < theNumberOfTimeSteps; i++) {
	fprintf (fp, "%lf %lf\n", i * src.dt, src.displacement[i]);
    }

    free (p_disp);
    fclose (fp);

    return 0;
}


/**
 * Print the filter values (for the given parameters) and the values
 * of the "average" source function.
 */
int
print_filter_and_signal(
			int    signalsize,
			double samplingfrequency,
			double thresholdfrequency,
			int    m
			)
{
    int r1, r2;

    r1 = print_filter ("filter.dat", signalsize, samplingfrequency,
		       thresholdfrequency, m);
    if (0 != r1) {
	ABORT_PROGRAM ("Could not print filter data");
	return -1;
    }

    r2 = print_slip_function ("slipfunc.dat");
    if (0 != r2) {
	ABORT_PROGRAM ("Could not print source slip function");
	return -1;
    }

    return 1;
}


/**
 * print_myForce_filepercycle:
 *
 */
static void print_myForces_filepercycle(char *filenameroot, int cycle) {
    
    int32_t lnid, i, j;        
    char forcefile[256];  
    int32_t numberOfNodes=0,tmpNumberOfNodes=0;
    FILE *fp;
    
    sprintf(forcefile, "%s.%d", filenameroot, cycle);
  
    fp = fopen(forcefile,"w");

    for ( lnid = 0; lnid <  myMesh->nharbored; lnid++){      
	if ( myForcesCycle [ lnid ] == cycle ){ 
	    numberOfNodes+=1;
	}
    }

    /* Checking */
    tmpNumberOfNodes=0;
    for ( lnid = 0; lnid <  myMesh->nharbored; lnid++){      
	if ( myForces[ lnid ] != NULL ){ 
	    tmpNumberOfNodes+=1;
	}
    }

    /* to debug */
    if(tmpNumberOfNodes != numberOfNodes){
	fprintf(stdout,"\n n =%d l=%d", tmpNumberOfNodes,numberOfNodes);

	EXITPROGRAM;
    }

  
    fwrite( &numberOfNodes, sizeof(int32_t),1,fp);

    for ( j = 0; j < theNumberOfTimeSteps; j++ ) {
	for ( lnid = 0; lnid <  myMesh->nharbored; lnid++){      
	    if ( myForcesCycle [ lnid ] ==  cycle ){ 
		for ( i = 0; i < 3;  i++ ){			
		    fwrite( &(myForces[lnid][j].x[i]) , sizeof(double), 1,fp);
		}
	    }
	}
    }
  
    fclose(fp);  
    return;
}

/**
 * print_myForce:
 *
 */
static void print_myForces(FILE *fp) {

    int32_t lnid, i, j;      
    
    for ( lnid = 0; lnid <  myMesh->nharbored; lnid++){    
	if ( myForces [ lnid ] !=  NULL ){      
	    for ( i = 0; i < 3;  i++ ){	
		fprintf( fp,"\n" );	
		for ( j = 0; j < theNumberOfTimeSteps; j++ )
		    fprintf ( fp," %f ", myForces [ lnid ] [ j ].x[ i ] ); 	   
	    }
	}
    }

    return;

}


/**
 * print_myForce:
 *
 */
static void print_myForces_transposed(FILE *fp) {

    int32_t lnid, iDir, j;      
  
    for ( j = 0; j < theNumberOfTimeSteps; j++ ) {
	for ( lnid = 0; lnid <  myMesh->nharbored; lnid++){    
	    if ( myForces [ lnid ] !=  NULL ){      
		for ( iDir = 0; iDir < 3;  iDir++ ){
		    fwrite( &(myForces[lnid][j].x[iDir]), sizeof(double), 1,fp);
		}
	    }
	}
    }

    return;

}



/**
 * free_myForce:
 *
 */
static void free_myForces( ) {

    int32_t lnid, i;      
  
    for ( lnid = 0; lnid <  myMesh->nharbored; lnid++){
    
	if ( myForces [ lnid ] !=  NULL ){
	    free(myForces[ lnid ]);      
	}      	 
    }
    /* Initialize myForces with NULL */  
    for ( i = 0; i < myMesh->nharbored; i++)
	myForces[i] = NULL; 
    return;

}

 
/*
 *  Filters a signal we add zeros such that newsize is 2^n 
 *
 */

int FilterSignal ( double *signal, int signalsize,
		   double samplingfrequency, double thresholdfrequency, 
		   int m ) {
  
  
    double *filter = NULL, *signalReal = NULL,  *signalImaginary = NULL;  
    double *signalFilteredReal = NULL , *signalFilteredImaginary = NULL;  
    double frequency;
  
    int i, n,  newSize;
  
    n = ( int ) ( log ( signalsize ) / log( 2 ) ) + 2 ;  
    newSize = pow ( 2 , n );
  
    /* Allocate and check */
    signalReal =              (double *) malloc( sizeof(double) * newSize);  
    signalImaginary =         (double *) malloc( sizeof(double) * newSize);
    signalFilteredReal =      (double *) malloc( sizeof(double) * newSize);
    signalFilteredImaginary = (double *) malloc( sizeof(double) * newSize);  
    filter =                  (double *) malloc( sizeof(double) * newSize); 
    
    if ( !signalReal || !signalImaginary || !signalFilteredReal || 
	 !signalFilteredImaginary || !filter ) {    
	perror("Failed to allocate space for the filter");
	MPI_Abort(MPI_COMM_WORLD, ERROR);
	exit(1);        
    }
  

    /* Init imaginary with zeroes and move signal to signalReal*/
    /* We will filter the derivative of the signal, notice that this is not 
       exactly correct but it allows you to introduce zeroes without a big 
       jump in the slip function */  
    /* Take derivative of the signal onsistent with the 2 order polynomial 
       you are using for the time interpolation */ 
  
    for ( i = 0 ; i < signalsize; i++ ){
    
	if ( i == 0 ){     
	    signalReal [ i ] = .5 * samplingfrequency * ( -3 * signal[i] +
							  4 * signal[i+1] -
							  1 * signal[i+2]);
	}      
	else if( i == signalsize-1)     
	    signalReal [ i ] = .5 * samplingfrequency * (      signal[i-2] -
							       4*signal[i-1] + 
							       3*signal[i  ]);   
	else     
	    signalReal [ i ] = .5 * samplingfrequency * (     -signal[i-1] +
							      signal[i+1]);         
    }
 
    for ( i = 0 ; i < newSize; i++)
	signalImaginary [ i ] = 0;

    for ( i = signalsize ; i < newSize; i++)    
	signalReal [ i ] = 0;

    FFT( 1,( long ) n, signalReal , signalImaginary);
  
    /* Calculate the filter, multiply and obtain the conjugate*/

    /* Frequency=0 */

    frequency = 0; 
    filter [ 0 ] = 1;  
    signalFilteredReal [ 0 ] = filter [ 0 ] * signalReal [ 0 ];  
    signalFilteredImaginary [ 0 ] = filter [ 0 ] * signalImaginary [ 0 ];   
  
    for ( i = 1; i <= newSize / 2 ; i++ ){    
	frequency = samplingfrequency * i / newSize;    
	filter [ i ]= 1 / ( 1 + pow ( ( frequency / thresholdfrequency ),
				      2 * m ));    
	filter [ i ]=pow( filter [ i ] , .5 );   
	signalFilteredReal [ i ] = filter [ i ] * signalReal [ i ];      
	signalFilteredImaginary [ i ] = filter [ i ] * signalImaginary [ i ];      
	signalFilteredReal [ newSize-i ] = signalFilteredReal [ i ];      
	signalFilteredImaginary [ newSize-i ] = -signalFilteredImaginary [ i ];             
    }
    
    FFT( -1, ( long ) n, signalFilteredReal, signalFilteredImaginary);
    
    /* Move signalFilteredReal to signal and integrate*/
  
    for ( i = 0 ; i < signalsize; i++ ){    
	if ( i == 0 )
	    signal [ i ] = 0;
	else if ( i == 1 )
	    signal [ i ] = .5 * ( 1/samplingfrequency ) *
		( signalFilteredReal [ i-1 ]+
		  signalFilteredReal [ i ]    );
	else      
	
	    signal [ i ] = signal [i-1] + .5*( 1/samplingfrequency ) * 
		( signalFilteredReal [ i - 1 ] + signalFilteredReal [ i     ] ); 
    }
    
    free(signalReal);  
    free(signalImaginary);  
    free(signalFilteredReal);
    free(signalFilteredImaginary);  
    free(filter);
  
    return 1;
  
}

  
/**
 * filter_myForce:
 *
 */
static void filter_myForce() {

    int32_t lnid, i, j;      

    double *signal;
  
    signal = ( double * ) malloc( sizeof ( double ) * theNumberOfTimeSteps );

    if ( !signal ) {     
	perror("Failed to allocate space for the signal");
	MPI_Abort(MPI_COMM_WORLD, ERROR);
	exit(1);     
    }    

    for ( lnid = 0; lnid <  myMesh->nharbored; lnid++ ) {     
	if ( myForces [ lnid ] !=  NULL ){       
	    for ( i = 0; i < 3;  i++ ) {      	
		for ( j = 0; j < theNumberOfTimeSteps; j++ )
		    signal [ j ] = myForces [ lnid ] [ j ].x[ i ]; 	
		FilterSignal(signal, theNumberOfTimeSteps,
			     1/theDeltaT,theThresholdFrequency, theNumberOfPoles);
		for ( j = 0; j < theNumberOfTimeSteps; j++)
		    myForces [ lnid ] [ j ].x [ i ] = signal [ j ]; 	     
	    }
	}     
    }   
    free(signal);

    return;
}
   

/**
 * load_myForces_with_point_source: computes the force vector for a 
 *                                  point source and updates myForeces
 * 
 *   input: octant where the source is located
 *          the pointer of the displacement array (time dependant)
 *  output: 0 fails 1 success
 *
 */ 
static int load_myForces_with_point_source (octant_t *octant, 
					    ptsrc_t *pointSource, int32_t cycle,
					    int32_t iForce){   
  
    int32_t eindex;       
    int j,iTime, iCoord, iNode;  
    double mu, center_x, center_y, center_z,nodalForceArea;  
    elem_t *elemp;  
    edata_t *edata;  
    tick_t edgeticks;  

    edgeticks = (tick_t)1 << (PIXELLEVEL - octant->level);
      
    /* Calculate the center coordinate of the element */
    center_x = myMesh->ticksize * (octant->lx + edgeticks / 2);
    center_y = myMesh->ticksize * (octant->ly + edgeticks / 2);  
    center_z = myMesh->ticksize * (octant->lz + edgeticks / 2);

    /* Go through my local elements to find which one matches the 
       containing octant. I should have a better solution than this.*/  
    for ( eindex = 0; eindex < myMesh->lenum; eindex++ ) {
    
	int32_t lnid0;    
	lnid0 = myMesh->elemTable[eindex].lnid[0];
    
	if ( ( myMesh->nodeTable[lnid0].x == octant->lx ) &&
	     ( myMesh->nodeTable[lnid0].y == octant->ly ) &&
	     ( myMesh->nodeTable[lnid0].z == octant->lz ) ) {
      
	    /* Sanity check */
	    if ( myMesh->elemTable[eindex].level != octant->level ) {
		fprintf(stderr, "Thread %d: source_init: error\n",myID);
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }            
      
	    /* Fill in the local node ids of the containing element */
	    memcpy ( pointSource->lnid, myMesh->elemTable[eindex].lnid,
		     sizeof(int32_t) * 8 );	
	    break;	
	}
      
    }  /* for all the local elements */
    
    if ( eindex == myMesh->lenum ) {    
	fprintf ( stderr, "Thread %d: source_init: ", myID );     
	fprintf ( stderr, "No element matches the containing octant.\n" );
	MPI_Abort ( MPI_COMM_WORLD, ERROR );
	exit(1);
    }
    
    /* Derive the local coordinate of the source inside the element */  
    pointSource->x =  pointSource -> domainCoords.x[ 0 ] - center_x;
    pointSource->y =  pointSource -> domainCoords.x[ 1 ] - center_y;
    pointSource->z =  pointSource -> domainCoords.x[ 2 ] - center_z;
  
    /* obtain the value of mu to get the moment in the case of extended*/
    /* fault */  
    elemp = &myMesh->elemTable[eindex];  
    edata = (edata_t *)elemp->data;  
    mu = edata->rho * edata->Vs * edata->Vs;    

    if ( theTypeOfSource == POINT ){    
	if ( pointSource->M0 == 0 )
	    pointSource->M0 = pow( 10, 1.5 * theMomentMagnitude + 9.1);    
	myM0+=pointSource->M0;
	pointSource->muArea = pointSource->M0; 
	pointSource->maxSlip = 1; /* Notice that this is a fake slip to
				     obtain the M0 given as data */
    }
    else{
	pointSource->muArea =  mu * pointSource->area;     
	/*    fprintf(stdout,"\n vs=%f mu=%e area=%e\n", edata->Vs,mu, pointSource->area);*/
	fflush(stdout);
	if(cycle == 0)
	    myM0+=fabs(pointSource->muArea*pointSource->maxSlip);  
    }  

    pointSource->edgesize = myMesh->ticksize * edgeticks;
  
    /* Set the nodal forces */  
    source_initnodalforce(pointSource);  /* get the weight for each node */
  
    int isinprocessor=0; /* if 8 will update the search */
    for ( iNode = 0; iNode < 8; iNode++){    
	j = pointSource -> lnid [ iNode ];
         
	if ( myForces [ j ] == NULL  && myForcesCycle[j] == cycle){            
	    myForces [ j ] = 
		( vector3D_t * ) malloc ( sizeof ( vector3D_t ) * 
					  theNumberOfTimeSteps);            
	    myNumberOfNodesLoaded +=1;      
	    myMemoryAllocated +=  sizeof ( vector3D_t ) * theNumberOfTimeSteps;      

	    if ( myForces [ j ] == NULL ) {	
		fprintf(stderr,"Thread %d: load_myForces...: out of mem\n",myID);	
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }

	    for ( iCoord = 0; iCoord < 3; iCoord++)	 
		for ( iTime = 0; iTime < theNumberOfTimeSteps; iTime++) 		  
		    myForces [j][iTime].x[iCoord]=0; 	                

	    isinprocessor+=1;       
	}
    }

    if( isinprocessor == 8 ){
	update_forceinprocessor( iForce, isforceinprocessor, 0 );
    }
  
    for ( iNode = 0; iNode < 8; iNode++){
	j = pointSource -> lnid [ iNode ];
	if ( myForcesCycle[j] == cycle ){  
	    for ( iCoord = 0; iCoord < 3; iCoord++){
		nodalForceArea = (pointSource->nodalForce[iNode][iCoord])*(pointSource->muArea);
		for ( iTime = 0; iTime < theNumberOfTimeSteps; iTime++)	      /* shifted the loops */	  	
		    myForces [j][iTime].x[iCoord]+= 
			(pointSource->displacement[iTime])* nodalForceArea;
	    }
	}
    }      
  
    return 1;
  
}

/**
 * load_myForces_with_point_source: computes the force vector for a 
 *                                  point source and updates myForeces
 * 
 *   input: octant where the source is located
 *          the pointer of the displacement array (time dependant)
 *  output: 0 fails 1 success
 *
 */ 
static int load_myForces_with_point_source_green (octant_t *octant, 
					    ptsrc_t *pointSource, int32_t cycle,
					    int32_t iForce){   
  
    int32_t eindex;       
    int j,iTime, iCoord, iNode;  
    double mu, center_x, center_y, center_z,nodalForce;  
    elem_t *elemp;  
    edata_t *edata;  
    tick_t edgeticks;  

    edgeticks = (tick_t)1 << (PIXELLEVEL - octant->level);
      
    //    fprintf(stdout,"Entro a load_myForces_with_point_source_green");


    /* Calculate the center coordinate of the element */
    center_x = myMesh->ticksize * (octant->lx + edgeticks / 2);
    center_y = myMesh->ticksize * (octant->ly + edgeticks / 2);  
    center_z = myMesh->ticksize * (octant->lz + edgeticks / 2);

    /* Go through my local elements to find which one matches the 
       containing octant. I should have a better solution than this.*/  
    for ( eindex = 0; eindex < myMesh->lenum; eindex++ ) {
    
	int32_t lnid0;    
	lnid0 = myMesh->elemTable[eindex].lnid[0];
    
	if ( ( myMesh->nodeTable[lnid0].x == octant->lx ) &&
	     ( myMesh->nodeTable[lnid0].y == octant->ly ) &&
	     ( myMesh->nodeTable[lnid0].z == octant->lz ) ) {
      
	    /* Sanity check */
	    if ( myMesh->elemTable[eindex].level != octant->level ) {
		fprintf(stderr, "Thread %d: source_init: error\n",myID);
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }            
      
	    /* Fill in the local node ids of the containing element */
	    memcpy ( pointSource->lnid, myMesh->elemTable[eindex].lnid,
		     sizeof(int32_t) * 8 );	
	    break;	
	}
      
    }  /* for all the local elements */
    
    if ( eindex == myMesh->lenum ) {    
	fprintf ( stderr, "Thread %d: source_init: ", myID );     
	fprintf ( stderr, "No element matches the containing octant.\n" );
	MPI_Abort ( MPI_COMM_WORLD, ERROR );
	exit(1);
    }
    
    /* Derive the local coordinate of the source inside the element */  
    pointSource->x =  pointSource -> domainCoords.x[ 0 ] - center_x;
    pointSource->y =  pointSource -> domainCoords.x[ 1 ] - center_y;
    pointSource->z =  pointSource -> domainCoords.x[ 2 ] - center_z;
  
    /* obtain the value of mu to get the moment in the case of extended*/
    /* fault */  
    elemp = &myMesh->elemTable[eindex];  
    edata = (edata_t *)elemp->data;  
    mu = edata->rho * edata->Vs * edata->Vs;    


    pointSource->edgesize = myMesh->ticksize * edgeticks;
  
    /* Set the nodal forces */  
    source_initnodalforce_greensfunction(pointSource);  /* get the weight for each node */
  
    int isinprocessor=0; /* if 8 will update the search */
    for ( iNode = 0; iNode < 8; iNode++){    
	j = pointSource -> lnid [ iNode ];
         
	if ( myForces [ j ] == NULL  && myForcesCycle[j] == cycle){            
	    myForces [ j ] = 
		( vector3D_t * ) malloc ( sizeof ( vector3D_t ) * 
					  theNumberOfTimeSteps);            
	    myNumberOfNodesLoaded +=1;      
	    myMemoryAllocated +=  sizeof ( vector3D_t ) * theNumberOfTimeSteps;      

	    if ( myForces [ j ] == NULL ) {	
		fprintf(stderr,"Thread %d: load_myForces...: out of mem\n",myID);	
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }

	    for ( iCoord = 0; iCoord < 3; iCoord++)	 
		for ( iTime = 0; iTime < theNumberOfTimeSteps; iTime++) 		  
		    myForces [j][iTime].x[iCoord]=0; 	                

	    isinprocessor+=1;       
	}
    }

    if( isinprocessor == 8 ){
	update_forceinprocessor( iForce, isforceinprocessor, 0 );
    }


    //    fprintf(stdout,"LLego");
  
    for ( iNode = 0; iNode < 8; iNode++){
	j = pointSource -> lnid [ iNode ];
	if ( myForcesCycle[j] == cycle ){  
	    for ( iCoord = 0; iCoord < 3; iCoord++){
		nodalForce = (pointSource->nodalForce[iNode][iCoord]);
		for ( iTime = 0; iTime < theNumberOfTimeSteps; iTime++)	{      /* shifted the loops */	  	
		    myForces [j][iTime].x[iCoord]+= 
			(pointSource->displacement[iTime])* nodalForce;
		    //		    		    fprintf(stdout,"\n %lf ",pointSource->displacement[iTime]);
		    
		}
		//	fprintf(stdout,"\n ");
	    }
	}
    }      
  
    return 1;  
}


/*
 * update_point_source: updates a source in the context of an extended 
 *                      fault 
 *   
 *
 */
static void update_point_source (ptsrc_t *pointSource, 
				 int32_t downdipindex, 
				 int32_t alongstrikeindex,
				 double area){
    int iWindow, iTime;
    double maxSlip = 0;

    for(iTime=0; iTime < theNumberOfTimeSteps; iTime++)
	pointSource->displacement[iTime]=0.; 
  
  
    for ( iWindow=0; iWindow<theNumberOfTimeWindows; iWindow++ ){	  

	pointSource->rake = 
	    theRakeMatrix [iWindow][ downdipindex] [ alongstrikeindex ];   

	pointSource->maxSlip = theSlipMatrix[ iWindow ] [ downdipindex] [ alongstrikeindex ];  
	maxSlip +=  pointSource->maxSlip;
	pointSource->area = area;     
	pointSource->delayTime = 
	    compute_initial_time ( pointSource->localCoords,
				   theExtendedHypLocCart,
				   theExtendedAverageRuptureVelocity);   
	pointSource->delayTime += theWindowDelay [ iWindow ];   


	compute_source_function( pointSource );       

    }

    pointSource->maxSlip= maxSlip;
     
}

static void 
compute_point_source_strike_srfh (ptsrc_t* ps, int32_t iSrc)
{

    vector3D_t pivot, pointInNorth, unitVec;

    double norm,fi ;


    if( theLonlatOrCartesian == 1 )return;

  
    if( theLonlatOrCartesian == 0 ){

	pivot = compute_domain_coords_linearinterp(theSourceLonArray[iSrc],
						   theSourceLatArray[iSrc],
						   theSurfaceCornersLong , 
						   theSurfaceCornersLat,
						   theRegionLengthEastM,
						   theRegionLengthNorthM );
    
	pointInNorth = compute_domain_coords_linearinterp(theSourceLonArray[iSrc],
							  theSourceLatArray[iSrc]+.1,
							  theSurfaceCornersLong , 
							  theSurfaceCornersLat,
							  theRegionLengthEastM,
							  theRegionLengthNorthM );

   
	/* Compute Unit Vector */
	unitVec.x[1]=pointInNorth.x[1]-pivot.x[1];
	unitVec.x[0]=pointInNorth.x[0]-pivot.x[0];

	norm = pow( pow( unitVec.x[0], 2)+ pow(unitVec.x[1],2), .5);
   
	unitVec.x[1]= unitVec.x[1]/norm;
	unitVec.x[0]= unitVec.x[0]/norm;
  
	/* Compute the Angle North-X axis */
	fi = atan( unitVec.x[0]/unitVec.x[1]);
    
	if(  unitVec.x[1] < 0 ) /* in rad*/
	    fi = fi + PI;

	/* Compute the strike */ 
	ps->strike =90+ ps->strike-( 180*fi/PI);
   
    
    }  
  
}



/*
 * update_point_source: updates a source in the context of an extended 
 *                      fault 
 *   
 *
 */
static void update_point_source_srfh (ptsrc_t *pointSource, int32_t isource){
    int iTime;

    for(iTime=0; iTime < theNumberOfTimeSteps; iTime++)
	pointSource->displacement[iTime]=0.; 

    pointSource->strike =  theSourceStrikeArray[isource];

    compute_point_source_strike_srfh(pointSource,isource);


    pointSource->dip    =  theSourceDipArray[isource];  
    pointSource->rake   =  theSourceRakeArray[isource];
    pointSource->area   =  theSourceAreaArray[isource];     
    pointSource->maxSlip=  theSourceSlipArray[isource];
    pointSource->nt1    =  theSourceNt1Array[isource];
    pointSource->dtfunction = theSourceDtArray[isource];
    pointSource->delayTime  = theSourceTinitArray[isource];
    pointSource->slipfundiscrete    = theSourceSlipFunArray[isource];
    pointSource->sourceFunctionType = theSourceFunctionType;


    compute_source_function( pointSource );


    return;
}
    

/*
 * update_point_green: updates a source in the context of an extended 
 *                     fault 
 *   
 *
 */
static void update_point_source_green (ptsrc_t *pointSource, int32_t isource){
    int iTime;

    //    fprintf(stdout,"Entro a  update_point_source_green");
    for(iTime=0; iTime < theNumberOfTimeSteps; iTime++)
	pointSource->displacement[iTime]=0.; 


    pointSource->weightX    =  theSourceWXArray[isource];  
    pointSource->weightY    =  theSourceWYArray[isource];
    pointSource->weightZ    =  theSourceWZArray[isource]; 
    
    pointSource->nt1    =  theSourceNt1Array[isource];
    pointSource->dtfunction = theSourceDtArray[isource];
    pointSource->delayTime  = theSourceTinitArray[isource];
    pointSource->sourcefunctiondiscrete = theSourceSourceFunArray[isource];
    pointSource->sourceFunctionType = theSourceFunctionType;

    //    fprintf(stdout,"\nXX %lf %lf %lf %d %lf %lf ", pointSource->weightX,
    //	                                             pointSource->weightY, 
    //	                                             pointSource->weightZ,
    //	    pointSource->nt1,pointSource->dtfunction,pointSource->delayTime);
	    

    //    fprintf(stdout,"\nAntes compute_source_function");
    compute_source_function_green( pointSource );
    //fprintf(stdout,"\nDespues compute_source_function");

    return;
}
    


/* 
 * init_planewithkinks_mapping:
 *
 */
static int init_planewithkinks_mapping(){
   
    int iKink, iDir;
  
    double distance;
    double module=0; 
  
    theTraceLengthAccumulated = 
	(double *) malloc ( sizeof( double ) * theNumberOfKinks );
  
  
    theTraceVectors = 
	( vector3D_t * ) malloc( sizeof(vector3D_t) * (theNumberOfKinks-1) );
  
    if ( theTraceLengthAccumulated == NULL || theTraceVectors == NULL ) { 
    
	fprintf ( stderr, 
		  "Thr. %d: Initializing planewithkinks mapping: out of memory\n",
		  myID);
    
	return -1;
    
    }
  
    /* Compute the accumulated distance along the fault trace */  
    theTraceLengthAccumulated[0] = 0;

    /*------------------------------------------------------------------------*/
  
    /*Transform theKinkLonArray and theKinkLatArray into domain coords X and Y */    
    for ( iKink=0 ; iKink < theNumberOfKinks; iKink++ ){
     
	theKinksDomainCoord[iKink] = 
	    compute_domain_coords_linearinterp(theKinkLonArray[iKink],
					       theKinkLatArray[iKink],
					       theSurfaceCornersLong , 
					       theSurfaceCornersLat,   
					       theRegionLengthEastM,
					       theRegionLengthNorthM );

    }
 

    /*----------------------------------------------------------------------*/
    for ( iKink=1 ; iKink < theNumberOfKinks; iKink++ ){
    
	distance= 
	    pow(theKinksDomainCoord[iKink].x[0]-theKinksDomainCoord[iKink-1].x[0],2)+
	    pow(theKinksDomainCoord[iKink].x[1]-theKinksDomainCoord[iKink-1].x[1],2);

	theTraceLengthAccumulated [ iKink ] = 
	    theTraceLengthAccumulated [ iKink-1 ] +	sqrt( distance ); 
    
    }

    theTotalTraceLength = theTraceLengthAccumulated [ theNumberOfKinks -1];
  
  
    /* Compute the trace vectors given in Global Coordinate */
    for ( iKink=0 ; iKink < theNumberOfKinks-1; iKink++) {        
	theTraceVectors[iKink].x[0] =
	    theKinksDomainCoord[iKink+1].x[0]-
	    theKinksDomainCoord[iKink].x[0];      
	theTraceVectors [ iKink ].x[1] =
	    theKinksDomainCoord[iKink+1].x[1]-
	    theKinksDomainCoord[iKink].x[1];      
	theTraceVectors [ iKink ].x[2] = 0;


	/* make it unitary */
      
	for ( iDir = 0; iDir < 3; iDir++) 
	    module += pow(theTraceVectors [ iKink ].x[iDir],2);
      
	module = sqrt(module);
    
	for ( iDir = 0; iDir < 3; iDir++)
	    theTraceVectors [ iKink ].x[iDir] = 
		theTraceVectors [ iKink ].x[iDir] / module;
          
    }
  

    /* Normalize theTraceLengthAccumulated*/

    for ( iKink=0 ; iKink < theNumberOfKinks; iKink++)
	theTraceLengthAccumulated [ iKink ] = 
	    theTraceLengthAccumulated [ iKink ] / theTotalTraceLength;
  
    return 1;

}




/*
 *
 *  compute_global_coords_mapping: it assumes that the x[0] component is along 
 *                                 strike  and x[1] is downdip 
 *
 */
static vector3D_t compute_global_coords_mapping(vector3D_t point){

    int iKink=-1, iKinkOrigin=-1;
  
    double normalizedDistance;

    vector3D_t origin;
  
    vector3D_t pointMapped;
  
    /* Normalize x[0]  value with the totalLength */
  
    normalizedDistance = point.x[0]/theTotalTraceLength;    
  
    /* Check between which elements is this length */
  
    /* do not make it more complex unless you really want to do something
       more complex in geometry */
  
    while ( iKinkOrigin < 0 ){
    
	iKink=iKink +1;
    
	if ( theTraceLengthAccumulated[iKink] <= normalizedDistance &&
	     theTraceLengthAccumulated[iKink+1] >= normalizedDistance )
	    iKinkOrigin = iKink;
    
    
    }     

  
    /* Using the value of the vector associated to this initial point 
       of the interval where it is located and the vector compute the 
       global coordinates */
  

    /* Transform the long lat to X Y global */

    origin.x[0] = theKinksDomainCoord[ iKinkOrigin].x[0];
  
    origin.x[1] = theKinksDomainCoord[ iKinkOrigin].x[1];
  
    origin.x[2] = point.x[1] + theExtendedDepthM;
  
  
    /* Compute the global coordinates using the unit vector */

  
    pointMapped.x[0] = origin.x[0] + ( point.x[0]-
				       theTraceLengthAccumulated[ iKinkOrigin ] * 
				       theTotalTraceLength ) *
	theTraceVectors[iKinkOrigin].x[0];
  
    pointMapped.x[1] = origin.x[1] + ( point.x[0]-
				       theTraceLengthAccumulated[ iKinkOrigin ] * 
				       theTotalTraceLength ) *
	theTraceVectors[iKinkOrigin].x[1];

    pointMapped.x[2] =  origin.x[2];
  
  
  
    return pointMapped;
  
}


/*
 *  
 *  compute_strike_planewithkinks:
 *
 */
static double compute_strike_planewithkinks(vector3D_t point){

    int iKink=-1, iKinkOrigin=-1;
  
    double normalizedDistance, strike;

    /* Normalize x[0]  value with the totalLength */

    normalizedDistance = point.x[0]/theTotalTraceLength;    

    /* Check between which elements is this length */
  
    /* do not make it more complex unless you really want to do something
       more complex in geometry */
  
    while ( iKinkOrigin < 0 ){

	iKink=iKink +1;

	if ( theTraceLengthAccumulated[iKink] <= normalizedDistance &&
	     theTraceLengthAccumulated[iKink+1] >= normalizedDistance )      

	    iKinkOrigin = iKink;
	
    
    }     
   
    if( theTraceVectors[iKinkOrigin].x[1] >= 0) 

	strike = acos(theTraceVectors[iKinkOrigin].x[0]);

    else if ( theTraceVectors[iKinkOrigin].x[0] < 0 &&
	      theTraceVectors[iKinkOrigin].x[1] <= 0 ) 

	strike = ( 3*PI/2 ) - acos(theTraceVectors[iKinkOrigin].x[0]);

    else if ( theTraceVectors[iKinkOrigin].x[0] > 0 &&
	      theTraceVectors[iKinkOrigin].x[1] < 0 ) 

	strike = acos(theTraceVectors[iKinkOrigin].x[0]) + 3*PI/2;

  
    /* compute deg that is what is take from the code */


    strike = 180* strike / PI;


    return strike;

}


/**
 * tag_myForcesCycle: it assigns a number to each element of
 *                    myForces in myForcesCycle. myForcesCycle
 *                    will be modified to choose in wich cycle
 *                    it should be written.
 *
 *  *   input: octant where the source is located
 *             the pointer of the displacement array (time dependant)
 *  output: -1 fails 1 success
 *
 */ 
static int tag_myForcesCycle(octant_t *octant, ptsrc_t *pointSource){   
  
    int32_t eindex;     
  
    int j, iNode;
     
    tick_t edgeticks;  
    edgeticks = (tick_t)1 << (PIXELLEVEL - octant->level);
     
    /* Go through my local elements to find which one matches the 
       containing octant.*/  
    for ( eindex = 0; eindex < myMesh->lenum; eindex++ ) {    
	int32_t lnid0;    
	lnid0 = myMesh->elemTable[eindex].lnid[0];    
	if ( ( myMesh->nodeTable[lnid0].x == octant->lx ) &&
	     ( myMesh->nodeTable[lnid0].y == octant->ly ) &&
	     ( myMesh->nodeTable[lnid0].z == octant->lz ) ) {      
	    /* Sanity check */
	    if ( myMesh->elemTable[eindex].level != octant->level ) {
		fprintf(stderr, "Thread %d: source_init: internal error\n",
			myID);
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }
      
	    /* Fill in the local node ids of the containing element */
	    memcpy ( pointSource->lnid, myMesh->elemTable[eindex].lnid,
		     sizeof(int32_t) * 8 );      
	    break;
	}

    }  /* for all the local elements */
    
    if ( eindex == myMesh->lenum ) {
	fprintf ( stderr, "Thread %d: source_init: ", myID );    
	fprintf ( stderr, "No element matches the containing octant.\n" );
	return -1;
    }
         
    for ( iNode = 0; iNode < 8; iNode++){    
	j = pointSource -> lnid [ iNode ];           
	if ( myForcesCycle [ j ] == -1 ){            	
	    myForcesCycle [ j ] = 1;      	  
	    myNumberOfNodesLoaded +=1;
	}

    }
  
    return 1;
  
}

void update_forceinprocessor(int32_t iForce, char *inoutprocessor, int onoff){

    int32_t whichByte, whichBit;

    char mask;

    whichByte = iForce/8;
    whichBit = 7 - iForce % 8; 

    mask = ( char )(pow(2,whichBit)*onoff);

    inoutprocessor[whichByte] =  inoutprocessor[whichByte] | mask;

    return;

}


int is_forceinprocessor(int32_t iForce, char *inoutprocessor){

    int32_t whichByte, whichBit;

    char mask,test;

    whichByte = iForce/8;
    whichBit = 7 - iForce % 8; 


    mask = ( char )pow(2,whichBit);

    test = inoutprocessor[whichByte] & mask;

    if ( test == mask )

	return 1;

    else 
    
	return 0;

}


/**
 * Close a file descriptor (FILE*) and set it to NULL.
 */
static int
close_file (FILE** fp)
{
    if (fp == NULL) {
	/* this is a programming error */
	ABORT_PROGRAM ("Passed NULL file pointer to close_file()");
	return -1;
    }

    if (*fp != NULL) {
	return fclose (*fp);
    }

    return 0;
}


/* 
 * open_filesforcomputemyForces: 
 */
static int
open_filescomputemyForces( const char*	physicsin, FILE** fpsourcedescription,
			   FILE** fptemporalsource, char*  temporalsource ){
    char sourceDirectoryOutput[256];
  
    *fpsourcedescription = NULL;
    *fptemporalsource    = NULL;
  
    if (myID == 0) { /* only PE 0 reads physics.in */
	FILE* fpphysics;
	char  sourceDescription[256];
	int   pt_ret;
    
	if ( (fpphysics = fopen (physicsin, "r")) == NULL) {
	    perror ("Could not open physics.in");
	    fprintf (stderr, "Filename for physics.in: %s\n", physicsin);
	    return -1;
	}
    
	pt_ret = parsetext (fpphysics, "source_directory_output", 's',
			    sourceDirectoryOutput);
	fclose (fpphysics);
    
	if (pt_ret != 0) {
	    fprintf (stderr, "Could not get 'source_directory_output' from "
		     "physics.in file: \'%s\'\n", physicsin);
	    return -1;
	}

	sprintf (sourceDescription, "%s/sourcedescription.out" , 
		 sourceDirectoryOutput);    

	*fpsourcedescription = fopen(sourceDescription, "w");
	if (*fpsourcedescription == NULL) {
	    perror ("Could not open source description file");
	    fprintf(stderr, "Filename for source description: '%s'\n",
		    sourceDescription );
	    return -1;
	}
    }
  
    /* all PEs execute this */
    MPI_Bcast (sourceDirectoryOutput, 256, MPI_CHAR, 0, MPI_COMM_WORLD);
    sprintf (temporalsource, "%s/force_process.%d", sourceDirectoryOutput,myID);
    *fptemporalsource = fopen (temporalsource, "w");
  
    if (*fptemporalsource == NULL) {
	perror ("Could not open temporary file for source data");
	fprintf (stderr, "temporary source file: %s\n", temporalsource);
	close_file (fpsourcedescription);  /* cleanup, only affect PE 0 */
	return -1;
    }
  
    return 1;
}


/*
 *  read_domain: reread the domain specifications from physics.in
 *               if prefered pass it 
 *
 *      output :  -1 fail  
 *                 1 success
 *
 */
static int read_domain(FILE *fp){ 
  

    double region_origin_latitude_deg, region_origin_longitude_deg;
    double region_azimuth_leftface_deg;  
    double region_depth_shallow_m, region_length_east_m; 
    double region_length_north_m, region_depth_deep_m;
  
  
    if ((parsetext(fp, "region_origin_latitude_deg", 'd', 
		   &region_origin_latitude_deg) != 0) ||
	(parsetext(fp, "region_origin_longitude_deg", 'd',
		   &region_origin_longitude_deg) != 0) ||
	(parsetext(fp, "region_depth_shallow_m", 'd',
		   &region_depth_shallow_m) != 0) ||
	(parsetext(fp, "region_length_east_m", 'd',
		   &region_length_east_m) != 0) ||
	(parsetext(fp, "region_length_north_m", 'd',
		   &region_length_north_m) != 0) ||
	(parsetext(fp, "region_depth_deep_m", 'd',
		   &region_depth_deep_m) != 0) ||
	(parsetext(fp, "region_azimuth_leftface_deg", 'd',
		   &region_azimuth_leftface_deg) )) {

	fprintf(stderr, "Error domain and source dir from physics.in");
	fflush(stdout);
	return -1;
    }
  

    /*Assign to Global variables*/   
  
    theRegionOriginLatDeg  = region_origin_latitude_deg;  
    theRegionOriginLongDeg = region_origin_longitude_deg;  
    theRegionDepthShallowM = region_depth_shallow_m;  
    theRegionLengthEastM   = region_length_east_m;  
    theRegionLengthNorthM  = region_length_north_m;
    theRegionDepthDeepM    = region_depth_deep_m;  
    theRegionAzimuthDeg    = region_azimuth_leftface_deg;

    return 1;

}



/*
 * read_filter :  if myForce vector is going to be read extra data is required 
 *     
 *       input :  fp - to source.in 
 *             
 *     
 *      output :  -1 fail  
 *                 1 success
 *
 *      Notes  :  a low pass butterworth filter is used
 *
 */
static int read_filter(FILE *fp){

   
    int source_is_filtered, number_of_poles;

    double threshold_frequency;  
  
    if ( parsetext( fp, "source_is_filtered", 'i', &source_is_filtered) != 0){
	fprintf(stderr, "Error parsing source.in reading filter parameters\n");	
	return -1;
    }
    if ( source_is_filtered == 1 ){    

	if ( parsetext( fp, "threshold_frequency", 'd', &threshold_frequency) != 0){
	    fprintf(stderr, "Error parsing source.in reading filter parameters\n");	
	    return -1;
	}
    
	if ( (parsetext(fp, "threshold_frequency", 'd', 
			&threshold_frequency) != 0) ||
	     (parsetext(fp, "number_of_poles", 'i',
			&number_of_poles)     != 0) ) {
	    fprintf(stderr, "Error parsing source.in reading filter parameters\n");	
	    return -1;
	}

    }


    theSourceIsFiltered   = source_is_filtered;
    theThresholdFrequency = threshold_frequency;
    theNumberOfPoles      = number_of_poles;
    
    return 1;

}


/*
 *  read_common_all_formats : read the common varibles to all formats
 *
 *           input : fp  - to source.in
 *                   fpw - to windows.in 
 *
 *          output :  -1 fail  
 *                     1 success
 *
 *           
 *           Notes :   window.in has the delay times of every pulse or shape
 *                     function that we use to represent the slip function.
 *
 *
 */
static int read_common_all_formats ( FILE *fp ){


    double average_risetime_sec, ricker_Ts, ricker_Tp;

    char  source_function_type[64];

    source_function_t sourceFunctionType;
  
    /* 
     *  From source.in slip function parameters 
     */
  
  
    if ( (parsetext(fp, "source_function_type", 's',
		    &source_function_type) != 0) ) {    
	fprintf(stderr, "Error parsing files from source.in");    
	return -1;
      
    }
    
    if ( strcasecmp(source_function_type, "ramp") == 0 )       
	sourceFunctionType = RAMP;      
    else if ( strcasecmp(source_function_type, "sine") == 0 )             
	sourceFunctionType = SINE;            
    else if ( strcasecmp(source_function_type, "quadratic") == 0 )        
	sourceFunctionType = QUADRATIC;        
    else if ( strcasecmp(source_function_type, "ricker") == 0 )             
	sourceFunctionType = RICKER;       
    else if ( strcasecmp(source_function_type, "exponential") == 0 )       
	sourceFunctionType = EXPONENTIAL;  
    else if ( strcasecmp(source_function_type, "discrete") == 0 )       
	sourceFunctionType = DISCRETE;  
      
    else {          
	fprintf(stderr, "Unknown excitation type %s\n", source_function_type);      
	return -1;      
    }
        
    if ( sourceFunctionType == RAMP       ||
	 sourceFunctionType == SINE       ||    
	 sourceFunctionType == QUADRATIC  ||
	 sourceFunctionType == EXPONENTIAL ) {      
    
	if ( (parsetext(fp, "average_risetime_sec", 'd', 
			&average_risetime_sec) != 0)){	
	    fprintf(stderr, "Cannot get time function parameters\n");
	    return -1;
	}
    }
    
    if (  sourceFunctionType == RICKER ) {
	if ( (parsetext(fp, "ricker_Ts", 'd', &ricker_Ts) != 0)||
	     (parsetext(fp, "ricker_Tp", 'd', &ricker_Tp) != 0)) {	
	    fprintf(stderr, "Cannot get parameters to ricker's function\n");
	    return -1;
	}
    }
         
    theSourceFunctionType = sourceFunctionType;
    theAverageRisetimeSec = average_risetime_sec;
    theRickerTs = ricker_Ts;
    theRickerTp = ricker_Tp;;

    return 1;

}

/*
 *  read_point_source :
 *
 *                     input : fp  - to source.in                    
 *
 *                    output :  -1 fail  
 *                               1 success
 *
 */

static int read_point_source(FILE *fp){

    double hypocenter_lat_deg, hypocenter_long_deg, hypocenter_depth_m;
    double source_strike_deg, source_dip_deg, source_rake_deg;  
    double  moment_magnitude, moment_amplitude;
    double auxiliar[8];
    int iCorner;

    /* You can either give M0 or Mw for a point Source */  
    moment_magnitude = 0;
    moment_amplitude = 0;
    if ( (parsetext(fp, "moment_magnitude", 'd', &moment_magnitude) != 0) )  
	if ( ( parsetext(fp,"moment_amplitude",'d',&moment_amplitude ) != 0) ){
	    fprintf(stderr,"Error moment:read_point_source\n");
	    return -1;
	}
                 
    if ( (parsetext(fp,"lonlat_or_cartesian", 'i',&theLonlatOrCartesian) != 0) ){
	fprintf(stderr,
		"Err lonlat_or_cartesian in source.in point source parameters\n");
	return -1;
    }

    if( theLonlatOrCartesian == 0 ){
        
	if ( (parsetext(fp,"hypocenter_lat_deg", 'd',&hypocenter_lat_deg) != 0) ||
	     (parsetext(fp,"hypocenter_long_deg",'d',&hypocenter_long_deg)!= 0)){
	    fprintf(stderr, "Err in hypocenter lon or lat:read_point_source\n");
	    return -1;
	}    

	parsedarray( fp, "domain_surface_corners", 8 ,auxiliar);
	for ( iCorner = 0; iCorner < 4; iCorner++){    
	    theSurfaceCornersLong[ iCorner ] = auxiliar [ iCorner * 2 ];
	    theSurfaceCornersLat [ iCorner ] = auxiliar [ iCorner * 2 +1 ];      
	}
    }

    /* Here I am not defining another variable for the cartesian it will be 
       controled with lonlatorcartesian */ 
    if( theLonlatOrCartesian == 1 ) 
	if ( (parsetext(fp,"hypocenter_x",'d',&hypocenter_lat_deg ) != 0) ||
	     (parsetext(fp,"hypocenter_y",'d',&hypocenter_long_deg)!= 0)  ){
	    fprintf(stderr, "Err hypocenter x or y:read_point_source\n");
	    return -1;
	}
  
    if ( (parsetext(fp,"hypocenter_depth_m", 'd',&hypocenter_depth_m) != 0) ||
	 (parsetext(fp,"source_strike_deg",  'd',&source_strike_deg)  != 0) ||
	 (parsetext(fp, "source_dip_deg",    'd',&source_dip_deg)     != 0) ||
	 (parsetext(fp, "source_rake_deg",   'd',&source_rake_deg)    != 0)) {

	fprintf(stderr, "Err fields:read_point_source\n");
	return -1;
    }
    
    theMomentMagnitude   = moment_magnitude;  
    theMomentAmplitude   = moment_amplitude;  
    theHypocenterLatDeg  = hypocenter_lat_deg;  
    theHypocenterLongDeg = hypocenter_long_deg;  
    theHypocenterDepthM  = hypocenter_depth_m;  
    theSourceStrikeDeg   = source_strike_deg;  
    theSourceDipDeg      = source_dip_deg;  
    theSourceRakeDeg     = source_rake_deg;
            
    return 1;

}


/*
 *  read_plane_source :
 *
 *           input :    fp  -  to source.in                    
 *                  fpslip  -  to slip table file
 *                  fprake  -  to rake table file
 *
 *
 *          output :  -1 fail  
 *                     1 success
 *
 */

static int  read_plane_source(FILE *fp, FILE *fpslip, FILE *fprake){
  
    int iWindow;      
    int cells_down_dip,cells_along_strike;
    int isminimumedgeautomatic, number_of_time_windows;
    int iCorner;
  
    double cell_size_along_strike_m, cell_size_down_dip_m;
    double lat_deg, long_deg, depth_m;
    double hypocenter_along_strike_m, hypocenter_down_dip_m;
    double average_rupture_velocity, strike_deg,dip_deg, minimum_edge_m;  

    double *auxiliar;

    if ( (parsetext(fp,"number_of_time_windows",'i',&number_of_time_windows) != 0) ) {    
	fprintf( stderr, "Error parsing fields from source.in num of windows\n");
	return -1;
    }
  
    theNumberOfTimeWindows = number_of_time_windows; 
    theWindowDelay= dvector(0,theNumberOfTimeWindows);
  
    parsedarray( fp, "time_windows", theNumberOfTimeWindows,theWindowDelay);
  
    /* Mw and M0  will be calculated */        
    if ((parsetext(fp,"extended_cell_size_down_dip_m",  'd', &cell_size_down_dip_m)   != 0) ||
	(parsetext(fp,"extended_isminimumedgeautomatic",'i', &isminimumedgeautomatic) != 0) ||
	(parsetext(fp,"extended_depth_m"      ,         'd', &depth_m)                != 0) ||
	(parsetext(fp,"extended_cells_along_strike",    'i', &cells_along_strike)     != 0) ||
	(parsetext(fp,"extended_cells_down_dip",        'i', &cells_down_dip)         != 0) ||
	(parsetext(fp, "extended_hypocenter_along_strike_m",
		   'd',&hypocenter_along_strike_m) != 0) ||
	(parsetext(fp, "extended_hypocenter_down_dip_m",'d',  &hypocenter_down_dip_m) != 0) ||
	(parsetext(fp, "extended_average_rupture_velocity", 'd', 	 
		   &average_rupture_velocity) != 0) ){      
	fprintf(stderr, "Cannot query plane source\n");
	return -1;      
    }

    if( theTypeOfSource == PLANE )
	if ((parsetext(fp,"extended_cell_size_along_strike_m",'d',&cell_size_along_strike_m) != 0)||
	    (parsetext(fp,"extended_lat_deg",               'd', &lat_deg)              != 0) ||
	    (parsetext(fp,"extended_long_deg",              'd', &long_deg)             != 0) ||
	    (parsetext(fp, "extended_strike_deg",           'd',  &strike_deg)          != 0) ||
	    (parsetext(fp, "extended_dip_deg"   ,           'd',  &dip_deg)             != 0)){
	    fprintf(stderr, "Cannot query plane source\n");
	    return -1;      
	}
    
    if ( isminimumedgeautomatic == 0 ){
	if ( parsetext(fp,"extended_minimum_edge_m",'d',&minimum_edge_m) != 0 ) return -1;    
    }  
    else    
	minimum_edge_m = theMinimumEdge/1; /* Denominator should be >= 1 */
  
    /* Read slip and rake matrices */
    theExtendedColumns = cells_along_strike;  
    theExtendedRows    = cells_down_dip;
    theSlipMatrix = malloc( sizeof(double **) * theNumberOfTimeWindows);  
    theRakeMatrix = malloc( sizeof(double **) * theNumberOfTimeWindows);

    for ( iWindow=0; iWindow < theNumberOfTimeWindows; iWindow++){    
	theSlipMatrix[iWindow] = dmatrix(0,theExtendedRows-1,0,
					 theExtendedColumns-1);
	theRakeMatrix[iWindow] = dmatrix(0,theExtendedRows-1,0,
					 theExtendedColumns-1);
	read_double_matrix_from_file(fpslip,theSlipMatrix[iWindow],
				     theExtendedRows,theExtendedColumns);
	read_double_matrix_from_file(fprake,theRakeMatrix[iWindow],
				     theExtendedRows,theExtendedColumns);      
    }



    /* corners of the surface */
    auxiliar = (double *)malloc(sizeof(double)*8);
    if ( auxiliar == NULL ) {   
	perror("Failed to allocate space in read planes");
	MPI_Abort ( MPI_COMM_WORLD, ERROR );
	return -1;
    }
  
    parsedarray( fp, "domain_surface_corners", 8 ,auxiliar);
    for ( iCorner = 0; iCorner < 4; iCorner++){    
	theSurfaceCornersLong[ iCorner ] = auxiliar [ iCorner * 2 ];
	theSurfaceCornersLat [ iCorner ] = auxiliar [ iCorner * 2 +1 ];      
    }
  
    free(auxiliar);



    fflush(stdout);
     
    if ( theTypeOfSource == PLANEWITHKINKS )
	if(read_planewithkinks(fp)==-1)return -1;

    theMomentMagnitude = 0;  
    theMomentAmplitude = 0; 
  
    theExtendedCellSizeAlongStrikeM   = cell_size_along_strike_m;  
    theExtendedCellSizeDownDipM       = cell_size_down_dip_m;  
    theExtendedMinimumEdge            = minimum_edge_m;   
    theExtendedLatDeg                 = lat_deg;  
    theExtendedLongDeg                = long_deg;  
    theExtendedDepthM                 = depth_m;  
    theExtendedHypocenterAlongStrikeM = hypocenter_along_strike_m;  
    theExtendedHypocenterDownDipM     = hypocenter_down_dip_m;  
    theExtendedAverageRuptureVelocity = average_rupture_velocity;  
    theExtendedStrikeDeg              = strike_deg;  
    theExtendedDipDeg                 = dip_deg;  
 

    fprintf(stdout,"\n Lat=%f Lon=%f Depth=%f",
	    theExtendedLatDeg,theExtendedLongDeg,theExtendedDepthM)  ;

    return 1;
  
}



/**
 *  read_planewithkinks :
 *
 *           input :    fp  -  to source.in                    
 *                  fpslip  -  to slip table file
 *                  fprake  -  to rake table file
 *           output :  1 - ok -1 fail 
 *
 */
static int
read_planewithkinks (FILE *fp)
{

    int iKink, iCorner;
   
    if ( (parsetext(fp, "extended_number_of_kinks", 'i', &theNumberOfKinks) != 0) ){  
	fprintf(stderr, "Cannot query plane source\n");    
	return -1;
    }
  
    /* read Fault Trace description */
    theKinkLonArray = (double *) malloc( sizeof(double) * theNumberOfKinks);
    theKinkLatArray = (double *) malloc( sizeof(double) * theNumberOfKinks); 		  
    if ( theKinkLonArray == NULL || theKinkLatArray == NULL ) {   
	perror("Failed to allocate space for the fault trace description");
	return -1;
    }
  
    double *auxiliar;
  
    /* trace of the "strike-slip" fault */
    auxiliar = (double *)malloc(sizeof(double)*theNumberOfKinks*2);
    if ( auxiliar == NULL ) {   
	perror("Failed to allocate space in read planes");
	return -1;
    }
  
    parsedarray( fp, "extended_kinks", theNumberOfKinks*2,auxiliar);
    for ( iKink = 0; iKink < theNumberOfKinks; iKink++){    
	theKinkLonArray[ iKink ] = auxiliar [ iKink * 2 ];
	theKinkLatArray[ iKink ] = auxiliar [ iKink * 2 +1 ];      
    }
  
    free(auxiliar);
  
    /* corners of the surface */
    /*  auxiliar = (double *)malloc(sizeof(double)*8); */
    /*   if ( auxiliar == NULL ) {    */
    /*     perror("Failed to allocate space in read planes"); */
    /*     MPI_Abort ( MPI_COMM_WORLD, ERROR ); */
    /*     return -1; */
    /*   } */
  
    /*   parsedarray( fp, "domain_surface_corners", 8 ,auxiliar); */
    /*   for ( iCorner = 0; iCorner < 4; iCorner++){     */
    /*     theSurfaceCornersLong[ iCorner ] = auxiliar [ iCorner * 2 ]; */
    /*     theSurfaceCornersLat [ iCorner ] = auxiliar [ iCorner * 2 +1 ];       */
    /*   } */
  
    /*   free(auxiliar); */
  
    return 1;
  
}





/**
 *  read_srfh_source: read all the parameters fot the standard 
 *                    rupture format hercules version.A matlab
 *                    file is provided (converter.m) to do the
 *                    conversion  from the scec srf to srfh.
 *
 *           input :      fp  -  to source.in                    
 *                  fpcoords  -  to coordinates of the fault
 *                  fpstrike  -   
 *                     fpdip  -
 *                    fprake  -
 *                    fpslip  -
 *                 fpslipfun  -  slip function discretized
 *
 *           output :  1 - ok -1 fail 
 *
 */
static int
read_srfh_source (FILE *fp, FILE *fpcoords, FILE *fparea, FILE *fpstrike,
		  FILE *fpdip, FILE *fprake, FILE *fpslip, FILE *fpslipfun){
    int32_t iSrc;
    double *auxiliar;
    int iCorner, iTime;

    if ( (parsetext(fp, "number_of_point_sources", 'i', &theNumberOfPointSources) != 0) ){  
	fprintf(stderr, "Cannot query number of point source\n");    
	return -1;
    }
  
    theSourceLonArray     = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceLatArray     = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceDepthArray   = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceAreaArray    = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceStrikeArray  = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceDipArray     = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceRakeArray    = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceSlipArray    = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceNt1Array     = malloc( sizeof( int )   * theNumberOfPointSources );
    theSourceTinitArray   = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceDtArray      = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceSlipFunArray = malloc( sizeof( double* ) * theNumberOfPointSources );
  
    if ( (theSourceLonArray    == NULL) || (theSourceLatArray     == NULL) || 
	 (theSourceDepthArray  == NULL) || (theSourceStrikeArray  == NULL) ||
	 (theSourceStrikeArray == NULL) || (theSourceDipArray     == NULL) ||
	 (theSourceRakeArray   == NULL) || (theSourceSlipArray    == NULL) ||
	 (theSourceNt1Array    == NULL) || (theSourceTinitArray   == NULL) ||
	 ( theSourceDtArray    == NULL) || (theSourceSlipFunArray == NULL) ){
	ABORT_PROGRAM ("Error source srfh matrices: read_srfh_source");
    }

    /* How to read the source coordinates (Note 2) */

    /* read fault points description */  
    for ( iSrc = 0; iSrc < theNumberOfPointSources; iSrc++ ){
	fscanf(fpcoords," %lf %lf %lf ", &(theSourceLatArray[iSrc]),
	       &(theSourceLonArray[iSrc]),
	       &(theSourceDepthArray[iSrc])); 
	fscanf(fparea,   " %lf ", &(theSourceAreaArray[iSrc]));
	fscanf(fpstrike, " %lf ", &(theSourceStrikeArray[iSrc]));
	fscanf(fpdip,    " %lf ", &(theSourceDipArray[iSrc]));
	fscanf(fprake,   " %lf ", &(theSourceRakeArray[iSrc]));
	fscanf(fpslip,   " %lf ", &(theSourceSlipArray[iSrc]));
	fscanf(fpslipfun," %d ",  &(theSourceNt1Array[iSrc]));
	
	if(theSourceNt1Array[iSrc]==0){
	  fprintf(stdout," Source %d has %d samples",iSrc,theSourceNt1Array[iSrc]);
	  exit(1);
	}


	fscanf(fpslipfun," %lf ", &(theSourceTinitArray[iSrc]));
	fscanf(fpslipfun," %lf ", &(theSourceDtArray[iSrc]));
    
	theSourceSlipFunArray[iSrc]=malloc(sizeof(double)*theSourceNt1Array[iSrc]);

	for ( iTime = 0; iTime < theSourceNt1Array[iSrc]; iTime++)
	    fscanf(fpslipfun," %lf ", &(theSourceSlipFunArray[iSrc][iTime]));      
    }
 
    /* corners of the surface */
    auxiliar = (double *)malloc(sizeof(double)*8);
    if ( auxiliar == NULL ) {   
	perror(" Alloc auxiliar: read_srfh_source");
	MPI_Abort ( MPI_COMM_WORLD, ERROR );
	return -1;
    }  
    parsedarray( fp, "domain_surface_corners", 8 ,auxiliar);
    for ( iCorner = 0; iCorner < 4; iCorner++){    
	theSurfaceCornersLong[ iCorner ] = auxiliar [ iCorner * 2 ];
	theSurfaceCornersLat [ iCorner ] = auxiliar [ iCorner * 2 +1 ];      
    }
    
    free(auxiliar);

    return 1;
  
}

/**
 *  read_green_source: 
 *
 *           output :  1 - ok -1 fail 
 *
 */
static int
read_green_source (FILE *fp, FILE *fpcoords, FILE *fpweights, FILE *fpsourcefun){

    int32_t iSrc;
    double *auxiliar;
    int iCorner, iTime;

    if ( (parsetext(fp, "number_of_point_sources", 'i', &theNumberOfPointSources) != 0) ){  
	fprintf(stderr, "Cannot query number of point source\n");    
	return -1;
    }
  
    theSourceLonArray       = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceLatArray       = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceDepthArray     = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceWXArray        = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceWYArray        = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceWZArray        = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceNt1Array       = malloc( sizeof( int )   * theNumberOfPointSources );
    theSourceTinitArray     = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceDtArray        = malloc( sizeof( double) * theNumberOfPointSources );
    theSourceSourceFunArray = malloc( sizeof( double* ) * theNumberOfPointSources );
  
    if ( (theSourceLonArray    == NULL) || (theSourceLatArray       == NULL) || 
	 (theSourceDepthArray  == NULL) || (theSourceWXArray        == NULL) ||
	 (theSourceWYArray     == NULL) || (theSourceWZArray        == NULL) ||
	 (theSourceNt1Array    == NULL) || (theSourceTinitArray     == NULL) ||
	 ( theSourceDtArray    == NULL) || (theSourceSourceFunArray == NULL) ){
	ABORT_PROGRAM ("Error source green matrices: read_green_source");
    }

    /* read fault points description */  
    for ( iSrc = 0; iSrc < theNumberOfPointSources; iSrc++ ){
	fscanf(fpcoords," %lf %lf %lf ", &(theSourceLatArray[iSrc]),
	       &(theSourceLonArray[iSrc]),
	       &(theSourceDepthArray[iSrc]));


	fscanf(fpweights,  " %lf ", &(theSourceWXArray [iSrc]));
	fscanf(fpweights,  " %lf ", &(theSourceWYArray [iSrc]));
	fscanf(fpweights,  " %lf ", &(theSourceWZArray [iSrc]));


	//	fprintf(stdout,"\n %lf %lf %lf",theSourceWXArray [iSrc],
	//	theSourceWYArray [iSrc],
	//	theSourceWZArray [iSrc]);

	fscanf(fpsourcefun," %d ",  &(theSourceNt1Array[iSrc]));
	fscanf(fpsourcefun," %lf ", &(theSourceTinitArray[iSrc]));
	fscanf(fpsourcefun," %lf ", &(theSourceDtArray[iSrc]));
      
	//fprintf(stdout,"\n %d %lf %lf",theSourceNt1Array[iSrc],
	//	theSourceTinitArray[iSrc],theSourceDtArray[iSrc]);


 
	theSourceSourceFunArray[iSrc]=malloc(sizeof(double)*theSourceNt1Array[iSrc]);
      

	for ( iTime = 0; iTime < theSourceNt1Array[iSrc]; iTime++)
	    fscanf(fpsourcefun," %lf ", &(theSourceSourceFunArray[iSrc][iTime]));      
    }
  
    /* corners of the surface */
    auxiliar = (double *)malloc(sizeof(double)*8);
    if ( auxiliar == NULL ) {   
	perror(" Alloc auxiliar: read_green_source");
	MPI_Abort ( MPI_COMM_WORLD, ERROR );
	return -1;
    }  
    parsedarray( fp, "domain_surface_corners", 8 ,auxiliar);
    for ( iCorner = 0; iCorner < 4; iCorner++){    
	theSurfaceCornersLong[ iCorner ] = auxiliar [ iCorner * 2 ];
	theSurfaceCornersLat [ iCorner ] = auxiliar [ iCorner * 2 +1 ];      
    }
  
    free(auxiliar);

    fprintf(stdout,"\n exit read green");
  
    return 1;
}


/**
 * search_nodes_octant: computes the force vector for a 
 *                                  point source and updates myForeces
 * 
 */ 
static int search_nodes_octant (octant_t *octant, 
				ptsrc_t *pointSource){   
  
    int32_t eindex;     
  
    tick_t edgeticks;  
    edgeticks = (tick_t)1 << (PIXELLEVEL - octant->level);
  
    
    /* Go through my local elements to find which one matches the 
       containing octant. I should have a better solution than this.*/
  
    for ( eindex = 0; eindex < myMesh->lenum; eindex++ ) {    
	int32_t lnid0;    
	lnid0 = myMesh->elemTable[eindex].lnid[0];
    
	if ( ( myMesh->nodeTable[lnid0].x == octant->lx ) &&
	     ( myMesh->nodeTable[lnid0].y == octant->ly ) &&
	     ( myMesh->nodeTable[lnid0].z == octant->lz ) ) {
      
	    /* Sanity check */
	    if ( myMesh->elemTable[eindex].level != octant->level ) {
		fprintf(stderr, "Thread %d: source_init: internal error\n",
			myID);
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }            
      
	    /* Fill in the local node ids of the containing element */
	    memcpy ( pointSource->lnid, myMesh->elemTable[eindex].lnid,
		     sizeof(int32_t) * 8 );      
	    break;

	}

    }  /* for all the local elements */
  
  
    if ( eindex == myMesh->lenum ) {
	fprintf ( stderr, "Thread %d: source_init: ", myID );    
	fprintf ( stderr, "No element matches the containing octant.\n" );
	MPI_Abort ( MPI_COMM_WORLD, ERROR );
	exit(1);
    }   
    
    return 1;
  
}

static void print_header_myForces(FILE *fptmpsrc){

    int32_t iNode;

    fwrite( &myNumberOfNodesLoaded, sizeof(int32_t),1,fptmpsrc);
    //fprintf(fptmpsrc,"\n%d \n",myNumberOfNodesLoaded);

    for ( iNode = 0; iNode <  myMesh->nharbored; iNode++){
	if( myForcesCycle[iNode] !=-1 ){
	    fwrite(&iNode, sizeof(int32_t),1,fptmpsrc);
	    //fprintf(fptmpsrc,"\n%d \n",iNode);
	}
    }
  
  
}

/* 
 * 
 *  concat_filespercylce_in_fptmpsrc:
 *  
 *
 */ 

static void concat_filespercycle_in_fptmpsrc(char *filenameroot, FILE *fptmpsrc,
					     int numberofcycles){  
    FILE **fp;
    char forcefile[256];
  
    int i,iDir, iCycle,iCycleNode, iDisk, diskAccess;

    int64_t totalMemory;
    int32_t *nodesPerCycle, nodesCheckSum=0;
    int32_t timeSteps, *cycleTimeSteps;  
    int32_t iNode, iTime;
    int32_t *nodesAcumulated;
  
    vector3D_t **buffer;
 
    fp =  malloc(sizeof(FILE *) * numberofcycles);

    nodesPerCycle = malloc( sizeof(int32_t) * numberofcycles);
    nodesAcumulated = malloc( sizeof(int32_t) * numberofcycles);
  

    for (iCycle = 0; iCycle < numberofcycles; iCycle++){
    
	sprintf(forcefile, "%s.%d", filenameroot, iCycle);    
	fp[iCycle]=fopen(forcefile,"r");
    
	if(fp[iCycle] == NULL){      
	    fprintf(stdout,"Err opening concat files");      
	    MPI_Abort(MPI_COMM_WORLD, ERROR);       
	    exit(1);    
	}       
    
	if(fread(&(nodesPerCycle[iCycle]),sizeof(int32_t) , 1, fp[iCycle])==0){
	    fprintf(stdout,"\nError reading files per cycle in concat function");
	    fprintf(stdout,"\nMyID = %d",myID);
	    fprintf(stdout,"\nnumber of cycles = %d Cycle = %d\n",numberofcycles, iCycle);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}
            
	/*fscanf(fp[iCycle],"%d",&(nodesPerCycle[iCycle]));    */
	nodesCheckSum += nodesPerCycle[iCycle];
  
    }

    nodesAcumulated[0]=0;
  
    for (iCycle = 1; iCycle < numberofcycles; iCycle++)  
	nodesAcumulated[iCycle]=nodesPerCycle[iCycle-1]+nodesAcumulated[iCycle-1];
    
    if(nodesCheckSum != myNumberOfNodesLoaded){  
	fprintf(stdout,"Err: the number of Nodes loaded does not match");    
	MPI_Abort(MPI_COMM_WORLD, ERROR);    
	exit(1);
    }
    
    totalMemory= myNumberOfNodesLoaded*theNumberOfTimeSteps*sizeof(vector3D_t);

    fprintf(stdout,"\n %d In concat totalMemoryReq= %ld  BUFFERSIZE=%d",
	    myID,totalMemory,BUFFERSIZE );

  
    if(totalMemory < BUFFERSIZE){
	diskAccess =1;
	timeSteps = theNumberOfTimeSteps;
	buffer = malloc( sizeof( vector3D_t * ) * timeSteps);
	if(buffer == NULL){
	    fprintf(stdout,"Err: totalMemory < BUFFERSIZE");
	    fprintf(stdout,"\nErr: cannot alloc buffer for concat timeSteps=%d",i, timeSteps);
	    fprintf(stdout,"\nErr: %dMem Available  = %dMB", myID,  compute_available_mem());
	    fprintf(stdout,"\nErr: %ddiskAccess=%d timeSteps=%d ",myID,diskAccess,timeSteps);
	    fprintf(stdout,"\nErr: %dtotalMemory=%ld ",myID,totalMemory);
	    fprintf(stdout,"Err: cannot alloc buff for concat");
	    fflush(stdout);
	    MPI_Abort(MPI_COMM_WORLD, ERROR);
	    exit(1);
	}
	for( i = 0; i < timeSteps; i++){
	    buffer[i] = malloc( sizeof( vector3D_t ) * myNumberOfNodesLoaded );
	    if(buffer[i] == NULL){
		fprintf(stdout,"Err: cannot alloc buffer[i] for concat, myNumberOfNodesLoaded=%d",
			myNumberOfNodesLoaded );
		fprintf(stdout,"Err: totalMemory= %d < BUFFERSIZE=%d",totalMemory,BUFFERSIZE );
		fprintf(stdout,"\nErr: %dMem Available  = %dMB", myID,  compute_available_mem());
		fprintf(stdout,"\nErr: %ddiskAccess=%d timeSteps=%d ",myID,diskAccess,timeSteps);
		fprintf(stdout,"\nErr: %dtotalMemory=%d ",myID,totalMemory);	
		fflush(stdout);
		MPI_Abort(MPI_COMM_WORLD, ERROR);
	
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }
	}
    
	cycleTimeSteps = malloc( sizeof(int32_t)*diskAccess );
	cycleTimeSteps[0]=timeSteps;    
    }
  
    else
	{              
	    diskAccess =  (int)floor(totalMemory/BUFFERSIZE)+1;
	    timeSteps = (int)floor(theNumberOfTimeSteps/diskAccess);
      
	    buffer = malloc( sizeof( vector3D_t * ) * (timeSteps));
	    if(buffer == NULL){
		fprintf(stdout,"Err: cannot alloc buff for concat,timeSteps %d ", timeSteps);
		fprintf(stdout,"Err: diskAccess=%d timeSteps=%d",diskAccess,timeSteps);

		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	    }
	    for( i = 0; i < timeSteps; i++){
		buffer[i] = malloc( sizeof( vector3D_t ) * myNumberOfNodesLoaded );
		if( buffer[i] == NULL ){	  
		    fprintf(stdout,"Err: cannot alloc buffer[i] for concat, myNumberOfNodesLoaded=%d",
			    myNumberOfNodesLoaded );
		    fprintf(stdout,"Err: totalMemory= %ld > BUFFERSIZE=%d",totalMemory,BUFFERSIZE );
		    fprintf(stdout,"\nErr: %dMem Available  = %dMB", myID,  compute_available_mem());
		    fprintf(stdout,"\nErr: %ddiskAccess=%d timeSteps=%d ",myID,diskAccess,timeSteps);
		    fprintf(stdout,"\nErr: %dtotalMemory=%d ",myID,totalMemory);
		    fflush(stdout);
		    MPI_Abort(MPI_COMM_WORLD, ERROR);
		    exit(1);
		}
	    }
    
	    cycleTimeSteps = malloc( sizeof(int32_t ) * ( diskAccess + 1 ));
      
	    for ( i = 0; i < diskAccess; i++)
		cycleTimeSteps[i]=timeSteps;
      
	    if ( theNumberOfTimeSteps % diskAccess != 0 ) {
		cycleTimeSteps [ diskAccess ] =theNumberOfTimeSteps % diskAccess;
		diskAccess += 1;
	    }
	}
  
  
    for (iDisk = 0; iDisk < diskAccess; iDisk++){
	for (iCycle = 0; iCycle < numberofcycles; iCycle++){
	    for ( iTime = 0; iTime < cycleTimeSteps[iDisk]; iTime++ ){		
		for ( iCycleNode = 0; iCycleNode < nodesPerCycle[iCycle] ; iCycleNode++){
		    iNode=nodesAcumulated[iCycle]+iCycleNode;	  
		    for (iDir=0; iDir<3; iDir++)

			if(fread(&(buffer[iTime][iNode].x[iDir]),sizeof(double) , 1, fp[iCycle])==0){
			    fprintf(stdout,"Error reading files per cycle in concat function");
			    fprintf(stdout,"\nMyID = %d",myID);
			    fprintf(stdout,"\nnumber of cycles = %d Cycle = %d Time= %d \n",numberofcycles, iCycle, iTime);
			    MPI_Abort(MPI_COMM_WORLD, ERROR);
			    exit(1);
			}	
		    /*fscanf(fp[iCycle],"%lf",&( buffer[iTime][iNode].x[iDir]));*/
		}
	    }
	}

	for ( iTime = 0; iTime < cycleTimeSteps[iDisk]; iTime++ ){
	    for (iCycle = 0; iCycle < numberofcycles; iCycle++){	
		for ( iCycleNode = 0; iCycleNode < nodesPerCycle[iCycle] ; iCycleNode++ ){
		    iNode=nodesAcumulated[iCycle]+iCycleNode;
		    for (iDir=0; iDir<3; iDir++){
			fwrite( &(buffer[iTime][iNode].x[iDir]), sizeof(double), 1,fptmpsrc);
			//fprintf(fptmpsrc," %lf ",buffer[iTime][iNode].x[iDir]);
		    }
		}	

	    }  
	}     
    }

    for( i = 0; i < timeSteps; i++)
	free(buffer[i]);
  
    free(buffer);
    free(nodesPerCycle);
  
    free(nodesAcumulated);
  
  
    fflush(stdout);


    return;
  
}

/**
 * fill_myForces_cycle
 *
 */
static void fill_myForces_cycle( int cycle) {
    
    int32_t lnid, i, j;        
  
    for ( j = numStepsNecessary; j < theNumberOfTimeSteps; j++ ) {
	for ( lnid = 0; lnid <  myMesh->nharbored; lnid++){ 
	    if ( myForcesCycle [ lnid ] ==  cycle ){ 
		for ( i = 0; i < 3;  i++ ){			
		    myForces[lnid][j].x[i]= myForces[lnid][numStepsNecessary-1].x[i];
		}
	    }
	}
    }
  
    return;
} 









/* 
 * compute_myForces_planes: compute the force vector for each 
 *                          node in myMesh
 *
 *
 */
static int  compute_myForces_planes(const char *physicsin){
  
    FILE *fpdsrc, *fptmpsrc;
    
    double *grdStrk, *grdDip;
    double minEdge, minEdgeAlongStrike, minEdgeDownDip, area;
  
    int iWindow,i,j,currentDownDip,currentAlongStrike,iTri;    
    int strkCells, dpCells;
    int cellPntsStrk, cellPntsDip, ndsStrk, ndsDip, iCorner,jCorner;

    int quadortri=1;

    int32_t iForce=0;
    octant_t *octant;
    ptsrc_t pntSrc;          
    vector3D_t origin, p [ 4 ],localCoords,globalCoords;

    /* i/o related vars */
    int32_t nodesPerCycle, memoryRequiredPerNode, memoryPerCycle;
    int32_t iCycle, iInternal, iNode; 

    /* convinient variables */
    double lat   = theExtendedLatDeg;
    double lon   = theExtendedLongDeg;
    double depth = theExtendedDepthM;

    char fileNameRoot[256];
      
    if(open_filescomputemyForces(physicsin, &fpdsrc, &fptmpsrc,fileNameRoot) == -1)
	return -1;

    if( theTypeOfSource == PLANE ){
	/*    origin=compute_cartesian_coords(lat,lon,depth );   */
	origin=compute_domain_coords_linearinterp(lon, lat,
						  theSurfaceCornersLong , 
						  theSurfaceCornersLat,
						  theRegionLengthEastM,
						  theRegionLengthNorthM );         
	origin.x[2]=  theExtendedDepthM;
	/*pntSrc.domainCoords.x[2]= depth;*/

	/*     fprintf(stdout,"\n ID = %d ", myID,origin.x[0],origin.x[1],origin.x[2]);*/
	/*fprintf(stdout,"\n ID = %d origin=(%lf,%lf,%lf)", myID,origin.x[0],origin.x[1],origin.x[2]);
	  fflush(stdout);*/
     
    }

    MPI_Barrier( MPI_COMM_WORLD );
  

    theExtendedHypLocCart.x[0]=theExtendedHypocenterAlongStrikeM;    
    theExtendedHypLocCart.x[1]=theExtendedHypocenterDownDipM;
    theExtendedHypLocCart.x[2]=0;
 
    strkCells = theExtendedColumns;
    dpCells   = theExtendedRows;            
  
    /* Variables that all poinSources will share */
    if ( theTypeOfSource == PLANEWITHKINKS ){          
	init_planewithkinks_mapping();
	theExtendedCellSizeAlongStrikeM=theTotalTraceLength/theExtendedColumns;
	pntSrc.dip = 90;  /* Fixed for strikeslip */    
    } else{    
	pntSrc.dip = theExtendedDipDeg; 
	pntSrc.strike = theExtendedStrikeDeg;
    }
    
    pntSrc.sourceFunctionType = theSourceFunctionType;
    pntSrc.T0 = theAverageRisetimeSec;  

    pntSrc.displacement = 
	(double *) malloc ( sizeof( double ) * theNumberOfTimeSteps );  
    if(pntSrc.displacement == NULL){
	fprintf(stdout, "Err alloc displacments");
	return -1;
    } 
  
    pntSrc.dt = theDeltaT;  
    pntSrc.numberOfTimeSteps = theNumberOfTimeSteps;
  
    /* Compute the minimum edge*/
    minEdge = theExtendedMinimumEdge ;   
  
    if ( theExtendedCellSizeAlongStrikeM < minEdge )
	minEdge = theExtendedCellSizeAlongStrikeM;               
    cellPntsStrk = ( int )(theExtendedCellSizeAlongStrikeM/minEdge);
    /*   cellPntsStrk =  1 + ( int )(theExtendedCellSizeAlongStrikeM/minEdge);   */
    minEdgeAlongStrike=theExtendedCellSizeAlongStrikeM/(double)cellPntsStrk;
      
    if ( theExtendedCellSizeDownDipM < minEdge )
	minEdge = theExtendedCellSizeDownDipM;           
    cellPntsDip =  ( int )(theExtendedCellSizeDownDipM / minEdge );
    /*  cellPntsDip = 1 + ( int )(theExtendedCellSizeDownDipM / minEdge ); */
    minEdgeDownDip=theExtendedCellSizeDownDipM/(double)cellPntsDip;

    /* Compute the grids of the fault */    
    theExtendedAlongStrikeDistanceM = 
	theExtendedColumns*theExtendedCellSizeAlongStrikeM;    
    theExtendedDownDipDistanceM =theExtendedRows*theExtendedCellSizeDownDipM;
        
    ndsStrk= 1+(theExtendedAlongStrikeDistanceM/minEdgeAlongStrike);
    ndsDip = 1+(theExtendedDownDipDistanceM/minEdgeDownDip);
    
    grdStrk = dvector(0, ndsStrk);
    grdDip  = dvector(0, ndsDip);  
      
    compute_1D_grid( theExtendedCellSizeAlongStrikeM, strkCells, 
		     cellPntsStrk, minEdgeAlongStrike, grdStrk);
    compute_1D_grid( theExtendedCellSizeDownDipM, dpCells  , 
		     cellPntsDip, minEdgeDownDip, grdDip);  
  
    WAIT;
    if(myID == 0 )
	fprintf(stdout,"\nStart force generation\n");
  
    theForceGenerationTime = -MPI_Wtime();
    theIOInitTime = -MPI_Wtime();
  
    /* Initialize the io */
    /* Array which indicates the loop a force will be written to a file */
    myForcesCycle =  malloc( sizeof(int32_t) * ( myMesh->nharbored ) );  
    for ( iForce = 0 ; iForce < myMesh->nharbored; iForce++ )
	myForcesCycle[ iForce ] = -1;
  
    isforceinprocessor=( char *)calloc((ndsDip*ndsStrk*2 / sizeof(char))+1 , sizeof(char));
    if(isforceinprocessor == NULL){
	fprintf(stdout,"Err allocating isforceinprocessor");
	return -1;
    }

    myNumberOfNodesLoaded =0;
    iForce=0;
    myNumberOfForces=0;


    /* Compute the necessary time steps to represent the source, after these, it will 
       repeat the last value */

    double timeDueWindows=0, timeDueSize=0, tempTime, totalExtendedRaiseTime;
    vector3D_t corner;
  
    /*  Compute the time due to the windows */
    /* Compute the max ruputure time */
    /* Add both */
  
    for ( iWindow=0; iWindow<theNumberOfTimeWindows; iWindow++ )
	timeDueWindows += theWindowDelay [ iWindow ];   
  
    for ( iCorner = 0; iCorner < 2; iCorner++)
	for ( jCorner = 0; jCorner <2; jCorner++){
      
	    corner.x[0] = grdStrk[ iCorner * (ndsStrk-1) ];
	    corner.x[1] = grdDip [ jCorner * (ndsDip-1)  ];
	    corner.x[2] = 0;
      
	    tempTime =compute_initial_time ( corner,
					     theExtendedHypLocCart,
					     theExtendedAverageRuptureVelocity);  
      
	    if ( tempTime > 	timeDueSize  )	
		timeDueSize = tempTime;
      
	}
    
    totalExtendedRaiseTime = timeDueWindows + timeDueSize;
  
    numStepsNecessary = (1.1*totalExtendedRaiseTime) / theDeltaT;
    if( myID == 0){
	fprintf(stdout," \n timeDueWindows = %f" , timeDueWindows );    
	fprintf(stdout," \n timeDueSize = %f" , 	timeDueSize  );  
	fprintf(stdout," \n numStepsNecessary = %d" , 	numStepsNecessary  );
    }
  
    WAIT;
    /* EXITPROGRAM;*/
  
    /* Go through the fault */
    for ( i = 0; i < ndsDip - 1; i++){          

	currentDownDip = i / cellPntsDip;		  
	for ( j = 0; j < ndsStrk - 1; j++){
	    currentAlongStrike = j / cellPntsStrk;        
	    /*Left(0) Right(1) triangle nodes*/      
	    for ( iTri = 0 ; iTri < quadortri; iTri++){	    
	
		p[0]=compute_local_coords(j+i*ndsStrk+iTri,ndsStrk,grdStrk,grdDip);
		p[1]=compute_local_coords(j+(i+iTri)*ndsStrk+1,ndsStrk,grdStrk,grdDip);
		p[2]=compute_local_coords(j+(i+1)*ndsStrk,ndsStrk,grdStrk,grdDip);     	
	  
	
		if ( quadortri == 1 ) /* rectangle */{
		    p[3]=compute_local_coords(j+(i+1)*ndsStrk+1,ndsStrk,grdStrk,grdDip);	

		    localCoords.x[0] = 0.5 * (p[0].x[0]+p[1].x[0]);
		    localCoords.x[1] = 0.5 * (p[0].x[1]+p[2].x[1]);
		    localCoords.x[2] = 0;
	  

		}
		else
		    localCoords = compute_centroid( &p[0] );


	
		if ( theTypeOfSource == PLANE ){ 

		    globalCoords = compute_global_coords(origin,localCoords,pntSrc.dip,
							 0 ,pntSrc.strike);		
		    /* Transform if neccesary to domain coordinates */
		    if ( theRegionAzimuthDeg == 0 )		
			pntSrc.domainCoords=globalCoords;	    	    
		    else
			pntSrc.domainCoords=compute_domain_coords(globalCoords ,
								  theRegionAzimuthDeg);
	  
		}  

		if ( theTypeOfSource == PLANEWITHKINKS ){
		    pntSrc.domainCoords = compute_global_coords_mapping( localCoords );
		}

		/* Check if this force is contained in this processor */
		if ( search_point( pntSrc.domainCoords, &octant ) == 1){ 
		    tag_myForcesCycle ( octant, &pntSrc );	  	  	    	 
		    update_forceinprocessor( iForce, isforceinprocessor,1);
		    myNumberOfForces +=1;
		} 
		iForce+=1;
	    }
	}
    }
  
    WAIT;
    if(myID == 0)
	fprintf(stdout, "\n Total forces = %d ", ndsDip*ndsStrk*2);

  
    memoryRequiredPerNode = sizeof( vector3D_t ) * theNumberOfTimeSteps;
 
    if(memoryRequiredPerNode > BUFFERSIZE ){
	fprintf(stdout,"\n Err not enough memory in the BUFFERSIZE to create the source");
	return -1;
    }
  
    nodesPerCycle = ( int32_t ) floor( BUFFERSIZE / memoryRequiredPerNode );
    memoryPerCycle = nodesPerCycle * memoryRequiredPerNode;
    myNumberOfCycles = (int32_t)floor(myNumberOfNodesLoaded / nodesPerCycle );    
    if ( myNumberOfNodesLoaded % nodesPerCycle != 0 )
	myNumberOfCycles += 1;


    /* go through myForceCycle and put the cycle for each node loaded */
    iCycle = 0;
    iInternal = 0;
    for ( iNode = 0; iNode <  myMesh->nharbored; iNode++){
	if( myForcesCycle[iNode] !=-1 ){
	    myForcesCycle[ iNode ] = iCycle;            
	    iInternal +=1;
	    if(iInternal == nodesPerCycle) {
		iCycle +=1;
		iInternal = 0;
	    }
	}
    }

  
    if(myID == 0){
	fprintf(stdout,"\n Mininum Edge = %f",minEdge);
	fprintf(stdout,"\n Minimum Memory Requiered to build the source = %d\n",memoryRequiredPerNode);

    }

    print_header_myForces(fptmpsrc);  

    WAIT;
    theIOInitTime += MPI_Wtime();


    if(myID==0) fprintf(stdout,"\n\nI/O Force Intitialization done: Time = %.2f seconds", theIOInitTime);      

    /* END I/O initialization*/
  
    myNumberOfNodesLoaded =0;

    for ( iCycle = 0; iCycle < myNumberOfCycles; iCycle++){ 
	iForce = 0;
	for ( i = 0; i < ndsDip - 1; i++){          
	    currentDownDip = i / cellPntsDip;	    
	    for ( j = 0; j < ndsStrk- 1; j++){
		currentAlongStrike = j / cellPntsStrk;       	
		/*Left(0) Right(1) triangle nodes*/     
		for ( iTri = 0 ; iTri < quadortri; iTri++){	    	  
		    if( is_forceinprocessor( iForce, isforceinprocessor ) == 1 ){	    

			p[0]=compute_local_coords(j+i*ndsStrk+iTri,ndsStrk,grdStrk,grdDip);	 
			p[1]=compute_local_coords(j+(i+iTri)*ndsStrk+1,ndsStrk,grdStrk,grdDip);      
			p[2]=compute_local_coords(j+(i+1)*ndsStrk,ndsStrk,grdStrk,grdDip);	   	    

			if ( quadortri == 1 ) /* rectangle */{
			    p[3]=compute_local_coords(j+(i+1)*ndsStrk+1,ndsStrk,grdStrk,grdDip);	
	      
			    localCoords.x[0] = 0.5 * (p[0].x[0]+p[1].x[0]);
			    localCoords.x[1] = 0.5 * (p[0].x[1]+p[2].x[1]);
			    localCoords.x[2] = 0;

			    area = fabs( (p[0].x[0]-p[1].x[0])*(p[0].x[1]-p[2].x[1]) );
	      
			    /*	      if (myID == 0 ){
			      fprintf(stdout," \n Centroide x=%f  y=%f z=%f,", 
			      localCoords.x[0],localCoords.x[1],localCoords.x[2]);
			      fprintf(stdout," \n area = %f ", area );
			      EXITPROGRAM;*/
			    /* }*/
	      
			}
			else{
			    area=compute_area( &p[0] );
	    
			    localCoords = compute_centroid( &p[0] );

			    /*      if (myID == 0 ){
				    fprintf(stdout," \n Centroide x=%f  y=%f z=%f,", 
				    localCoords.x[0],localCoords.x[1],localCoords.x[2]);
				    fprintf(stdout," \n area = %f ", area );
				    EXITPROGRAM;
				    }*/
	      
			}

			pntSrc.localCoords = localCoords;
	    
			if ( theTypeOfSource == PLANE ){
			    globalCoords = compute_global_coords( origin, localCoords,pntSrc.dip ,
								  0 ,theExtendedStrikeDeg);
			    /* Transform if neccesary to domain coordinates just for plane*/	  
			    if ( theRegionAzimuthDeg == 0 )		
				pntSrc.domainCoords = globalCoords;	    	    
			    else{
				pntSrc.domainCoords =  
				    compute_domain_coords( globalCoords, theRegionAzimuthDeg );
				pntSrc.strike = theExtendedStrikeDeg - theRegionAzimuthDeg;
			    }	 
			}
	    	    
	    
			if ( theTypeOfSource == PLANEWITHKINKS ){
			    pntSrc.strike = compute_strike_planewithkinks( localCoords );            
			    pntSrc.domainCoords = compute_global_coords_mapping( localCoords );
			}
	    


			if ( search_point( pntSrc.domainCoords, &octant ) == 1){

	      
			    if( myID == 0 && iCycle == 0) /* useful just with one processor */
				fprintf(fpdsrc,"\n%f %f %f %f %f %f %f %f",
					pntSrc.domainCoords.x[0], pntSrc.domainCoords.x[1],
					pntSrc.domainCoords.x[2], pntSrc.strike,
					pntSrc.dip, pntSrc.rake, pntSrc.maxSlip,
					pntSrc.tinit);    
	      
			    /* This is the part that has to be reviewed */
			    update_point_source( &pntSrc,currentDownDip,currentAlongStrike,area);
			    load_myForces_with_point_source(octant,&pntSrc,iCycle,iForce);

	
			}	 
		    }
		    iForce += 1;
		}         
	    }
	}

	fill_myForces_cycle(iCycle);

	if ( theSourceIsFiltered == 1 ) filter_myForce();
	print_myForces_filepercycle(fileNameRoot,iCycle);     
	free_myForces();

    }

    free(pntSrc.displacement);
  
    /* Concatenation process */
    WAIT;
    theForceGenerationTime += MPI_Wtime();
    theConcatTime = -MPI_Wtime();
    if( myID == 0 ){
	fprintf(stdout,"\nSource files generated: Time = %.2f seconds", theForceGenerationTime);
	fprintf(stdout,"\nConcatenation process starts");
    }

    /* EXITPROGRAM; */
    if(myNumberOfCycles != 0)
	concat_filespercycle_in_fptmpsrc(fileNameRoot,fptmpsrc,myNumberOfCycles);  

    WAIT;
    theConcatTime += MPI_Wtime();
  
  
    if( myID ==0 )
	fprintf(stdout,"\nConcatenation process done: Time = %.2f seconds", theConcatTime );

    /*  if(myID == 0){      
	fprintf(stdout,"\n cycles = %d ", myNumberOfCycles);
	}*/
      
    /*  WAIT;*/

    if(myID == 0)      
	fprintf(stdout,"\n cycles = %d ", myNumberOfCycles);
      
    /*EXITPROGRAM;*/

    free(myForces);
    
    if(myID == 0 )fclose(fpdsrc);
  
    fclose(fptmpsrc);
  
    return 1;  
}


/**
 * Use various global configuration variables to compute the coordinates
 * of a point source.
 */
static void 
fill_point_source_coordinates (ptsrc_t* ps)
{
    /* Coordinates */    
    /*  ps->globalCoords.x[0] = ((theHypocenterLatDeg - theRegionOriginLatDeg) */
    /* 			     *DIST1LAT); */
    /*     ps->globalCoords.x[1] = ((theHypocenterLongDeg - theRegionOriginLongDeg) */
    /* 			     *DIST1LON); */
    /*     ps->globalCoords.x[2] = theHypocenterDepthM - theRegionDepthShallowM */;

    /* Transform if neccesary to the domain coordinates */
    /*   if (theRegionAzimuthDeg == 0) { */
    /* 	ps->domainCoords.x[0] = ps->globalCoords.x[0]; */
    /* 	ps->domainCoords.x[1] = ps->globalCoords.x[1]; */
    /* 	ps->domainCoords.x[2] = ps->globalCoords.x[2]; */
    /*     } else { */
    /* 	ps->domainCoords */
    /* 	    = compute_domain_coords (ps->globalCoords, theRegionAzimuthDeg); */
    /*     } */

    if( theLonlatOrCartesian == 1 ){

	ps->domainCoords.x[0] = theHypocenterLatDeg;
	ps->domainCoords.x[1] = theHypocenterLongDeg;
	ps->domainCoords.x[2] = theHypocenterDepthM;

    }

  
    if( theLonlatOrCartesian == 0 ){

	printf("\n %f %f", theHypocenterLongDeg,theHypocenterLatDeg);

	ps->domainCoords =compute_domain_coords_linearinterp(theHypocenterLongDeg,
							     theHypocenterLatDeg,
							     theSurfaceCornersLong , 
							     theSurfaceCornersLat,
							     theRegionLengthEastM,
							     theRegionLengthNorthM );


	/*    printf("\n %f %f",ps->domainCoords.x[0],ps->domainCoords.x[1]);
	      printf("\n %f %f",theSurfaceCornersLong[0],theSurfaceCornersLat[0]);
	      printf("\n %f %f",theSurfaceCornersLong[1],theSurfaceCornersLat[1]);
	      printf("\n %f %f",theSurfaceCornersLong[2],theSurfaceCornersLat[2]);
	      printf("\n %f %f",theSurfaceCornersLong[3],theSurfaceCornersLat[3]);
	      exit(1); */

	ps->domainCoords.x[2] = theHypocenterDepthM;

    }



  
}

/**
 * Use various global configuration variables to compute the strike
 * of a point source.
 * 
 *
 * Note: 
 *       The strike is computed asumming linear interpolation for the domain
 *       A point having the same longitude and an increment in the latitude
 *       is used to obtain a vector pointing the north. After,  the rake 
 *       is computed considering the vector pointing towards north.
 */
static void 
compute_point_source_strike (ptsrc_t* ps)
{

    vector3D_t pivot, pointInNorth, unitVec;

    double norm,fi ;


    if( theLonlatOrCartesian == 1 )return;

  
    if( theLonlatOrCartesian == 0 ){

	printf("\nin strike correction = %f %f", theHypocenterLongDeg,theHypocenterLatDeg);

	pivot = compute_domain_coords_linearinterp(theHypocenterLongDeg,
						   theHypocenterLatDeg,
						   theSurfaceCornersLong , 
						   theSurfaceCornersLat,
						   theRegionLengthEastM,
						   theRegionLengthNorthM );
    
	pointInNorth = compute_domain_coords_linearinterp(theHypocenterLongDeg,
							  theHypocenterLatDeg+.1,
							  theSurfaceCornersLong , 
							  theSurfaceCornersLat ,
							  theRegionLengthEastM,
							  theRegionLengthNorthM );
    

	/* Compute Unit Vector */
	unitVec.x[1]=pointInNorth.x[1]-pivot.x[1];
	unitVec.x[0]=pointInNorth.x[0]-pivot.x[0];

	norm = pow( pow( unitVec.x[0], 2)+ pow(unitVec.x[1],2), .5);
   
	unitVec.x[1]= unitVec.x[1]/norm;
	unitVec.x[0]= unitVec.x[0]/norm;
    
	/* Compute the Angle North-X axis */
	fi = atan( unitVec.x[0]/unitVec.x[1]);
    
	if(  unitVec.x[1] < 0 ) /* in rad*/
	    fi = fi + PI;

	/* Compute the strike */ 
    
	ps->strike =90+ ps->strike-( 180*fi/PI);
   
    } 

}

/**
 * Compute force due to a dipole.
 *  
 * \return 1 OK, -1 on error.
 */
static int compute_myForces_point (const char *physicsin) 
{
    int32_t   iForce;
    octant_t* octant;
    ptsrc_t   pntSrc;
    FILE*     fpdsrc   = NULL;
    FILE*     fptmpsrc = NULL;
    char      fileNameRoot[256];
    int	      rt;
  
    fill_point_source_coordinates (&pntSrc);
  
    rt = open_filescomputemyForces (physicsin, &fpdsrc, &fptmpsrc,fileNameRoot);
    if (rt == -1) return -1;
  
    if( myID == 0 ){
	fprintf(fpdsrc," %lf %lf %lf",pntSrc.domainCoords.x[0]
		,pntSrc.domainCoords.x[1]
		,pntSrc.domainCoords.x[2]);
	close_file (&fpdsrc);
    }
  
    /* search for source in this processor */
    if ( search_point (pntSrc.domainCoords, &octant) != 1) {
	close_file (&fptmpsrc);
	/* not in this processor, nothing else to do here!, return OK */
	return 1;
    }
  
    isforceinprocessor = (char*)calloc (1, sizeof (char));
  
    if (isforceinprocessor == NULL) {
	fprintf (stdout,"Err allocating forEcinprocessbis");
	return -1;
    }
  
    iForce=0;
    update_forceinprocessor (0, isforceinprocessor,1);

    /* Array which indicates the cycle a force will written to a file in
     * this case is one cycle
     */
    myForcesCycle = malloc (sizeof(int32_t) * (myMesh->nharbored));
    for (iForce = 0 ; iForce < myMesh->nharbored; iForce++) {
	myForcesCycle[ iForce ] = -1;
    }
  
    tag_myForcesCycle (octant, &pntSrc);
    myNumberOfForces += 1;
  
    /* This variable is introduced thinking in large faults. It is approx.
       the time it takes to reach the maximum value of slip. We compute every
       thing in this type of source*/
    numStepsNecessary = theNumberOfTimeSteps;

    /* We have to load the force vector */

    /*first we build the source */  
    pntSrc.strike             = theSourceStrikeDeg;
  
    compute_point_source_strike (&pntSrc);

    pntSrc.dip                = theSourceDipDeg;
    pntSrc.rake               = theSourceRakeDeg;
    pntSrc.sourceFunctionType = theSourceFunctionType;
    pntSrc.T0                 = theAverageRisetimeSec;
    pntSrc.delayTime          = 0;
    pntSrc.maxSlip            = 1;
    pntSrc.M0                 = theMomentAmplitude;
    pntSrc.displacement
	= (double*)calloc (theNumberOfTimeSteps, sizeof (double));
  
    if (pntSrc.displacement == NULL) {
	fprintf(stderr, "Err allocating pntSrc.displacement array");
	return -1;
    }

    pntSrc.dt = theDeltaT;  
    pntSrc.numberOfTimeSteps = theNumberOfTimeSteps;
  
    int32_t memoryRequiredPerNode, nodesPerCycle,myNumberOfCycles, iCycle,
	iInternal, iNode, memoryPerCycle;
  
    memoryRequiredPerNode = 3 * sizeof (vector3D_t) * theNumberOfTimeSteps;
    if (memoryRequiredPerNode > BUFFERSIZE) {
	fputs ("\n ERROR: Not enough memory in the BUFFERSIZE to create the "
	       "source", stderr);
	return -1;
    }
  
    /* note: these are done using integer arithmetic */
    nodesPerCycle    = BUFFERSIZE / memoryRequiredPerNode;
    memoryPerCycle   = nodesPerCycle * memoryRequiredPerNode;
    myNumberOfCycles = ((myNumberOfNodesLoaded - 1) / nodesPerCycle) + 1;
  
    /* go through myForceTag and put the cycle for each node loaded */
    iCycle    = 0;
    iInternal = 0;
  
    for (iNode = 0; iNode <  myMesh->nharbored; iNode++) {
	if (myForcesCycle[iNode] !=-1) {
	    myForcesCycle[ iNode ] = iCycle;            
	    iInternal +=1;
      
	    if (iInternal == nodesPerCycle) {
		iCycle +=1;
		iInternal = 0;
	    }
	}
    }

    /* fprintf(stdout,"\n num time steps =%d", theNumberOfTimeSteps); */
            
    compute_source_function (&pntSrc);
    myNumberOfNodesLoaded = 0;
    load_myForces_with_point_source (octant, &pntSrc, 0, 0);
    print_header_myForces (fptmpsrc);
    print_myForces_transposed (fptmpsrc);

    /* cleanup */
    free (pntSrc.displacement);
    close_file (&fptmpsrc);

    return 1;
}


/* 
 * compute_myForces_srfh: compute the force vector for each 
 *                        node in myMesh
 *
 *                return -1 fail 1 ok
 */
static int  compute_myForces_srfh(const char *physicsin){
  
    FILE *fpdsrc, *fptmpsrc;
    int32_t iSrc,iForce=0;
    octant_t *octant;
    ptsrc_t pntSrc;          
  
    /* i/o related vars */
    int32_t nodesPerCycle, memoryRequiredPerNode, memoryPerCycle;
    int32_t iCycle, iInternal, iNode; 

    /* convinient variables */
    char fileNameRoot[256]; 

    numStepsNecessary = theNumberOfTimeSteps;
  
    theForceGenerationTime = -MPI_Wtime();
    theIOInitTime          = -MPI_Wtime();
 
    if(open_filescomputemyForces(physicsin,&fpdsrc,&fptmpsrc,fileNameRoot)==-1)
	return -1;

    pntSrc.displacement = (double*)malloc(sizeof(double)*theNumberOfTimeSteps );  
    if(pntSrc.displacement == NULL){
	fprintf(stdout, "Err alloc displacments");
	return -1;
    } 
  
    pntSrc.dt = theDeltaT;  
    pntSrc.numberOfTimeSteps = theNumberOfTimeSteps;
      
    WAIT;

    if(myID == 0 ){
	fprintf(stdout,"\n\nForce generation:starts\n");
    }
    
    theForceGenerationTime = -MPI_Wtime();
    theIOInitTime = -MPI_Wtime();
  
    /* Initialize the IO */
    /* Array which indicates the loop a force will be written to a file */
    myForcesCycle =  malloc( sizeof(int32_t) * ( myMesh->nharbored ) );  
    for ( iForce = 0 ; iForce < myMesh->nharbored; iForce++ )
	myForcesCycle[ iForce ] = -1;
  
    isforceinprocessor=( char *)calloc(theNumberOfPointSources*2 / sizeof(char)+1 , 
				       sizeof(char));  
    if(isforceinprocessor == NULL){
	fprintf(stdout,"Err allocating isforceinprocessor");
	return -1;
    }
  
    myNumberOfNodesLoaded =0;
    iForce=0;
    myNumberOfForces=0;  
  
    /* Go through the fault */
    for ( iSrc = 0; iSrc <theNumberOfPointSources ; iSrc++){          
    
	pntSrc.domainCoords = compute_domain_coords_linearinterp(theSourceLonArray[iSrc],
								 theSourceLatArray[iSrc],
								 theSurfaceCornersLong , 
								 theSurfaceCornersLat,
								 theRegionLengthEastM,
								 theRegionLengthNorthM );    
	pntSrc.domainCoords.x[2]= theSourceDepthArray[iSrc];

	/* Check if this force is contained in this processor */
	if ( search_point( pntSrc.domainCoords, &octant ) == 1){ 
	    tag_myForcesCycle ( octant, &pntSrc );	  	  	    	 
	    update_forceinprocessor( iForce, isforceinprocessor,1);
	    myNumberOfForces +=1;
	} 
	iForce+=1;    
    }
  
    WAIT;
  
    if(myID == 0) fprintf(stdout, "\nTotal forces = %d ", iForce);

    memoryRequiredPerNode = sizeof( vector3D_t ) * theNumberOfTimeSteps;
  
    if(memoryRequiredPerNode > BUFFERSIZE ){
	fprintf(stdout,"\n Memory in the BUFFERSIZE source: compute_myForces_srfh");
	return -1;
    }

    nodesPerCycle = ( int32_t ) floor( BUFFERSIZE / memoryRequiredPerNode );
    memoryPerCycle = nodesPerCycle * memoryRequiredPerNode;
    myNumberOfCycles = (int32_t)floor(myNumberOfNodesLoaded / nodesPerCycle );    
    if ( myNumberOfNodesLoaded % nodesPerCycle != 0 ) myNumberOfCycles += 1;
  
    /* go through myForceCycle and put the cycle for each node loaded */
    iCycle = 0;
    iInternal = 0;
    for ( iNode = 0; iNode <  myMesh->nharbored; iNode++){
	if( myForcesCycle[iNode] !=-1 ){
	    myForcesCycle[ iNode ] = iCycle;            
	    iInternal +=1;
	    if(iInternal == nodesPerCycle) {
		iCycle +=1;
		iInternal = 0;
	    }
	}
    }

    print_header_myForces(fptmpsrc);  

    WAIT;

    theIOInitTime += MPI_Wtime();

    if(myID==0)
	fprintf(stdout,"\n I/O Force Intitialization ...%.2f s", theIOInitTime);      
  
    WAIT; 

    /* END I/O initialization*/
    if(myID == 0)
	fprintf(stdout,"\n Source files ...........");      
  

    myNumberOfNodesLoaded =0;
  
    for ( iCycle = 0; iCycle < myNumberOfCycles; iCycle++){     
	for ( iSrc = 0; iSrc < theNumberOfPointSources ; iSrc++){          
	    if( is_forceinprocessor( iSrc, isforceinprocessor ) == 1 ){	    
	
		pntSrc.domainCoords = compute_domain_coords_linearinterp(theSourceLonArray[iSrc],
									 theSourceLatArray[iSrc],
									 theSurfaceCornersLong , 
									 theSurfaceCornersLat,
									 theRegionLengthEastM,
									 theRegionLengthNorthM );
		pntSrc.domainCoords.x[2]= theSourceDepthArray[iSrc];       	      	       
	
		if ( search_point( pntSrc.domainCoords, &octant ) == 1){
	  
		    update_point_source_srfh( &pntSrc,iSrc);
	  
		    load_myForces_with_point_source(octant,&pntSrc,iCycle,iSrc);
	  
		    if( myID == 0 && iCycle == 0){ /*useful just with one processor*/
			fprintf(fpdsrc,"\n source=%f %f %f %f %f %f %f %f ",
				pntSrc.domainCoords.x[0], pntSrc.domainCoords.x[1],
				pntSrc.domainCoords.x[2], pntSrc.strike,
				pntSrc.dip, pntSrc.rake, pntSrc.maxSlip,
				pntSrc.delayTime);   
			/*	  fprintf(stdout,"\n%f %f %f %f %f %f %f %f ",
			  pntSrc.domainCoords.x[0], pntSrc.domainCoords.x[1],
			  pntSrc.domainCoords.x[2], pntSrc.strike,
			  pntSrc.dip, pntSrc.rake, pntSrc.maxSlip,
			  pntSrc.delayTime); }*/}
	  
		}	



	    }
	}
    
	fill_myForces_cycle(iCycle);
    
	if ( theSourceIsFiltered == 1 ) filter_myForce();
	print_myForces_filepercycle(fileNameRoot,iCycle);     
	free_myForces();    
    }

    free(pntSrc.displacement);
  
    /* Concatenation process */
    WAIT;
    theForceGenerationTime += MPI_Wtime();
    theConcatTime = -MPI_Wtime();
    if( myID == 0 ){
	fprintf(stdout,"done...%.2f s", theForceGenerationTime);
	fprintf(stdout,"\n Concatenation process...");
    }

    if(myNumberOfCycles != 0)
	concat_filespercycle_in_fptmpsrc(fileNameRoot,fptmpsrc,myNumberOfCycles);
  
    WAIT;
    theConcatTime += MPI_Wtime();
  
  
    if( myID ==0 )
	fprintf(stdout,"done...%.2f s", theConcatTime );
  
    free(myForces);
    
    if(myID == 0 )fclose(fpdsrc);
  
    fclose(fptmpsrc);
  
    return 1;  

}



/* 
 * compute_myForces_green: compute the force vector for each 
 *                        node in myMesh
 *
 *                return -1 fail 1 ok
 */
static int  compute_myForces_green(const char *physicsin){
  
    FILE *fpdsrc, *fptmpsrc;
    int32_t iSrc,iForce=0;
    octant_t *octant;
    ptsrc_t pntSrc;          
  
    /* i/o related vars */
    int32_t nodesPerCycle, memoryRequiredPerNode, memoryPerCycle;
    int32_t iCycle, iInternal, iNode; 

    /* convinient variables */
    char fileNameRoot[256]; 

    numStepsNecessary = theNumberOfTimeSteps;
  
    theForceGenerationTime = -MPI_Wtime();
    theIOInitTime          = -MPI_Wtime();
 
    //    fprintf(stdout,"Entro a compute_myForces_green");
    fflush(stdout);

    if(open_filescomputemyForces(physicsin,&fpdsrc,&fptmpsrc,fileNameRoot)==-1)
	return -1;

    pntSrc.displacement = (double*)malloc(sizeof(double)*theNumberOfTimeSteps );  
    if(pntSrc.displacement == NULL){
	fprintf(stdout, "Err alloc displacments");
	return -1;
    } 
  
    pntSrc.dt = theDeltaT;  
    pntSrc.numberOfTimeSteps = theNumberOfTimeSteps;
      
    WAIT;

    if(myID == 0 ){
	fprintf(stdout,"\n\nForce generation:starts\n");
    }
    
    theForceGenerationTime = -MPI_Wtime();
    theIOInitTime = -MPI_Wtime();
  
    /* Initialize the IO */
    /* Array which indicates the loop a force will be written to a file */
    myForcesCycle =  malloc( sizeof(int32_t) * ( myMesh->nharbored ) );  
    for ( iForce = 0 ; iForce < myMesh->nharbored; iForce++ )
	myForcesCycle[ iForce ] = -1;
  
    isforceinprocessor=( char *)calloc(theNumberOfPointSources*2 / sizeof(char)+1 , 
				       sizeof(char));  
    if(isforceinprocessor == NULL){
	fprintf(stdout,"Err allocating isforceinprocessor");
	return -1;
    }
  
    myNumberOfNodesLoaded =0;
    iForce=0;
    myNumberOfForces=0;  
  
    /* Go through the fault */
    for ( iSrc = 0; iSrc <theNumberOfPointSources ; iSrc++){          
    
	pntSrc.domainCoords = compute_domain_coords_linearinterp(theSourceLonArray[iSrc],
								 theSourceLatArray[iSrc],
								 theSurfaceCornersLong , 
								 theSurfaceCornersLat,
								 theRegionLengthEastM,
								 theRegionLengthNorthM );    
	pntSrc.domainCoords.x[2]= theSourceDepthArray[iSrc];

	/* Check if this force is contained in this processor */
	if ( search_point( pntSrc.domainCoords, &octant ) == 1){ 
	    tag_myForcesCycle ( octant, &pntSrc );	  	  	    	 
	    update_forceinprocessor( iForce, isforceinprocessor,1);
	    myNumberOfForces +=1;
	} 
	iForce+=1;    
    }
  
    WAIT;
  
    if(myID == 0) fprintf(stdout, "\nTotal forces = %d ", iForce);

    memoryRequiredPerNode = sizeof( vector3D_t ) * theNumberOfTimeSteps;
  
    if(memoryRequiredPerNode > BUFFERSIZE ){
	fprintf(stdout,"\n Memory in the BUFFERSIZE source: compute_myForces_srfh");
	return -1;
    }

    nodesPerCycle = ( int32_t ) floor( BUFFERSIZE / memoryRequiredPerNode );
    memoryPerCycle = nodesPerCycle * memoryRequiredPerNode;
    myNumberOfCycles = (int32_t)floor(myNumberOfNodesLoaded / nodesPerCycle );    
    if ( myNumberOfNodesLoaded % nodesPerCycle != 0 ) myNumberOfCycles += 1;
  
    /* go through myForceCycle and put the cycle for each node loaded */
    iCycle = 0;
    iInternal = 0;
    for ( iNode = 0; iNode <  myMesh->nharbored; iNode++){
	if( myForcesCycle[iNode] !=-1 ){
	    myForcesCycle[ iNode ] = iCycle;            
	    iInternal +=1;
	    if(iInternal == nodesPerCycle) {
		iCycle +=1;
		iInternal = 0;
	    }
	}
    }

    print_header_myForces(fptmpsrc);  

    WAIT;

 
    theIOInitTime += MPI_Wtime();

    if(myID==0)
	fprintf(stdout,"\n I/O Force Intitialization ...%.2f s", theIOInitTime);      
  
    WAIT; 

    /* END I/O initialization*/
    if(myID == 0)
	fprintf(stdout,"\n Source files ...........");      
  
    myNumberOfNodesLoaded =0;
  
    //fprintf(stdout,"\n My numberOfCycles=%d\n", myNumberOfCycles);
    //    fprintf(stdout,"\n My theNumberOfPointSources=%d\n", theNumberOfPointSources);


    for ( iCycle = 0; iCycle < myNumberOfCycles; iCycle++){     
	for ( iSrc = 0; iSrc < theNumberOfPointSources ; iSrc++){          
	    if( is_forceinprocessor( iSrc, isforceinprocessor ) == 1 ){	    
	

		pntSrc.domainCoords = compute_domain_coords_linearinterp(theSourceLonArray[iSrc],
									 theSourceLatArray[iSrc],
									 theSurfaceCornersLong , 
									 theSurfaceCornersLat,
									 theRegionLengthEastM,
									 theRegionLengthNorthM );
		pntSrc.domainCoords.x[2]= theSourceDepthArray[iSrc];       	      	       

		//		fprintf(stdout,"\n myID %d  Lon=%f Lat=%f \n",myID,theSourceLonArray[iSrc],theSourceLatArray[iSrc]);


		//fprintf(stdout,"\n %f %f %f \n",pntSrc.domainCoords.x[0],pntSrc.domainCoords.x[1],pntSrc.domainCoords.x[2]);



		if ( search_point( pntSrc.domainCoords, &octant ) == 1){
	  
		    update_point_source_green( &pntSrc,iSrc);
	  
		    load_myForces_with_point_source_green(octant,&pntSrc,iCycle,iSrc);
	  
		    if( myID == 0 && iCycle == 0){ /*useful just with one processor*/
		      fprintf(stdout,"\n%f %f %f %f %f %f %f %f ",
			      pntSrc.domainCoords.x[0], pntSrc.domainCoords.x[1],
			      pntSrc.domainCoords.x[2], pntSrc.weightX,
			      pntSrc.weightY, pntSrc.weightZ,pntSrc.delayTime);   
		    }	  
		}	

	    }
	}
    
	fill_myForces_cycle(iCycle);
    
	if ( theSourceIsFiltered == 1 ) filter_myForce();
	print_myForces_filepercycle(fileNameRoot,iCycle);     
	free_myForces();    
    }

    free(pntSrc.displacement);
  
    /* Concatenation process */
    WAIT;
    theForceGenerationTime += MPI_Wtime();
    theConcatTime = -MPI_Wtime();
    if( myID == 0 ){
	fprintf(stdout,"done...%.2f s", theForceGenerationTime);
	fprintf(stdout,"\n Concatenation process...");
    }

    if(myNumberOfCycles != 0)
	concat_filespercycle_in_fptmpsrc(fileNameRoot,fptmpsrc,myNumberOfCycles);
  
    WAIT;
    theConcatTime += MPI_Wtime();
  
  
    if( myID ==0 )
	fprintf(stdout,"done...%.2f s", theConcatTime );
  
    free(myForces);
    
    if(myID == 0 )fclose(fpdsrc);
  
    fclose(fptmpsrc);
  
    return 1;  

}



/**
 * Broadcast plane parameters.
 *
 * \note It initializes various global variables (too many to list here).
 */
static void
broadcast_plane_parameters (void)
{
    /* The following #defines would not be necessary if the compiler were
     * to support automatic (stack) array allocation by specifying their
     * length with vars declared as constant.  We undefine these macros at
     * the end of this function.
     */
#define PLANE_PARAMS_DN		13
#define PLANE_PARAMS_IN		 2
  
    int     i;
    double  plane_params_d[PLANE_PARAMS_DN];
    int32_t plane_params_i[PLANE_PARAMS_IN];

    plane_params_d[ 0] = theExtendedCellSizeAlongStrikeM;
    plane_params_d[ 1] = theExtendedCellSizeDownDipM;
    plane_params_d[ 2] = theExtendedLatDeg;
    plane_params_d[ 3] = theExtendedLongDeg;
    plane_params_d[ 4] = theExtendedDepthM;
    plane_params_d[ 5] = theExtendedAlongStrikeDistanceM;
    plane_params_d[ 6] = theExtendedDownDipDistanceM;
    plane_params_d[ 7] = theExtendedHypocenterAlongStrikeM;
    plane_params_d[ 8] = theExtendedHypocenterDownDipM;
    plane_params_d[ 9] = theExtendedAverageRuptureVelocity;
    plane_params_d[10] = theExtendedStrikeDeg;
    plane_params_d[11] = theExtendedDipDeg;
    plane_params_d[12] = theExtendedMinimumEdge;

    plane_params_i[0] = theExtendedColumns;
    plane_params_i[1] = theExtendedRows;

    MPI_Bcast (plane_params_d, PLANE_PARAMS_DN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (plane_params_i, PLANE_PARAMS_IN, MPI_INT,    0, MPI_COMM_WORLD);

    if (0 != myID) {
	theExtendedCellSizeAlongStrikeM   = plane_params_d[ 0];
	theExtendedCellSizeDownDipM	  = plane_params_d[ 1];
	theExtendedLatDeg		  = plane_params_d[ 2];
	theExtendedLongDeg		  = plane_params_d[ 3];
	theExtendedDepthM		  = plane_params_d[ 4];
	theExtendedAlongStrikeDistanceM   = plane_params_d[ 5];
	theExtendedDownDipDistanceM	  = plane_params_d[ 6];
	theExtendedHypocenterAlongStrikeM = plane_params_d[ 7];
	theExtendedHypocenterDownDipM	  = plane_params_d[ 8];
	theExtendedAverageRuptureVelocity = plane_params_d[ 9];
	theExtendedStrikeDeg		  = plane_params_d[10];
	theExtendedDipDeg		  = plane_params_d[11];
	theExtendedMinimumEdge		  = plane_params_d[12];

	theExtendedColumns = plane_params_i[0];
	theExtendedRows    = plane_params_i[1];
    }

 

    /* broadcast slip and rake matrices */
    if (myID != 0) {
	/* other processes need to allocate memory, PE 0 has already done it */
	theSlipMatrix = malloc (sizeof(double **) * theNumberOfTimeWindows);
	theRakeMatrix = malloc (sizeof(double **) * theNumberOfTimeWindows);
	
	if ( (theSlipMatrix == NULL) || (theRakeMatrix == NULL) ) {
	    solver_abort ("broacast_plane_parameters", "malloc failed "
			  "allocating slip and/or rake  matrices", NULL);
	}
	
	for (i = 0; i < theNumberOfTimeWindows; i++) {
	    theSlipMatrix[i]
		= dmatrix (0, theExtendedRows - 1, 0, theExtendedColumns - 1);
	    theRakeMatrix[i]
		= dmatrix (0, theExtendedRows - 1, 0, theExtendedColumns - 1);
	}
    } /* if (myID != 0) */




    for (i = 0; i < theNumberOfTimeWindows; i++) {

	/* since dmalloc allocates a contiguous chunk for the matrix,
	 * we can get away with specifying the address of the first row
	 * and then the size of the matrix for the broadcast.
	 */
	MPI_Bcast (theSlipMatrix[i][0], theExtendedRows * theExtendedColumns,
		   MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Bcast (theRakeMatrix[i][0], theExtendedRows * theExtendedColumns, 
		   MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
 

    MPI_Bcast (theSurfaceCornersLat,  4, MPI_DOUBLE, 0, MPI_COMM_WORLD);      
    MPI_Bcast (theSurfaceCornersLong, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);




#undef PLANE_PARAMS_DN
#undef PLANE_PARAMS_IN
}

/**
 * Broadcast srfh parameters.
 *
 * \note It initializes various global variables (too many to list here).
 */
static void
broadcast_srfh_parameters (void)
{
    int32_t iSource;

    MPI_Bcast(&theNumberOfPointSources,  1, MPI_INT , 0, MPI_COMM_WORLD);
      
    if ( myID != 0) {
	
	theSourceLonArray   = malloc( sizeof( double) * theNumberOfPointSources );
	theSourceLatArray   = malloc( sizeof( double) * theNumberOfPointSources );
	theSourceDepthArray = malloc( sizeof( double) * theNumberOfPointSources );
	theSourceAreaArray  = malloc( sizeof( double) * theNumberOfPointSources );
	theSourceStrikeArray= malloc( sizeof( double) * theNumberOfPointSources );
	theSourceDipArray   = malloc( sizeof( double) * theNumberOfPointSources );
	theSourceRakeArray  = malloc( sizeof( double) * theNumberOfPointSources );	
	theSourceSlipArray  = malloc( sizeof( double) * theNumberOfPointSources );	
    
	theSourceSlipArray  = malloc( sizeof( double) * theNumberOfPointSources );
	theSourceNt1Array   = malloc( sizeof( int )   * theNumberOfPointSources );
	theSourceTinitArray = malloc( sizeof( double) * theNumberOfPointSources );
	theSourceDtArray    = malloc( sizeof( double) * theNumberOfPointSources );
    
	theSourceSlipFunArray  = malloc( sizeof( double * ) * theNumberOfPointSources );

	if ( (theSourceLonArray    == NULL) || (theSourceLatArray    == NULL) || 
	     (theSourceDepthArray  == NULL) || (theSourceStrikeArray == NULL) ||
	     (theSourceStrikeArray == NULL) || (theSourceDipArray    == NULL) ||
	     (theSourceRakeArray == NULL)   ) {
	    ABORT_PROGRAM (" Error source srfh matrices:broadcast_srfh_parameters ");
	}

	if ( ( theSourceNt1Array  == NULL) || ( theSourceTinitArray    == NULL) || 
	     ( theSourceDtArray  == NULL) || (theSourceSlipFunArray == NULL))
	    ABORT_PROGRAM (" Error source srfh matrices:broadcast_srfh_parameters ");
    
    }

    WAIT;
  
    MPI_Bcast (theSourceLonArray,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);	     
    MPI_Bcast (theSourceLatArray,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);
    MPI_Bcast (theSourceDepthArray,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);  
    MPI_Bcast (theSourceAreaArray,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);
    MPI_Bcast (theSourceStrikeArray,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);      
    MPI_Bcast (theSourceDipArray,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);
    MPI_Bcast (theSourceRakeArray,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);
    MPI_Bcast (theSourceSlipArray,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);  
    MPI_Bcast (theSourceNt1Array,theNumberOfPointSources , MPI_INT, 0,
	       MPI_COMM_WORLD);
    MPI_Bcast (theSourceTinitArray,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);
    MPI_Bcast (theSourceDtArray,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);

    if(myID == 0)
	fprintf(stdout,"\n  Broadcast small arrays source srfh done %d \n", myID);
  
    fflush(stdout);


    WAIT;

    if(myID != 0){

	for ( iSource = 0; iSource < theNumberOfPointSources; iSource++ ){ 
	    /* double check*/
	    if( theSourceNt1Array[iSource] == 0){
		fprintf(stdout,"\n error the size of the slipfunc = 0 in  source%d \n",iSource); 
		fflush(stdout);
		ABORTEXIT;

	    }

	    theSourceSlipFunArray[ iSource ] =
		malloc( sizeof( double ) * theSourceNt1Array[iSource] );
      
	}
    }

    WAIT;
  
    for ( iSource = 0; iSource < theNumberOfPointSources; iSource++ ){     
	MPI_Bcast( theSourceSlipFunArray[iSource],
		   theSourceNt1Array[iSource], MPI_DOUBLE, 0,
		   MPI_COMM_WORLD);    
    }
  
    WAIT;
  
    MPI_Bcast (theSurfaceCornersLat,  4, MPI_DOUBLE, 0, MPI_COMM_WORLD);      
    MPI_Bcast (theSurfaceCornersLong, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
}
/**
 * Broadcast green parameters.
 *
 * \note It initializes various global variables (too many to list here).
 */
static void
broadcast_green_parameters (void)
{
    int32_t iSource;

    MPI_Bcast(&theNumberOfPointSources,  1, MPI_INT , 0, MPI_COMM_WORLD);
      
    if ( myID != 0) {
	theSourceLonArray       = malloc( sizeof( double) * theNumberOfPointSources );
	theSourceLatArray       = malloc( sizeof( double) * theNumberOfPointSources );
	theSourceDepthArray     = malloc( sizeof( double) * theNumberOfPointSources );
	theSourceWXArray        = malloc( sizeof( double) * theNumberOfPointSources );
	theSourceWYArray        = malloc( sizeof( double) * theNumberOfPointSources );
	theSourceWZArray        = malloc( sizeof( double) * theNumberOfPointSources );
	theSourceNt1Array       = malloc( sizeof( int )   * theNumberOfPointSources );
	theSourceTinitArray     = malloc( sizeof( double) * theNumberOfPointSources );
	theSourceDtArray        = malloc( sizeof( double) * theNumberOfPointSources );
	theSourceSourceFunArray = malloc( sizeof( double* ) * theNumberOfPointSources );
	if ( (theSourceLonArray    == NULL) || (theSourceLatArray       == NULL) || 
	     (theSourceDepthArray  == NULL) || (theSourceWXArray        == NULL) ||
	     (theSourceWYArray     == NULL) || (theSourceWZArray        == NULL) ||
	     (theSourceNt1Array    == NULL) || (theSourceTinitArray     == NULL) ||
	     (theSourceDtArray     == NULL) || (theSourceSourceFunArray == NULL) ){
	    ABORT_PROGRAM ("Error source green matrices:broadcast_green_parameters ");
	}      
    
    }

    WAIT;
  
    MPI_Bcast (theSourceLonArray,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);	     
    MPI_Bcast (theSourceLatArray,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);
    MPI_Bcast (theSourceDepthArray,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);  
    MPI_Bcast (theSourceWXArray ,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);
    MPI_Bcast (theSourceWYArray ,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);      
    MPI_Bcast (theSourceWZArray ,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);
    MPI_Bcast (theSourceNt1Array,theNumberOfPointSources , MPI_INT, 0,
	       MPI_COMM_WORLD);
    MPI_Bcast (theSourceTinitArray,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);
    MPI_Bcast (theSourceDtArray,theNumberOfPointSources , MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);

    if(myID == 0)
	fprintf(stdout,"\n  Broadcast small arrays source srfh done %d \n", myID);
  
    fflush(stdout);

    WAIT;

    if(myID != 0){
	for ( iSource = 0; iSource < theNumberOfPointSources; iSource++ ){ 
	    /* double check*/
	    if( theSourceNt1Array[iSource] == 0){
		fprintf(stdout,"\n error the size of the slipfunc = 0 in  source%d \n",iSource); 
		fflush(stdout);
		ABORTEXIT;
	    }
	    theSourceSourceFunArray [ iSource ] =
		malloc( sizeof( double ) * theSourceNt1Array[iSource] );      
	}
    }

    WAIT;
  
    for ( iSource = 0; iSource < theNumberOfPointSources; iSource++ ){     
	MPI_Bcast(theSourceSourceFunArray[iSource],
		  theSourceNt1Array[iSource], MPI_DOUBLE, 0,
		  MPI_COMM_WORLD);    
    }
  
    WAIT;
  
    MPI_Bcast (theSurfaceCornersLat,  4, MPI_DOUBLE, 0, MPI_COMM_WORLD);      
    MPI_Bcast (theSurfaceCornersLong, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
}

/**
 * Pass (broadcast) source information to all PEs.
 */
static void broadcast_info (void)
{
    double  dbuf1[10], dbuf2[5], dbuf3[6];
    int32_t ibuf1[3];

    ibuf1[0] = theNumberOfTimeSteps;
    ibuf1[1] = theSourceIsFiltered;
    ibuf1[2] = theTypeOfSource;

    dbuf1[0] = theRegionOriginLatDeg;
    dbuf1[1] = theRegionOriginLongDeg;
    dbuf1[2] = theRegionAzimuthDeg;
    dbuf1[3] = theRegionDepthShallowM;
    dbuf1[4] = theRegionLengthEastM;
    dbuf1[5] = theRegionLengthNorthM;
    dbuf1[6] = theRegionDepthDeepM;
    dbuf1[7] = theDeltaT;
    dbuf1[8] = theValidFreq;

    dbuf2[0] = theAverageRisetimeSec;
    dbuf2[1] = theMomentMagnitude;
    dbuf2[2] = theMomentAmplitude;
    dbuf2[3] = theRickerTs;
    dbuf2[4] = theRickerTp;

    /* broadcast theRegion* vars (9)*/
    MPI_Bcast (dbuf1,		       9, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* broadcast theNumberOfTimeStpes* vars (3) */
    MPI_Bcast (ibuf1,		        3, MPI_INT,    0, MPI_COMM_WORLD);
    MPI_Bcast (&theThresholdFrequency,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);  
    MPI_Bcast (&theNumberOfPoles,       1, MPI_INT,    0, MPI_COMM_WORLD);  
    MPI_Bcast (&theSourceFunctionType,  1, MPI_INT   , 0, MPI_COMM_WORLD);    

    /* broadcast theAverageRisetimeSec, theMomentMagnitude ... (5 total) */
    MPI_Bcast (dbuf2,		      5, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&theNumberOfTimeWindows, 1, MPI_INT   , 0, MPI_COMM_WORLD);

    theRegionOriginLatDeg  = dbuf1[0];
    theRegionOriginLongDeg = dbuf1[1];
    theRegionAzimuthDeg	 = dbuf1[2];
    theRegionDepthShallowM = dbuf1[3];
    theRegionLengthEastM   = dbuf1[4];
    theRegionLengthNorthM  = dbuf1[5];
    theRegionDepthDeepM    = dbuf1[6];
    theDeltaT              = dbuf1[7];
    theValidFreq           = dbuf1[8]; 

    theNumberOfTimeSteps   = ibuf1[0];
    theSourceIsFiltered    = ibuf1[1];
    theTypeOfSource        = ibuf1[2];

    theAverageRisetimeSec  = dbuf2[0];
    theMomentMagnitude     = dbuf2[1];
    theMomentAmplitude     = dbuf2[2];
    theRickerTs		   = dbuf2[3];
    theRickerTp		   = dbuf2[4];

    if ( myID != 0 ) theWindowDelay = dvector(0, theNumberOfTimeWindows);
  
    MPI_Bcast (theWindowDelay, theNumberOfTimeWindows, MPI_DOUBLE,
	       0, MPI_COMM_WORLD);
  
    if (theTypeOfSource == POINT) {
      
	MPI_Bcast (&theLonlatOrCartesian, 1, MPI_INT, 0, MPI_COMM_WORLD);

	dbuf3[0] = theHypocenterLatDeg;
	dbuf3[1] = theHypocenterLongDeg;
	dbuf3[2] = theHypocenterDepthM;
	dbuf3[3] = theSourceStrikeDeg;
	dbuf3[4] = theSourceDipDeg;
	dbuf3[5] = theSourceRakeDeg;	
    
	/* broadcast the vars above (6) */
	MPI_Bcast (dbuf3, 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
	theHypocenterLatDeg  = dbuf3[0];
	theHypocenterLongDeg = dbuf3[1];
	theHypocenterDepthM  = dbuf3[2];
	theSourceStrikeDeg   = dbuf3[3];
	theSourceDipDeg      = dbuf3[4];
	theSourceRakeDeg     = dbuf3[5];

	if( theLonlatOrCartesian == 0 ){
	    MPI_Bcast (theSurfaceCornersLat,  4, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
	    MPI_Bcast (theSurfaceCornersLong, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);  
	}
    }

    if (theTypeOfSource == PLANE  || theTypeOfSource == PLANEWITHKINKS) {
	broadcast_plane_parameters();
    }
  
    if (theTypeOfSource == PLANEWITHKINKS) {
	MPI_Bcast(&theNumberOfKinks,  1, MPI_INT , 0, MPI_COMM_WORLD);    
	if ( myID != 0) {
	    theKinkLonArray = malloc( sizeof(vector3D_t) * theNumberOfKinks);
	    theKinkLatArray = malloc( sizeof(vector3D_t) * theNumberOfKinks);	  
	    if ( (theKinkLonArray == NULL) || (theKinkLatArray == NULL) ) {
		ABORT_PROGRAM (" Error allocating kink matrices");
	    }
	}            
	MPI_Bcast(theKinkLatArray,theNumberOfKinks,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(theKinkLonArray,theNumberOfKinks,MPI_DOUBLE,0,
		  MPI_COMM_WORLD);
      
	MPI_Bcast (theSurfaceCornersLat,  4, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
	MPI_Bcast (theSurfaceCornersLong, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    if(myID == 0)
	fprintf(stdout,"\n  Broadcast general source done %d ", myID);
  
    fflush(stdout);
  
    WAIT;

    if (theTypeOfSource == SRFH) {      
	broadcast_srfh_parameters ( );	
    }
  
    if (theTypeOfSource == GREEN) {      
	broadcast_green_parameters ( );	
    }


    /*fprintf(stdout,"\n  Broadcast done %d \n", myID);*/

    WAIT;
    
    return;

}


/*
 *  init_source_parameters: Read and initialize various static global 
 *                          variables. 
 *
 *               
 *                   input: physicsin contains the domain dimensions and some
 *                          other physical quantities and the path for the dir
 *                          where the files that describe the fault are. 
 *
 *                  return: 1 if OK, -1 on error.
 */
static int init_source_parameters(const char *physicsin){

    FILE *fp, *fparea, *fpstrike, *fpdip, *fprake, *fpslip,  *fpcoords, *fpslipfun;
    FILE *fpweights, *fpsourcefun ;

    char source_directory[256], sourcein[256], slipin[256], slipfunin[256];
    char coordsin[256], areain[256], strikein[256], dipin[256], rakein[256];
    char weightsxyz[256], sourcefunin[256];
    char type_of_source[32];
    
  
    if ( ( fp = fopen(physicsin, "r" ) ) == NULL) {  /*Read Domain and source path*/   
	fprintf(stderr, "Err opening %s\n", physicsin);    
	fflush(stdout);
	return -1;
    }
    if ( read_domain(fp) == -1 ){  /*Read Domain */   
	fprintf(stderr, "Err reading domain coords  %s\n", physicsin);    
	fflush(stdout);
	return -1;
    }


    if ( parsetext( fp, "source_directory", 's', &source_directory ) ){
	fprintf(stderr, "Err parsing src dir from %s",physicsin);
	fflush(stdout);
	return -1;
    }
    fclose(fp);  

    sprintf(sourcein, "%s/input.in" , source_directory); /*Read source properties*/ 
    fp = fopen( sourcein, "r" );     
    if ( ( fp = fopen( sourcein, "r" ) ) == NULL ) {    
	fprintf(stderr, "Error opening %s\n", sourcein);    
	return -1;    
    }  

    if ( read_filter(fp) == -1 ) return -1; /*Read filter*/
  
    /*Formats*/
    if ( parsetext(fp, "type_of_source", 's', &type_of_source) != 0 ) {
	fprintf( stderr, "Err from source.in for type_of_source\n" );
	return -1;
    }

    if (      strcasecmp(type_of_source, "point")   == 0 )  
	theTypeOfSource = POINT;      
    else if ( strcasecmp(type_of_source, "plane")   == 0 )           
	theTypeOfSource = PLANE;          
    else if ( strcasecmp(type_of_source, "planewithkinks") == 0)           
	theTypeOfSource = PLANEWITHKINKS; 
    else if ( strcasecmp(type_of_source, "srfh") == 0)           
	theTypeOfSource = SRFH;   
    else if ( strcasecmp(type_of_source, "green") == 0)           
	theTypeOfSource = GREEN;      
    else {        
	fprintf(stderr, "Err Unknown type_of_source %s\n", type_of_source);
	return -1;
    }
  
    if(myID==0){
	fprintf(stdout,"\nType of source:%s \n",type_of_source);    
    }

    if( read_common_all_formats(fp) == -1 ) return -1;
  
    fprintf(stdout,"\nCommon source parameters read\n");
    fflush(stdout);  
    /*  WAIT;*/



    if ( theTypeOfSource == POINT ){ 
	if( read_point_source(fp) == -1 )
	    return -1;
	fprintf(stdout,"\nSource Point parameters read");
	fflush(stdout);  
    }

    if ( theTypeOfSource == PLANE || theTypeOfSource == PLANEWITHKINKS ){     
  
	sprintf ( slipin,   "%s/slip.in" , source_directory );
	sprintf ( rakein,   "%s/rake.in" , source_directory );
  
	if ((fpslip = fopen(slipin, "r")) == NULL) {
	    fprintf(stderr, "Error opening %s\n", slipin);    
	    return -1;    
	}
     
	if ( ( fprake = fopen(rakein, "r" ) ) == NULL) {    
	    fprintf(stderr, "Error opening %s\n", rakein);    
	    return -1;    
	}
     
	if( read_plane_source(fp, fpslip, fprake) == -1 ) return -1 ;
   
	fprintf(stdout,"\nSource Plane parameters read");
	fflush(stdout);  

	fclose(fpslip);
	fclose(fprake);
    }
  
    if ( theTypeOfSource == SRFH ) {
	sprintf ( coordsin, "%s/coords.in", source_directory );
	sprintf ( areain,   "%s/area.in",   source_directory );
	sprintf ( strikein, "%s/strike.in", source_directory );
	sprintf ( dipin,    "%s/dip.in",  source_directory );
	sprintf ( rakein,   "%s/rake.in", source_directory );
	sprintf ( slipin,   "%s/slip.in", source_directory );
	sprintf ( slipfunin,"%s/slipfunction.in", source_directory );   
	if ( (fpcoords = fopen( coordsin, "r")) == NULL ||  
	     (fparea   = fopen( areain,   "r")) == NULL ||  
	     (fpstrike = fopen( strikein, "r")) == NULL ||  
	     (fpdip    = fopen( dipin,    "r")) == NULL ||  
	     (fprake   = fopen( rakein,   "r")) == NULL ||  
	     (fpslip   = fopen( slipin,   "r")) == NULL ||
	     (fpslipfun   = fopen( slipfunin,   "r")) == NULL  ) {
	    fprintf(stderr, "Error opening srfh files\n" );    
	    return -1;    
	}    
	if( read_srfh_source(fp,fpcoords, fparea, fpstrike, fpdip,
			     fprake, fpslip, fpslipfun) == -1 )
	    return -1;           
	fclose(fpcoords);
	fclose(fparea);
	fclose(fpstrike);
	fclose(fpdip);
	fclose(fprake);
	fclose(fpslip);
	fclose(fpslipfun);
    }
  
    if ( theTypeOfSource == GREEN ) {
	fprintf(stdout,"\nstart reading green");

	sprintf ( coordsin,   "%s/coords.in", source_directory );
	sprintf ( weightsxyz, "%s/weightsxyz.in",   source_directory );

	sprintf ( sourcefunin,"%s/sourcefunction.in", source_directory );   

	if ( (fpcoords     = fopen( coordsin,    "r")) == NULL ) {
	    fprintf(stderr,"Error opening coords.in green files\n" );    
	    return -1;    
	}  
	if ( (fpweights    = fopen( weightsxyz,  "r")) == NULL) {
	  fprintf(stderr,"Error opening %s green files\n",weightsxyz );    
	    return -1;    
	}  
	if ((fpsourcefun  = fopen( sourcefunin,  "r")) == NULL ) {
	    fprintf(stderr,"Error opening sourcefunction.in green files\n" );    
	    return -1;    
	}  

  
	if( read_green_source(fp,fpcoords, fpweights, fpsourcefun) == -1 )
	    return -1;           
	fclose(fpcoords);
	fclose(fpweights);
	fclose(fpsourcefun);
    }

    fclose(fp);  
    return 1;
}

/* 
 *  compute_print_source: Compute and print the forces needed in
 *                        psolve
 *    input: 
 *          informationdirectory - where source related files are located
 *                      myoctree - octree to be used
 *                        mymesh - information realtive to the mesh
 *                      myforces - the load vector
 *           numericsinformation - see quakesource.h for description
 *                mpiinformation - see quakesource.h for description 
 *
 *    return:     
 *  
 */
int compute_print_source( const char *physicsin, octree_t *myoctree, 
			  mesh_t *mymesh, numerics_info_t numericsinformation,
			  mpi_info_t mpiinformation ){    
    /*Mesh related */
    myOctree = myoctree; 
    myMesh   = mymesh;

    /*MPI Related*/
    myID         = mpiinformation.myid;  
    theGroupSize = mpiinformation.groupsize;

    /*Numerics related*/
    theDeltaT            = numericsinformation.deltat;
    theNumberOfTimeSteps = numericsinformation.numberoftimesteps;
    theValidFreq         = numericsinformation.validfrequency;
    theMinimumEdge       = numericsinformation.minimumh;

    theDomainX = numericsinformation.xlength;
    theDomainY = numericsinformation.ylength;
    theDomainZ = numericsinformation.zlength;
  
    myForces = calloc( myMesh->nharbored, sizeof(vector3D_t *) );  
    if ( myForces == NULL ) {    
	fprintf(stderr, "Thread %d: quakesource: out of memory\n",myID);    
	ABORTEXIT;
    }
  

    if ( myID == 0 ) 
	if ( init_source_parameters( physicsin ) == -1 ){      
	    fprintf(stdout,"Err init_source_parameters failed");
	    ABORTEXIT;
	}
  
  
    broadcast_info();
  
    /*  fprintf(stdout,"\nBroadcast done outside id %d", myID);*/
    /*  fflush(stdout); */


    if ( theTypeOfSource == POINT )
	if ( compute_myForces_point(physicsin) == -1 ){
	    fprintf(stdout,"Err compute_myForces_point failed");
	    ABORTEXIT;
	}
   
    
    if (  theTypeOfSource == PLANE || theTypeOfSource == PLANEWITHKINKS )
	if ( compute_myForces_planes(physicsin) == -1 ){       
	    fprintf(stdout,"Err compute_myForces_planes failed");
	    ABORTEXIT;
	}
  

    /* Standard Rupture Fault Hercules, variation of SRF by Graves */
    if(  theTypeOfSource == SRFH )
	if( compute_myForces_srfh(physicsin )==-1){
	    fprintf(stdout,"Err compute_myForces_point failed");
	    ABORTEXIT;
	}  
  
    if(  theTypeOfSource == GREEN ){
      //	fprintf(stdout,"Entro a Green");
	fflush(stdout);
	if( compute_myForces_green(physicsin )==-1){
	    fprintf(stdout,"Err compute_myForces_green failed");
	    ABORTEXIT;
	}  
    }
  
    MPI_Barrier( MPI_COMM_WORLD );
  
    source_compute_print_stat();
  
    MPI_Barrier( MPI_COMM_WORLD );
 
    return 0;

}
