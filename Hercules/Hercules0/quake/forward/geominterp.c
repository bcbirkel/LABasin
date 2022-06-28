#include "geodataq.h"
#include "geominterp.h"

/*-------------------------------------------------------------------------*
 *                                                                         *
 *   Search_Inpoly                                                         *
 *                                                                         *
 *-------------------------------------------------------------------------*/
int Search_Inpoly( polyg_t *polygon, nodedb_t *nodes, double xt, double yt){

    double xnew,ynew,xold,yold;
    double x1,y1,x2,y2;
    int i,inside=0;
    int *incid, npoints;
    
    incid = polygon->incidences;
    npoints = polygon -> numnodes;   
    
    if (npoints < 3) 	return(0);
    
    xold=nodes[incid[npoints-1]].x;
    yold=nodes[incid[npoints-1]].y;
    for (i=0 ; i < npoints ; i++) {
	xnew=nodes[incid[i]].x;
	ynew=nodes[incid[i]].y;
	if (xnew > xold) {
	    x1=xold;
	    x2=xnew;
	    y1=yold;
	    y2=ynew;
	}
	else {
	    x1=xnew;
	    x2=xold;
	    y1=ynew;
	    y2=yold;
	}
	if ((xnew < xt) == (xt <= xold)         /* edge "open" at left end */
	    && (yt-y1)*(x2-x1)
            < (y2-y1)*(xt-x1)) {
	    inside=!inside;
	}
	xold=xnew;
	yold=ynew;
    }
    return(inside);
}


/*-------------------------------------------------------------------------*
 *
 *
 *   Phi_Local: local shape function
 *
 *                 1----------------2
 *                 |       ^etha    |
 *                 |       |   csi  |
 *                 |        ---->   |
 *                 |                | 
 *                 |                | 
 *                 0----------------3
 *
 *   input:       i - shape function number
 *              csi - horizontal coord
 *            ethai - vertical coord
 *      ethai  csii - arrays with the convention of orintation in csi and etha
 *
 *--------------------------------------------------------------------------*/
double Phi_Local( int i, double csi,double etha, double *csii, double *ethai){

    double val;
  
    val = 0.25 * ( 1 + csii[i] * csi) * ( 1 + ethai[i]* etha);

    return val;
}




/*-------------------------------------------------------------------------*
 *
 *
 *   Interp_In_Ele_Global: interpolate in a regular grid using linear 
 *                         elements
 *
 *                  1----------------2
 *                  |                |
 *         Y        |                |
 *         ^        |                |
 *         |        |                | 
 *         |        |                | 
 *         |        0----------------3
 *         ----->X
 *
 *-------------------------------------------------------------------------*/
double Interp_In_Ele_Global( double x, double y, double **zgrid, double *xgrid, double *ygrid,
			     int nxinzgrid, int nyinzgrid ){
  
    int nodeX[4], nodeY[4], iNode;  
        double stepX, stepY,csi,etha, zInterp, phi;

    static double csiI[4]  = {-1, -1,  1,  1};
    static double ethaI[4] = {-1,  1,  1, -1} ;
  
    /* sanitiy check */ 
    if( x < xgrid[0] || x > xgrid[nxinzgrid-1] || y < ygrid[0] || y > ygrid[nyinzgrid-1] ){
	return 0;
    }     
  
    /* step in x and y */
    /*    stepX = xgrid[1]-xgrid[0];
	  stepY = ygrid[1]-ygrid[0];*/   
    stepX=(xgrid[nxinzgrid-1]-xgrid[0])/(nxinzgrid-1);
    stepY=(ygrid[nyinzgrid-1]-ygrid[0])/(nyinzgrid-1);

    /* search in the mesh for the node 1*/
    nodeX[0] = (int)floor( (x-xgrid[0]) / stepX ); 
    nodeY[0] = (int)floor( (y-ygrid[0]) / stepY );         
    if(x == xgrid[nxinzgrid-1]) nodeX[0]= nodeX[0]-1;
    if(y == ygrid[nyinzgrid-1]) nodeY[0]= nodeY[0]-1;
    nodeX[1] = nodeX[0];
    nodeY[1] = nodeY[0]+1;
    nodeX[2] = nodeX[0]+1;
    nodeY[2] = nodeY[0]+1;
    nodeX[3] = nodeX[0]+1;
    nodeY[3] = nodeY[0];

    /* compute csi and etha */  
    csi  = xgrid[nodeX[0]]+.5*stepX;
    etha = ygrid[nodeY[0]]+.5*stepY;  
    csi  = 2 * ( x-csi  ) / stepX;
    etha = 2 * ( y-etha ) / stepY;
  
    /* loop over all nodes and compute the value of z interpolated */
    zInterp =0;
  
  
    for ( iNode = 0; iNode < 4; iNode++){   
	/* to debug */
	if( nodeX[iNode] > nxinzgrid || nodeY[iNode] > nyinzgrid){
	    fprintf(stdout,"x=%lf y=%lf",x,y);
	    fprintf(stdout,"Error nodex=%d> %d  nodey=%d > %d ", nodeX[iNode], nxinzgrid,
		    nodeY[iNode], nyinzgrid  );
	    fflush(stdout);
	    exit(1);
	}
      
   
	phi = Phi_Local( iNode, csi, etha, &(csiI[0]), &(ethaI[0]));    
	zInterp += phi * zgrid[nodeX[iNode]][nodeY[iNode]];
    }
  
    return zInterp;
}

/*-------------------------------------------------------------------------*
 *
 *
 *   Interp_In_Ele_Global_N: interpolate in a regular grid using the nearest point 
 *                         elements
 *
 *                  1----------------2
 *                  |                |
 *         Y        |                |
 *         ^        |                |
 *         |        |                | 
 *         |        |                | 
 *         |        0----------------3
 *         ----->X
 *
 *-------------------------------------------------------------------------*/
double Interp_In_Ele_Global_N( double x, double y, double **zgrid, double *xgrid, double *ygrid,
			     int nxinzgrid, int nyinzgrid ){
  
   int nodeX[4], nodeY[4], iNode;  
        double stepX, stepY,csi,etha, zInterp, phi;

    static double csiI[4]  = {-1, -1,  1,  1};
    static double ethaI[4] = {-1,  1,  1, -1} ;
  
    /* sanitiy check */ 
    if( x < xgrid[0] || x > xgrid[nxinzgrid-1] || y < ygrid[0] || y > ygrid[nyinzgrid-1] ){
	return 0;
    }     
  
    /* step in x and y */
    /*    stepX = xgrid[1]-xgrid[0];
	  stepY = ygrid[1]-ygrid[0];*/   
    stepX=(xgrid[nxinzgrid-1]-xgrid[0])/(nxinzgrid-1);
    stepY=(ygrid[nyinzgrid-1]-ygrid[0])/(nyinzgrid-1);

    /* search in the mesh for the node 1*/
    nodeX[0] = (int)floor( (x-xgrid[0]) / stepX ); 
    nodeY[0] = (int)floor( (y-ygrid[0]) / stepY );         
    if(x == xgrid[nxinzgrid-1]) nodeX[0]= nodeX[0]-1;
    if(y == ygrid[nyinzgrid-1]) nodeY[0]= nodeY[0]-1;
    nodeX[1] = nodeX[0];
    nodeY[1] = nodeY[0]+1;
    nodeX[2] = nodeX[0]+1;
    nodeY[2] = nodeY[0]+1;
    nodeX[3] = nodeX[0]+1;
    nodeY[3] = nodeY[0];

    /* compute csi and etha */  
    csi  = xgrid[nodeX[0]]+.5*stepX;
    etha = ygrid[nodeY[0]]+.5*stepY;  
    csi  = 2 * ( x-csi  ) / stepX;
    etha = 2 * ( y-etha ) / stepY;
  
    /* loop over all nodes and compute the value of z interpolated */
    zInterp =0;
  
    if(csi!=0)csi=csi/fabs(csi);
    if(etha!=0)etha=etha/fabs(etha);
  
    for ( iNode = 0; iNode < 4; iNode++){   
	/* to debug */
	if( nodeX[iNode] > nxinzgrid || nodeY[iNode] > nyinzgrid){
	    fprintf(stdout,"x=%lf y=%lf",x,y);
	    fprintf(stdout,"Error nodex=%d> %d  nodey=%d > %d ", nodeX[iNode], nxinzgrid,
		    nodeY[iNode], nyinzgrid  );
	    fflush(stdout);
	    exit(1);
	}
      
   
	phi = Phi_Local( iNode, csi, etha, &(csiI[0]), &(ethaI[0]));    
	zInterp += phi * zgrid[nodeX[iNode]][nodeY[iNode]];
    }
  
    return zInterp;
}

/*-------------------------------------------------------------------------*
 *
 *
 * interpolate_linear: interpolate linearly a function 
 *                     if the time is larger than the one supported by the function
 *                     the last value will be assigned
 *
 *
 *-------------------------------------------------------------------------*/

static double  Interpolate_Linear ( double x,int numsamples, double samplingdx, 
				    double *discretefunction ){
    double maxx, m, b;

    int interval;

    maxx = (numsamples-1)*samplingdx;

    if ( x >= maxx ){    
	return discretefunction[numsamples-1];
    }
    else{

	/* locate the interval */
	interval = floor(x/samplingdx);
    
	/* y = mx +b */
	m =  discretefunction[interval+1] - discretefunction[interval] ;
	m =  m/samplingdx;

	b =  discretefunction[interval] - m * interval * samplingdx;

	return m*x+b;    
    }    
}





/*-------------------------------------------------------------------------*
 *
 *
 * interpolate_linear_general: naive implementation to interpolate linearly a function 
 *                             if the x is larger or smaller than the one 
 *                             supported by the function
 *                             the last or first value will be assigned
 *
 *
 *-------------------------------------------------------------------------*/

static float  Interpolate_Linear_General(float x, int numsamples,float *xdiscrete,
					 float *discretefunction){
  float maxx,minx, m, b,samplingDx;
  int iInterval;

    int interval;
    minx=  xdiscrete[0];
    maxx = xdiscrete[numsamples-1];

    if ( x >= maxx )
	return discretefunction[numsamples-1];
    if ( x <= minx)
	return discretefunction[0]; 

    else{
       
	/* locate the interval using a n complexity this has to change */

        /* go interval by interval */
	for (iInterval=0; iInterval<numsamples; iInterval++){
	  if( xdiscrete[iInterval] <= x )
	    if( x <= xdiscrete[iInterval+1]){
		interval=iInterval;
		samplingDx=xdiscrete[iInterval+1]-xdiscrete[iInterval];	    
	    }
	}	
    
	/* y = mx +b */
	m =  discretefunction[interval+1] - discretefunction[interval] ;
	m =  m/samplingDx;

	b =  discretefunction[interval] - m * xdiscrete[interval];

	return m*x+b;    
    }    
}







/*-------------------------------------------------------------------------*
 *
 *
 * Update_Stencil_Plane: when using a mesh based on polynomials the each node
 *                       can havea a variation in the normal direccion of the
 *                       2d mesh.
 *
 *-------------------------------------------------------------------------*/

void Update_Stencil_Plane(meshdb_t *current, double y){

    int iNode, iNodeInterp; 
    
    for ( iNode = 0; iNode < current->numnodesinterpolated; iNode++){

	iNodeInterp = current->nodesinterpolated[iNode];

	current->nodes[iNodeInterp].x=  
	    Interpolate_Linear ( y, current->howmanynodes[iNode],
				 current->nodesinterpcoordsy[iNode][1]-
				 current->nodesinterpcoordsy[iNode][0],
				 current->nodesinterpcoordsx[iNode]);
    }
    return;
}



double Interp_In_Ele_Global_Vel( double x, double y, int surfid, double **zgrid, double *xgrid, 
				 double *ygrid, int nxinzgrid, int nyinzgrid,
				 int numprofiles, velprofile_t *velprofiles, double *vp,double *vs,
				 double *rho, double depth){
    
    int nodeX[4], nodeY[4], iNode, iProfile;  
    double stepX, stepY,csi,etha, zInterp, phi;
    double velP, velS, density;
    static double csiI[4]  = {-1, -1,  1,  1};
    static double ethaI[4] = {-1,  1,  1, -1} ;
  
    /* sanitiy check */ 
    if( x < xgrid[0] || x > xgrid[nxinzgrid-1] ||
	y < ygrid[0] || y > ygrid[nyinzgrid-1] )
	return 0;
        
    /* step in x and y */
    /*    stepX = xgrid[1]-xgrid[0];
	  stepY = ygrid[1]-ygrid[0];*/
    stepX = (xgrid[nxinzgrid-1]-xgrid[0])/((double)(nxinzgrid-1));
    stepY = (ygrid[nyinzgrid-1]-ygrid[0])/((double)(nyinzgrid-1));
  

  
    /* search in the mesh for the node 1*/
    nodeX[0] = (int)floor( (x-xgrid[0]) / stepX ); 
    nodeY[0] = (int)floor( (y-ygrid[0]) / stepY );   
    
    if(x == xgrid[nxinzgrid-1]) nodeX[0]= nodeX[0]-1;
    if(y == ygrid[nyinzgrid-1]) nodeY[0]= nodeY[0]-1;
    
    nodeX[1] = nodeX[0];
    nodeY[1] = nodeY[0]+1;
    nodeX[2] = nodeX[0]+1;
    nodeY[2] = nodeY[0]+1;
    nodeX[3] = nodeX[0]+1;
    nodeY[3] = nodeY[0];
    
    /* compute csi and etha */  
    csi  = xgrid[nodeX[0]]+.5*stepX;
    etha = ygrid[nodeY[0]]+.5*stepY;  
    csi  = 2 * ( x-csi  ) / stepX;
    etha = 2 * ( y-etha ) / stepY;
    
    /* loop over all nodes and compute the value of z interpolated */
    zInterp=0; 
    *vs    =0;
    *vp    =0;
    *rho   =0;
 
    for ( iNode = 0; iNode < 4; iNode++){   
	/* to debug */
	if( nodeX[iNode] > nxinzgrid || nodeY[iNode] > nyinzgrid){
	    fprintf(stdout,"x=%lf y=%lf",x,y);
	    fprintf(stdout,"Error nodex=%d> %d  nodey=%d > %d ", 
		    nodeX[iNode], nxinzgrid, nodeY[iNode], nyinzgrid  );
	    fflush(stdout);
	    exit(1);
	}
      
	iProfile=(int)(zgrid[nodeX[iNode]][nodeY[iNode]]);
	
	if(iProfile >= numprofiles){
	    fprintf(stdout,
		    "\nError input database iProfile=%d > numProfiles=%d in Surf=%d", 
                     iProfile,numprofiles,surfid);
	    fflush(stdout);
	    exit(1);
	}
      


	velP    = (double)Interpolate_Linear_General(depth, velprofiles[iProfile].numpoints,
				     velprofiles[iProfile].depth,velprofiles[iProfile].vp );
	velS    = (double)Interpolate_Linear_General(depth, velprofiles[iProfile].numpoints,
				     velprofiles[iProfile].depth,velprofiles[iProfile].vs );
	density = (double)Interpolate_Linear_General(depth, velprofiles[iProfile].numpoints,
				     velprofiles[iProfile].depth,velprofiles[iProfile].rho );
	
	
	phi = Phi_Local( iNode, csi, etha, &(csiI[0]), &(ethaI[0]));    	
	*vs     += phi * velS;
	*vp     += phi * velP;
	*rho    += phi * density;
    }  
    return zInterp;
}


