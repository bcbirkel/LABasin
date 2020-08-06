#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int Search_Inpoly( polyg_t *polygon, nodedb_t *nodes, double xt, double yt);

double Phi_Local( int i, double csi,double etha, double *csii, double *ethai);

double Interp_In_Ele_Global( double x, double y, double **zgrid, double *xgrid, double *ygrid, int nxinzgrid, int nyinzgrid );

double Interp_In_Ele_Global_N( double x, double y, double **zgrid, double *xgrid, double *ygrid, int nxinzgrid, int nyinzgrid );

  
static double  Interpolate_Linear ( double x,int numsamples, double samplingdx, 
				    double *discretefunction );

static float  Interpolate_Linear_General(float x, int numsamples,float *xdiscrete,float *discretefunction);
void Update_Stencil_Plane(meshdb_t *current, double y);


double Interp_In_Ele_Global_Vel( double x, double y, int surfid, double **zgrid, double *xgrid, 
				 double *ygrid, int nxinzgrid, int nyinzgrid,
				 int numprofiles, velprofile_t *velprofiles, double *vp,double *vs,
				 double *rho, double depth);
