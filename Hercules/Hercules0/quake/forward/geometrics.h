


#define PI 3.14159265358979323846264338327  


typedef struct vector3D_t{
  
  double x [ 3 ];
  
}vector3D_t;




vector3D_t compute_global_coords( vector3D_t origin, vector3D_t local,
				  double dip, double rake ,double strike);


vector3D_t compute_centroid(vector3D_t *p);
  

double compute_area(vector3D_t *p);


int compute_1D_grid( double cellSize, int numberOfCells, int pointsInCell,
		     double minimumEdgeTriangle, double *grid);


vector3D_t compute_domain_coords(vector3D_t point, double azimuth);


vector3D_t compute_domain_coords_linearinterp(double lon , double lat, 
					      double *longcorner, 
					      double *latcorner,
					      double domainlengthetha,
					      double domainlengthcsi);

int compute_lonlat_from_domain_coords_linearinterp(double y , double x , 
						   double *lon,double *lat,
						   double *longcorner, 
						   double *latcorner,
						   double domainlengthX,
						   double domainlengthY);

