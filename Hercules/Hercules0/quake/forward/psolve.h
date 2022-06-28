/* -*- C -*-
 *
 * psolve.h: Generate an unstructured mesh and solve the linear system
 *           thereof derived and output the results.
 *
 * Input:    material database (cvm etree), physics.in, numerical.in.
 * Output:   mesh database (mesh.e) and 4D output.
 *
 * Copyright (c) 2005 Tiankai Tu, Hongfeng Yu, Leonardo Ramirez-Guzman
 *
 * All rights reserved.  May not be used, modified, or copied
 * without permission.
 *
 * Contact:
 * Tiankai Tu
 * Computer Science Department
 * Carnegie Mellon University
 * 5000 Forbes Avenue
 * Pittsburgh, PA 15213
 * tutk@cs.cmu.edu
 *
 * Notes:
 *   For a history of changes see the ChangeLog file.  Also, if you
 *   perform any changes to the source code, please document the changes
 *   by adding the appropriate entry in the ChangeLog file and a descriptive
 *   message during CVS commit.  Thanks!
 */

#ifndef PSOLVE_H
#define PSOLVE_H

#include <inttypes.h>
#include <stdio.h>


/* physics related */

typedef enum
{

  RAYLEIGH = 0, MASS, NONE

} damping_type_t;



typedef struct corner_ref_t
{
    int32_t corner_lnid[8];
    double min_value, max_value;

} corner_ref_t;



/**
 * VIS_GETVELOCITY: Get the velocity at a corner (cornerid) of an octant in
 *               a particular direction (axis = 0, 1, 2 for x, y, z).
 *
 * - Assume octant->appdata is not NULL. Otherwise, the macro causes crash.
 *
 */
#define VIS_LNID(octant,cornerid) (((corner_ref_t *)octant->appdata)->corner_lnid[cornerid])
#define VIS_GETMINVALUE(octant) (((corner_ref_t *)octant->appdata)->min_value)
#define VIS_GETMAXVALUE(octant) (((corner_ref_t *)octant->appdata)->max_value)

extern double theGlobalDeltaT; /* delta time (sec) per timestep */

#define VIS_GETVELOCITY(solver,octant,cornerid,axis) \
((solver->tm2[VIS_LNID(octant,cornerid)].f[axis] - \
  solver->tm1[VIS_LNID(octant,cornerid)].f[axis]) / theGlobalDeltaT)



/**
 * mdata_t: Mesh database (element) record payload.
 *
 */
typedef struct mdata_t
{
    int64_t nid[8];
    float edgesize, Vp, Vs, rho;

} mdata_t;



/**
 * fvector_t: 3-ary double vector
 *
 */
typedef struct fvector_t
{
    double f[3];

} fvector_t;



/**
 * out_hdr_t: 4D output file header.
 */
typedef struct out_hdr_t
{
    /** file type string identifier: "Hercules 4D output vnnn"  */
    char    file_type_str[29];
    int8_t  format_version;	/**< File format version		*/
    int8_t  endiannes;		/**< File endianess: 0=little, 1=big	*/

    /** Identifier of the platform where the file was generated */
    int8_t   platform_id;
    unsigned char ufid[16];	/**< "Unique" file identifier		*/
    int64_t  total_nodes;	/**< Node count.			*/
    int32_t  output_steps;	/**< Number of output time steps.	*/

    /** Number of components per (node) record, e.g., 1 vs. 3 */
    int32_t  scalar_count;
    int8_t   scalar_size;	/**< size (in bytes) of each scalar value  */

    /**<
     * Type of scalar, it can take one of the following values
     * - INVALID:		 0
     * - FLOAT32 (float):	 1
     * - FLOAT64 (double):	 2
     * - FLOAT128 (long double): 3
     * - INT8:			 4
     * - UINT8:			 5
     * - INT16:			 6
     * - UINT16:		 7
     * - INT32:			 8
     * - UINT32:		 9
     * - INT64:			 10
     * - UINT64:		 11
     */
    int8_t  scalar_type;

    /**<
     * scalar class, it can take one of the following values:
     * - INVALID:		 0
     * - FLOAT_CLASS		 1
     * - INT_CLASS		 2
     */
    int8_t  scalar_class;

    int8_t  quantity_type;	/* 0: unknown, 1: displacement, 2: velocity */


    /** Mesh parameters: extent of the simulated region in meters */
    double  domain_x, domain_y, domain_z;

    /**
     * Mesh parameter: tick size, factor for converting from domain units (m)
     * to etree units.
     */
    double  mesh_ticksize;

    /** Simulation parameter: delta t (time) */
    double  delta_t;

    /** Mesh parameter: total number of elements */
    int64_t total_elements;

    /** 4D Output parameter: how often is an output time step written out */
    int32_t output_rate;

    /** Simulation parameter: total number of simulation time steps */
    int32_t total_time_steps;

    int64_t generation_date;	/**< Time in seconds since the epoch	*/
} out_hdr_t;



/*---------------- Solver data structures -----------------------------*/

/**
 * e_t: Constants initialized element structure.
 *
 */
typedef struct e_t
{
    double c1, c2, c3, c4;

} e_t;



/**
 * n_t: Constants initialized node structure.
 *
 *      This structure has been changed as a result of the new algorithm
 *      for the terms involved with the damping in the solution for the
 *      next time step.  Old structure is commented and we will keep it
 *      like that until we are done with the Terashake runs.
 *
 *      This changed is being administered by RICARDO
 *      May and June 2006
 *
 */
typedef struct n_t
{
/*  Old version of the structure */
/*  double mass2x;               */
/*  double mplus[3], mminus[3];  */

/*  New version                  */
    double mass_simple;
    double mass2_minusaM[3];
    double mass_minusaM[3];

} n_t;



/**
 * fmatrix_t: stiffness matrices
 *
 */
typedef struct fmatrix_t
{
    double f[3][3];

} fmatrix_t;



/**
 * solver_t: Solver abstraction.
 *
 */
typedef struct solver_t
{
    e_t *eTable;                  /* Element computation-invariant table */
    n_t *nTable;                  /* Node computation-invariant table */

    fvector_t *tm1;               /* Displacements at timestep t - 1 */
    fvector_t *tm2;               /* Displacements at timestep t - 2 */
    fvector_t *force;             /* Force accumulation at timestep t */

} solver_t;




#ifdef __cplusplus
extern "C" {
#endif

extern int  solver_abort (const char* function_name, const char* error_msg,
			  const char* format, ...);
extern void solver_output_seq (void);
extern int  parsetext (FILE* fp, const char* querystring, const char type,
		       void* result);



#ifdef __cplusplus
}
#endif
#endif /* PSOLVE_H */
