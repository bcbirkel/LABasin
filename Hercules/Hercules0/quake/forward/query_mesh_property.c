/**
 * q4.c - Query the time series of an observation point.
 *
 * Copyright (c) 2010 Leonardo Ramirez-Guzman
 * All rights reserved.  May not be used, modified, or copied 
 * without permission.
 *
 * Leonardo Ramirez-Guzman
 *
 */

#include "query_mesh_property.h"


int mesh_point(double x, double y, double z, etree_t *mep,out_hdr_t result_hdr,
               double *vp,double *vs, double *rho)
{
    double edgesize, ldb[3], center[3], distance[3], phi[8];
    mdata_t mdata;
    etree_addr_t searchAddr, elemAddr;
    int out_step, which;
    off_t offset_base;

    searchAddr.x = (etree_tick_t)(x / result_hdr.mesh_ticksize);
    searchAddr.y = (etree_tick_t)(y / result_hdr.mesh_ticksize);
    searchAddr.z = (etree_tick_t)(z / result_hdr.mesh_ticksize);
    searchAddr.level = ETREE_MAXLEVEL;

    /* Search the mesh etree for the element that contains the point*/
    if (etree_search(mep, searchAddr, &elemAddr, NULL, &mdata) != 0) {
        fprintf(stderr, "%s\n", etree_strerror(etree_errno(mep)));
        return -1;
    }

    *vp=mdata.Vp;
    *vs=mdata.Vs;
    *rho=mdata.rho;
    
    return 0;
    
}
    
