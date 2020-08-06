/**
 * q4.h - Query the time series of an observation point.
 *
 * Copyright (c) 2005 Tiankai Tu 
 * All rights reserved.  May not be used, modified, or copied 
 * without permission.
 *
 * Tiankai Tu
 * Computer Science Department
 * Carnegie Mellon University
 * 5000 Forbes Avenue
 * Pittsburgh, PA 15213
 * tutk@cs.cmu.edu
 *
 */

#ifndef Q4_H
#define Q4_H

#include <stdio.h>
#include "psolve.h"
#include "etree.h"


/**
 * q4_point: query a point (x, y, z) from the dataset of mesh (mep) and
 *           simulation result file (result_fp). The meta data of the 
 *           dataset should be passed in with result_hdr. The output
 *           is printed in out_format (0 for binary, 1 for ASCII) into outfp.
 *
 *           Return 0 if successful; -1 if the point couldn't be found.
 */
extern int
q4_point(double x, double y, double z, etree_t *mep, FILE *result_fp,
         out_hdr_t result_hdr, FILE *outfp, int out_format);

#endif /* Q4_H */
