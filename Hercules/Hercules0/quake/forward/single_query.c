/**
 * single_query.c - Query the time series of a single observation point.
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

#include <stdio.h>
#include <stdlib.h>

#include "q4.h"


int main(int argc, char **argv)
{
    etree_t *mep;
    FILE *result_fp;
    out_hdr_t result_hdr;
    double x, y, z;

    if (argc != 6) {
        fprintf(stderr, "usage: single_query meshetree 4d-output x y z\n");
        fprintf(stderr, "meshetree: path to the output mesh etree.\n");
        fprintf(stderr, "4D-output: path to the 4D simulation results.\n");
        fprintf(stderr, 
                "x, y, z: coordinate (in meters) of the query point\n");
        fprintf(stderr, "\n"); 

        return -1;
    }

    if ((sscanf(argv[3], "%lf", &x) != 1) ||
        (sscanf(argv[4], "%lf", &y) != 1) ||
        (sscanf(argv[5], "%lf", &z) != 1)) {
        perror("sscanf");
        return -1;
    }

    /* Open mesh databases */
    if ((mep = etree_open(argv[1], O_RDONLY, 0, 0, 0)) == NULL) {
        fprintf(stderr, "Cannot open mesh element etree %s\n", argv[1]);
        exit(1);
    }

    /* Open 4D output */
    if ((result_fp = fopen(argv[2], "r")) == NULL) {
        fprintf(stderr, "Cannot open 4D output file %s\n", argv[2]);
        exit(1);
    }

    /* Get the metadata */
    if (fread(&result_hdr, sizeof(out_hdr_t), 1, result_fp) != 1) {
        fprintf(stderr, "Cannot read 4D-out result header info from %s.\n", 
                argv[2]);
        exit(1);
    }

    if (q4_point(x, y, z, mep, result_fp, result_hdr, stdout, 1) != 0) {
        fprintf(stderr, "Cannot find the query point %f %f %f\n", 
                x, y, z);
        exit(1);
    }

    fclose(result_fp);
    etree_close(mep);
    
    return 0;
}

