/*
 * Error.h
 *
 *
 * Copyright (c) 2005 Leonardo Ramirez-Guzman
 * All rights reserved.  May not be used, modified, or copied
 * without permission.
 *
 * Contact:
 * Leonardo Ramirez-Guzman
 * Civil and Environmental Engineering
 * Carnegie Mellon University
 * 5000 Forbes Avenue
 * Pittsburgh, PA 15213
 * lramirez@andrew.cmu.edu
 *
 */


static void Check_Error(int condition, const char *errormessage ){

    if ( condition ){
	
	fprintf(stdout," \nErr: %s\n ", errormessage);
	exit(1);
    }
    
    return;
	
}


