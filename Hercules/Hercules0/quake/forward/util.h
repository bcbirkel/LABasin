/* -*- C -*-
 *
 * Copyright (C) 2007 Julio Lopez. All Rights Reserved.
 * For copyright and disclaimer details see the file COPYING
 * included in the package distribution.
 *
 * Package:     CMU Quake project's Hercules numerical solver.
 * Description: Miscellaneous support functions.
 *
 * Author:      Julio Lopez (jclopez at ece dot cmu dot edu)
 *
 * File info:   $RCSfile: util.h,v $
 *              $Revision: 1.7 $
 *              $Date: 2007/07/16 16:42:59 $
 *              $Author: jclopez $
 */
#ifndef HERC_UTIL_H
#define HERC_UTIL_H

#include <sys/types.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#include "funcname.h"

/** Replace fseeko with fseek */
#ifdef BIGBEN
#define fseeko          fseek
#endif


/**
 * Macros with the format string for int64_t and uint64_t types.
 */
#if (defined __WORDSIZE) && (__WORDSIZE == 64)
#  define UINT64_FMT		"lu"
#  define INT64_FMT		"ld"
#  define MPI_INT64		MPI_LONG
#else /*  __WORDSIZE && __WORDSIZE == 64 */
#  define UINT64_FMT		"llu"
#  define INT64_FMT		"lld"
#  define MPI_INT64		MPI_LONG_LONG_INT
#endif

#if (defined _FILE_OFFSET_BITS) && (_FILE_OFFSET_BITS == 64)
#  define OFFT_FMT		UINT64_FMT
#else
#  define OFFT_FMT		"d"
#endif



#ifndef STRINGIZE
# ifdef HAVE_STRINGIZE
#  define STRINGIZE(x)          #x
# else
#define STRINGIZE(s)            "x"
# endif
#endif /* STRINGIZE */



#ifndef STRINGX
# ifdef HAVE_STRINGIZE
#   define STRINGX(x)           STRINGIZE(x)
# else /*  HAVE_STRINGIZE */
#   define STRINGX(x)           "x"
# endif /*  HAVE_STRINGIZE */
#endif /* ! STRINGX */


#define HU_VOID_CAST_EXP	(void)

/**
 * Execute a given statement only if the condition is true.
 * The goal of this macro is to simplify certain code constructions such as:
 *
 * if (foo) {
 *	bar();
 * }
 *
 * and replace them with
 *   SUL_COND_EXEC (foo, bar());
 *
 * It is mainly intended to be used inside other macros.
 */
#define HU_COND_EXEC(flag,statement)  \
    (HU_VOID_CAST_EXP ((flag) ? (statement) : 0 ))

 

#define HU_ASSERT_COND(p,c)						\
    do {                                                                \
        if ((c) && ! (p)) {                                             \
            hu_assert_abort( __FILE__,	__FUNCTION_NAME, __LINE__,	\
			     STRINGX(p) );				\
        }                                                               \
    } while (0)



#ifdef HERCULES_DISABLE_ASSERT

#define HU_ENABLE_ASSERT        0
#define HU_ASSERT(p)
#define HU_ASSERT_PTR(p)

#else  /* HERCULES_DISABLE_ASSERT */

#define HU_ENABLE_ASSERT	1
#define HU_ASSERT(p)            HU_ASSERT_COND(p,HU_ENABLE_ASSERT)
#define HU_ASSERT_PTR(p)        HU_ASSERT(NULL != (p))

#endif /* HERCULES_DISABLE_ASSERT */


/**
 * Unconditional assertions, these are not disabled at compile time.
 */
#define HU_ASSERT_ALWAYS(p)     HU_ASSERT_COND(p,1)
#define HU_ASSERT_PTR_ALWAYS(p) HU_ASSERT_ALWAYS(NULL != (p))


/**
 * Conditional global barrier.
 * Execute a \c MPI_Barrier( \c MPI_COMM_WORLD ) when flag is set.
 */
#define HU_COND_GLOBAL_BARRIER(flag)		\
    HU_COND_EXEC( flag, MPI_Barrier( MPI_COMM_WORLD ) )


#define HERC_ERROR		-100


#define XMALLOC_N(type,n)			\
    ((type*)hu_xmalloc((n * sizeof(type)), NULL))

#define XMALLOC(type)			XMALLOC(type,1)


#define XMALLOC_VAR_N(var,type,n)					\
    do {								\
	var = ((type*)hu_xmalloc((n * sizeof(type)), STRINGX(var)));	\
    } while (0)


#define XMALLOC_VAR(var,type)			\
    XMALLOC_VAR(var,type,1)


enum hu_config_option_t {
    HU_REQUIRED,
    HU_REQ_WARN,
    HU_OPTIONAL,
    HU_OPTIONAL_WARN
};

typedef enum hu_config_option_t hu_config_option_t;


/* ------------------------------------------------------------------------- *
 * Function declaration. 
 * Documentation along with function definition in .c file.
 * ------------------------------------------------------------------------- */

/** Compatibility function renaming macro */
#define solver_abort	hu_solver_abort

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Type of the parsing function for configuration variables in the
 * config file.
 */
typedef int (*hu_config_parse_fn)( const char* str, void* value );


extern void* hu_xmalloc( size_t size, const char* varname );

extern void
hu_assert_abort(
	const char*     file,
        const char*     function_name,
        int             line,
        const char*     assertion_text
        );


extern FILE*  hu_fopen( const char* filename, const char* mode );
extern int    hu_fclose( FILE* fp );
extern int    hu_fseeko( FILE* fp, off_t offset, int whence );
extern off_t  hu_ftello( FILE* fp );
extern size_t hu_fwrite( void* ptr, size_t size, size_t n, FILE* fp );
extern size_t hu_fread( void* ptr, size_t size, size_t n, FILE* fp );

extern int    hu_print_tee_va( FILE *fp, const char* format, va_list args );

extern int
hu_parsetext_v(	FILE* fp, const char* string, const char type, void* result );

extern int hu_solver_abort( const char* function_name, const char* error_msg,
			    const char* format, ... );


extern int hu_parse_double( const char* str, double* value );
extern int hu_parse_float( const char* str, float* value );
extern int hu_parse_int( const char* str, int base, int* value );

extern int
hu_config_read(
        FILE*              fp,
        const char*        key,
	hu_config_parse_fn parse_val,
	hu_config_option_t flag,
        void*		   val
        );


extern int
hu_config_get_string( FILE* fp, const char* key, hu_config_option_t flag,
		      char** val, size_t* len );

extern int
hu_config_get_string_unsafe( FILE* fp, const char* key, hu_config_option_t flag,
			     char* val );

extern int
hu_config_get_int(
	FILE*		   fp,
	const char*	   key,
	int*		   val,
	hu_config_option_t flag
	);

extern int
hu_config_get_uint(
	FILE*	           fp,
	const char*	   key,
	unsigned int*	   val,
	hu_config_option_t flag
	);


extern int
hu_config_get_double(
	FILE*	           fp,
	const char*	   key,
	double*		   val,
	hu_config_option_t flag
	);


extern int
hu_config_get_float(
	FILE*	           fp,
	const char*	   key,
	float*		   val,
	hu_config_option_t flag
	);


static inline int
hu_config_get_int_req( FILE* fp, const char* key, int* value )
{
    return hu_config_get_int( fp, key, value, HU_REQUIRED );
}


static inline int
hu_config_get_uint_req( FILE* fp, const char* key, unsigned int* value )
{
    return hu_config_get_uint( fp, key, value, HU_REQUIRED );
}


static inline int
hu_config_get_double_req( FILE* fp, const char* key, double* value )
{
    return hu_config_get_double( fp, key, value, HU_REQUIRED );
}


static inline int
hu_config_get_float_req( FILE* fp, const char* key, float* value )
{
    return hu_config_get_float( fp, key, value, HU_REQUIRED );
}


static inline int
hu_config_get_string_req( FILE* fp, const char* key, char** val, size_t* len )
{
    return hu_config_get_string( fp, key, HU_REQUIRED, val, len );
}

static inline int
hu_config_get_int_opt( FILE* fp, const char* key, int* value )
{
    return hu_config_get_int( fp, key, value, HU_OPTIONAL );
}


static inline int
hu_config_get_uint_opt( FILE* fp, const char* key, unsigned int* value )
{
    return hu_config_get_uint( fp, key, value, HU_OPTIONAL );
}


static inline int
hu_config_get_double_opt( FILE* fp, const char* key, double* value )
{
    return hu_config_get_double( fp, key, value, HU_OPTIONAL );
}


static inline int
hu_config_get_float_opt( FILE* fp, const char* key, float* value )
{
    return hu_config_get_float( fp, key, value, HU_OPTIONAL );
}


static inline int
hu_config_get_string_opt( FILE* fp, const char* key, char** val, size_t* len )
{
    return hu_config_get_string( fp, key, HU_OPTIONAL, val, len );
}


static inline int
read_config_string2( FILE* fp, const char* key, char* val, size_t len )
{
    return hu_config_get_string_req( fp, key, &val, &len );
}

/**
 * Free the memory pointed by *ptr and then set *ptr to NULL.
 */
static inline void
xfree( void** ptr )
{
    if (ptr != NULL) {
	if (NULL != *ptr) {
	    free( *ptr );
	}
	*ptr = NULL;
    }
}


static inline void
xfree_char( char** ptr )
{
    xfree( (void**)ptr );
}


#ifdef __cplusplus
} /* extern "C" { */
#endif

#endif /* HERC_UTIL_H */
