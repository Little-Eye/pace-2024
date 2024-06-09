#ifndef _SuiteOPT_H_
#define _SuiteOPT_H_
#include <math.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#ifdef MATLAB_MEX_FILE
#include "matrix.h"
#include "mex.h"
#endif

#define SOPT_OUT_OF_MEMORY                                      (901)
#define SOPT_ERROR_IN_INPUT_MATRIX                              (902)
#define SOPT_MATRIX_ELEMENT_WAS_ZERO                            (903)
#define SOPT_START_MESSAGES                                     (900)
#define SOPT_END_MESSAGES                                       (999)

/* some constants used by SOPT */

#define SuiteOPTfalse 0
#define SuiteOPTtrue 1

#ifndef NULL
#define NULL 0
#endif

#define SOPTZERO ((SOPTFLOAT) 0)

/* by default, sopt_timer (clock_gettime) is enabled */
//#define SUITEOPT_TIMER

/* define the long version of integers */
#define LONG long

/* define the integer precision for the BLAS */
#define BLAS_INT long

/* SOPTFLOAT is the default precision of floating point variables */
#define SOPTFLOAT double

/* DLONG is a compiler flag that is defined when compiling with MATLAB.
   MATLAB use doubles and long ints. */
#ifdef DLONG
#define SOPTFLOAT double
#define SOPTINT LONG
#define SuiteOPTinfint LONG_MAX
#define CHOLMODlong SuiteOPTtrue

#else
/* Otherwise, select precision by commenting out one pair of definitions and
   keeping the other.  SOPTINT is default precision of integers; CHOLMODlong
   is true if using long integers in SuitOPT, otherwise it is false */

/* standard ints */
#define SOPTINT int
#define CHOLMODlong SuiteOPTfalse
#define SuiteOPTinfint INT_MAX

/* long ints */
/*
#define SOPTINT long
#define CHOLMODlong SuiteOPTtrue
#define SuiteOPTinfint LONG_MAX
*/

#endif

/* ANSI C99 has a clean definition of IEEE infinity, defined in <math.h>.
   MATLAB has this as well.  With the ANSI (C90) version of C, there is no
   well-defined infinity, so DBL_MAX is used instead.

   You can override these defaults and define your own version of infinity
   with (for example):

   cc -ansi -DSuiteOPTinf=1e200 pproj.c ...
*/

/* infinite float */
#ifndef SuiteOPTinf
#ifdef INFINITY
/* ANSI C99 (gcc -std=c99) */
#define SuiteOPTinf INFINITY
#else
/* ANSI C90 (gcc -ansi) */
#define SuiteOPTinf DBL_MAX
#endif
#endif

#define SOPTMAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define SOPTMIN(a,b) ( ((a) < (b)) ? (a) : (b) )

#define EMPTY (SOPTINT) -1

/* debugging options */
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE
#define ASSERT(expression) (mxAssert ((expression), ""))
#else
#define ASSERT(expression) (assert (expression))
#endif
#else
#define ASSERT(expression)
#endif

/* When using long its, need to include "_l" when calling CHOLMOD routines */
#if CHOLMODlong
#define CHOLMOD(name) cholmod_l_ ## name
#else
#define CHOLMOD(name) cholmod_ ## name
#endif

/* By default, the BLAS are not used in SuiteOPT (they are used in
   SuiteSparse). To use the BLAS in SuiteOPT, comment out the next
   line and then below, specify whether your BLAS routines include an
   underscore, and set the threshhold dimension for using the BLAS.  */
#define NOBLAS

#ifndef NOBLAS

/* when BLAS are used, comment out the next statement if there is no
   underscore in the subroutine names */
#define BLAS_UNDERSCORE

#ifdef BLAS_UNDERSCORE

#define SOPT_DGEMV dgemv_
#define SOPT_DAXPY daxpy_
#define SOPT_DDOT ddot_
#define SOPT_DSCAL dscal_
#define SOPT_DCOPY dcopy_
#define SOPT_IDAMAX idamax_

#else

#define SOPT_DGEMV dgemv
#define SOPT_DAXPY daxpy
#define SOPT_DDOT ddot
#define SOPT_DSCAL dscal
#define SOPT_DCOPY dcopy
#define SOPT_IDAMAX idamax

#endif

/* define the starting size of vectors when BLAS are used */
#define MATVEC_START 1000000
#define DAXPY_START  1000000
#define DDOT_START   1000000
#define DSCAL_START  1000000
#define DCOPY_START  1000000
#define IDAMAX_START 1000000

/* protypes for the BLAS routines that could be used */
void SOPT_DGEMV (char *trans, BLAS_INT *m, BLAS_INT *n, SOPTFLOAT *alpha,
        SOPTFLOAT *A, BLAS_INT *lda, SOPTFLOAT *X, BLAS_INT *incx,
        SOPTFLOAT *beta, SOPTFLOAT *Y, BLAS_INT *incy) ;

void SOPT_DAXPY (BLAS_INT *n, SOPTFLOAT *DA, SOPTFLOAT *DX, BLAS_INT *incx,
        SOPTFLOAT *DY, BLAS_INT *incy) ;

SOPTFLOAT SOPT_DDOT (BLAS_INT *n, SOPTFLOAT *DX, BLAS_INT *incx, SOPTFLOAT *DY,
        BLAS_INT *incy) ;

void SOPT_DSCAL (BLAS_INT *n, SOPTFLOAT *DA, SOPTFLOAT *DX, BLAS_INT *incx) ;

void SOPT_DCOPY (BLAS_INT *n, SOPTFLOAT *DX, BLAS_INT *incx, SOPTFLOAT *DY,
        BLAS_INT *incy) ;

BLAS_INT SOPT_IDAMAX (BLAS_INT *n, SOPTFLOAT *DX, BLAS_INT *incx) ;

#endif

/* prototypes for codes in SIOPT.c */
void sopt_free
(
    void * p
) ;

void * sopt_malloc
(
    int *status,
    SOPTINT   n,
    int    size
) ;

void sopt_error
(
    int status,
    const char *file,
    int line,
    const char *message
) ;

double sopt_timer ( void ) ;

void sopt_print_TF
(
    int TF /* TRUE or FALSE */
) ;

void sopt_printA
(
    SOPTINT  ncol, /* number of cols in A */
    SOPTINT   *Ap, /* size ncol+1, column pointers */
    SOPTINT  *Anz, /* if NULL, A is packed; otherwise gives # nonzeros in cols*/
    SOPTINT   *Ai, /* size Ap [ncol], row indices for A */
    SOPTFLOAT *Ax, /* size Ap [ncol], numerical entries of A */
    char    *what  /* name of the matrix */
) ;

void sopt_printAMATLAB
(
    SOPTINT   const ncol, /* number of cols in A */
    SOPTINT   const  *Ap, /* size ncol+1, column pointers */
    SOPTINT   const  *Ai, /* size Ap [ncol], row indices for A */
    SOPTINT   const *Anz, /* if NULL, A packed; otherwise # nonzeros in cols */
    SOPTFLOAT const  *Ax, /* size Ap [ncol], numerical entries of A */
    char           *what  /* name of the matrix */
) ;

void sopt_printx
(
    SOPTFLOAT const *x, /* numerical entries in the vector */
    SOPTINT   const  n, /* dimension of the vector */
    char          *what  /* name of the vector */
) ;

void sopt_printxMATLAB
(
    SOPTFLOAT *x, /* numerical entries in the vector */
    SOPTINT    n, /* dimension of the vector */
    char   *what  /* name of the vector */
) ;

void sopt_printi
(
    SOPTINT *i, /* array of SOPTINT's */
    SOPTINT  n, /* dimension of i */
    char *what  /* name of array */
) ;

void sopt_printiMATLAB
(
    SOPTINT *i, /* array of SOPTINT's */
    SOPTINT  n, /* dimension of i */
    char *what  /* name of the array */
) ;

void sopt_print_int
(
    int     *i, /* array of int's */
    SOPTINT  n, /* dimension of i */
    char *what  /* name of array */
) ;

void sopt_print_intMATLAB
(
    int     *i, /* array of int's */
    SOPTINT  n, /* dimension of i */
    char *what  /* name of the array */
) ;

void sopt_add
(
    SOPTFLOAT       *x,  /* array to which s is added */
    SOPTFLOAT const  s,  /* scalar */
    SOPTINT   const  n   /* length of x */
) ;

void sopt_addi
(
    SOPTINT       *x,  /* array to which s is added */
    SOPTINT const  s,  /* scalar */
    SOPTINT const  n   /* length of x */
) ;

void sopt_scale
(
    SOPTFLOAT       *x,  /* array to be scaled */
    SOPTFLOAT const *y,  /* array used for the scaling */
    SOPTFLOAT const  s,  /* scale */
    SOPTINT   const  n   /* length of x */
) ;

SOPTFLOAT sopt_scale_max
(
    SOPTFLOAT       *x,  /* scaled array */
    SOPTFLOAT const *y,  /* array used for the scaling */
    SOPTFLOAT const  s,  /* scale */
    SOPTINT   const  n   /* length of x */
) ;

void sopt_step
(
    SOPTFLOAT       *xnew, /* updated x vector */
    SOPTFLOAT const    *x, /* current x */
    SOPTFLOAT const    *d, /* search direction */
    SOPTFLOAT const alpha, /* stepsize */
    SOPTINT   const     n  /* dimension */
) ;

SOPTFLOAT sopt_step_max
(
    SOPTFLOAT       *xnew, /* updated x vector */
    SOPTFLOAT const    *x, /* current x */
    SOPTFLOAT const    *d, /* search direction */
    SOPTFLOAT const alpha, /* stepsize */
    SOPTINT   const     n  /* dimension */
) ;

void sopt_daxpy
(
    SOPTFLOAT       *x, /* input and output vector */
    SOPTFLOAT const *d, /* direction vector */
    SOPTFLOAT const  s, /* stepsize */
    SOPTINT   const  n  /* length of the vectors */
) ;

void sopt_copyx
(
    SOPTFLOAT       *x, /* output of copy */
    SOPTFLOAT const *y, /* input of copy */
    SOPTINT   const  n  /* length of vectors */
) ;

void sopt_copyx_noblas
(
    SOPTFLOAT       *x, /* output of copy */
    SOPTFLOAT const *y, /* input of copy */
    SOPTINT   const  n  /* length of vectors */
) ;

void sopt_copyi
(
    SOPTINT       *x, /* output of copy */
    SOPTINT const *y, /* input of copy */
    SOPTINT const  n  /* length of vectors */
) ;

void sopt_copy_int
(
    int           *x, /* output of copy */
    int     const *y, /* input of copy */
    SOPTINT const  n  /* length of vectors */
) ;

SOPTFLOAT sopt_dot
(
    SOPTFLOAT const *x, /* first vector */
    SOPTFLOAT const *y, /* second vector */
    SOPTINT   const  n  /* length of vectors */
) ;

void sopt_initx
(
    SOPTFLOAT      *x,  /* array to be initialized */
    SOPTFLOAT const s,  /* scalar */
    SOPTINT   const n   /* length of x */
) ;

void sopt_initi
(
    SOPTINT      *x,  /* array to be initialized */
    SOPTINT const s,  /* scalar */
    SOPTINT const n   /* length of x */
) ;

void sopt_init_int
(
    int          *x,  /* array to be initialized */
    int     const s,  /* scalar */
    SOPTINT const n   /* length of x */
) ;

SOPTFLOAT sopt_sup_normx
(
    SOPTFLOAT const *x, /* vector */
    SOPTINT   const  n  /* length of vector */
) ;

SOPTINT sopt_supi
(
    SOPTINT const *x, /* vector */
    SOPTINT const  n  /* length of vector */
) ;

void sopt_transpose
(
    SOPTINT          *Bp, /* size nrow+1, column pointers (output) */
    SOPTINT          *Bi, /* size Ap [ncol], row indices of B (output) */
    SOPTFLOAT        *Bx, /* size Ap [ncol], numerical entries of B (output) */
    SOPTINT   const  *Ap, /* size ncol+1, column pointers */
    SOPTINT   const  *Ai, /* size Ap [ncol], row indices for A */
    SOPTFLOAT const  *Ax, /* size Ap [ncol], numerical entries of A */
    SOPTINT   const nrow, /* number of rows in A */
    SOPTINT   const ncol, /* number of cols in A */
    SOPTINT           *W  /* work array of size nrow */
) ;

int sopt_sort_cols
(
    SOPTINT    *Ap, /* column pointers */
    SOPTINT    *Ai, /* row indices */
    SOPTFLOAT  *Ax, /* numerical values */
    SOPTINT   *atp, /* row pointers for transpose */
    SOPTINT   *ati, /* column indices for transpose */
    SOPTFLOAT *atx, /* numerical values for transpose */
    SOPTINT   nrow, /* number of rows */
    SOPTINT   ncol  /* number of cols */
) ;

int sopt_convert_dense_to_sparse
(
    SOPTINT            **Ap, /* return column pointers */
    SOPTINT            **Ai, /* return row indices */
    SOPTFLOAT          **Ax, /* return numerical entries */
    SOPTFLOAT const      *A, /* matrix entries in dense format */
    SOPTINT   const    nrow, /* number of rows in A that are used */
    SOPTINT   const    ncol, /* number of columns in A that are used */
    int       const by_rows  /* T => matrix stored by rows, F => by columns */
) ;

int sopt_convert_triple_to_sparse /* returned integer:
                                   0 conversion successful
                                   1 out-of-memory
                                   2 error in input matrix */
(
    SOPTINT         **Ap, /* column pointers */
    SOPTINT         **Ai, /* row indices (increasing order in each column) */
    SOPTFLOAT       **Ax, /* numerical entries */
    SOPTINT        *nrow, /* return 1 + largest element of Ti */
    SOPTINT        *ncol, /* return 1 + largest element of Tj */
    SOPTINT   const  *Ti, /* row    indices   of nonzero entries */
    SOPTINT   const  *Tj, /* column indices   of nonzero entries */
    SOPTFLOAT const  *Tx, /* nonzero numerical values */
    SOPTINT   const  nnz, /* number of nonzeros in the matrix */
    int const order_cols, /* TRUE if column in increasing order */
    int const        sym  /* TRUE if matrix symmetric and only elements on the
                             diagonal and on one side of the diagonal given */
) ;

int sopt_check_matrix /* return 1 if an error was detected, otherwise return 0*/
(
    SOPTINT   const  *Ap, /* column pointers */
    SOPTINT   const  *Ai, /* row indices */
    SOPTFLOAT const  *Ax, /* numerical entries */
    SOPTINT   const ncol  /* number of columns in matrix */
) ;

void sopt_minsortx
(
    SOPTINT         *y, /* n-by-1 (output) */
    SOPTFLOAT const *x, /* n-by-1 (input not modified) */
    SOPTINT         *w, /* n-by-1, (input, working array) */
    SOPTINT          n  /* number of elements to sort */
) ;

void sopt_minsorti
(
    SOPTINT        *y, /* n-by-1 (output) */
    SOPTINT  const *x, /* n-by-1 (input not modified) */
    SOPTINT        *w, /* n-by-1, (input, working array) */
    SOPTINT         n  /* number of elements to sort */
) ;

#endif
