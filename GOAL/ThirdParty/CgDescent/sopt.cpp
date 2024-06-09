#include "sopt.h"
#define SOPTZERO ((SOPTFLOAT) 0)
#define SOPTONE  ((SOPTFLOAT) 1)

/* =========================================================================
   ============================== sopt_free ================================
   ========================================================================= */
void sopt_free
(
    void * p
)
{

#ifdef MATLAB_MEX_FILE
    mxFree (p) ;
#else
    free (p) ;
#endif

}

/* =========================================================================
   ============================== sopt_malloc ==============================
   ========================================================================= */
void * sopt_malloc
(
    int *status,
    SOPTINT   n,
    int    size
)
{
    void *p ;
    if ( n > 0 )
    {
#ifdef MATLAB_MEX_FILE
        p = mxMalloc ((LONG) n * size) ;
#else
        p = malloc ((LONG) n * size) ;
#endif
    }
    else
    {
        return (NULL) ;
    }
    if      ( p == NULL )                     *status = SOPT_OUT_OF_MEMORY ;
    else if ( *status != SOPT_OUT_OF_MEMORY ) *status = 0 ;
    return (p) ;
}

/* ========================================================================== */
/* ====== sopt_error ======================================================== */
/* ========================================================================== */
/* when -g compiler option is used, prints line number of error
   ========================================================================== */
void sopt_error
(
    int status,
    const char *file,
    int line,
    const char *message
)
{
    if (status < 0)
    {
        printf ("file: %s line: %d status: %d %s\n",
                 file, line, status, message) ;
        fflush (stdout) ;
#ifdef MATLAB_MEX_FILE
        mexErrMsgTxt (message) ;
#else
        ASSERT (0) ;
        abort ( ) ;
#endif
    }
}

/* ========================================================================== */
/* ====== sopt_timer ======================================================== */
/* ========================================================================== */
/* Returns current walltime in seconds.  */
/* ========================================================================== */

#ifdef SUITEOPT_TIMER
#include <time.h>

double sopt_timer ( void )
{
    struct timespec t ;
    clock_gettime (CLOCK_MONOTONIC, &t) ;
    return ( (double) t.tv_sec + 1.e-9*((double) t.tv_nsec) ) ;
}

#else

double sopt_timer ( void )
{
    return ( (double) 0 ) ;
}
#endif

/* ==========================================================================
   === sopt_print_TF ========================================================
   ==========================================================================
    Print TRUE if TF is TRUE, otherwise print FALSE
   ========================================================================== */
void sopt_print_TF 
(
    int TF /* TRUE or FALSE */
)
{
    if ( TF == SuiteOPTtrue )
    {
        printf ("TRUE\n") ;
    }
    else if ( TF == SuiteOPTfalse )
    {
        printf ("FALSE\n") ;
    }
    else
    {
        printf ("EMPTY\n") ;
    }
}

/* ==========================================================================
   === sopt_printA ==========================================================
   ==========================================================================
    Print a sparse matrix
   ========================================================================== */
void sopt_printA
(
    SOPTINT  ncol, /* number of cols in A */
    SOPTINT   *Ap, /* size ncol+1, column pointers */
    SOPTINT  *Anz, /* if NULL, A is packed; otherwise gives # nonzeros in cols*/
    SOPTINT   *Ai, /* size Ap [ncol], row indices for A */
    SOPTFLOAT *Ax, /* size Ap [ncol], numerical entries of A */
    char    *what  /* name of the matrix */
)
{
    SOPTINT j, p, q ;
    printf ("%s =\n", what) ;
    if ( Anz == NULL ) /* matrix is packed */
    {
        p = 0 ;
        for (j = 0; j < ncol; j++)
        {
            q = Ap [j+1] ;
            for (; p < q; p++)
            {
                printf ("%ld %ld %e\n", (LONG) Ai [p], (LONG) j, Ax [p]) ;
            }
        } 
    }
    else /* Anz gives the number of nonzeros in each column of A */
    {
        for (j = 0; j < ncol; j++)
        {
            p = Ap [j] ;
            q = Ap [j] + Anz [j] ;
            for (; p < q; p++)
            {
                printf ("%ld %ld %e\n", (LONG) Ai [p], (LONG) j, Ax [p]) ;
            }
        }
    }
}

/* ==========================================================================
   === sopt_printAMATLAB ====================================================
   ==========================================================================
    Print text that can be fed to MATLAB to construct a sparse matrix.
    The sparse matrix is constructed from triples and rows/columns start at 1.
   ========================================================================== */
void sopt_printAMATLAB
(
    SOPTINT   const ncol, /* number of cols in A */
    SOPTINT   const  *Ap, /* size ncol+1, column pointers */
    SOPTINT   const  *Ai, /* size Ap [ncol], row indices for A */
    SOPTINT   const *Anz, /* if NULL, A packed; otherwise # nonzeros in cols */
    SOPTFLOAT const  *Ax, /* size Ap [ncol], numerical entries of A */
    char           *what  /* name of the matrix */
)
{
    SOPTINT j, p, q ;
    printf ("Temp = [\n") ;
    if ( Anz != NULL )
    {
        for (j = 0; j < ncol; j++)
        {
            p = Ap [j] ;
            q = Ap [j] + Anz [j] ;
            for (; p < q; p++)
            {
                printf ("%ld %ld %25.15e\n",
                        (LONG) Ai [p]+1, (LONG) j+1, Ax [p]) ;
            }
        }
    }
    else
    {
        for (j = 0; j < ncol; j++)
        {
            p = Ap [j] ;
            q = Ap [j+1] ;
            for (; p < q; p++)
            {
                printf ("%ld %ld %25.15e\n",
                        (LONG) Ai [p]+1, (LONG) (j+1), Ax [p]) ;
            }
        }
    }
    printf ("] ;\n") ;
    printf ("%s = sparse (Temp(:, 1), Temp(:, 2), Temp(:, 3)) ;\n", what) ;
}

/* ==========================================================================
   === sopt_printx ==========================================================
   ==========================================================================
    Print a float vector.
   ========================================================================== */
void sopt_printx
(
    SOPTFLOAT const *x, /* numerical entries in the vector */
    SOPTINT   const  n, /* dimension of the vector */
    char          *what  /* name of the vector */
)
{
    SOPTINT i ;
    printf ("%s =\n", what) ;
    for (i = 0; i < n; i++)
    {
        printf ("%ld %25.15e\n", (LONG) i, x [i]) ;
    }
}

/* ==========================================================================
   === sopt_printxMATLAB ====================================================
   ==========================================================================
    Print text that can be fed to MATLAB to construct a dense vector.
   ========================================================================== */
void sopt_printxMATLAB
(
    SOPTFLOAT *x, /* numerical entries in the vector */
    SOPTINT    n, /* dimension of the vector */
    char   *what  /* name of the vector */
)
{
    SOPTINT j ;
    printf ("%s = [\n", what) ;
    for (j = 0; j < n; j++)
    {
        printf ("%26.16e\n", x [j]) ;
    }
    printf ("] ;\n") ;
}

/* ==========================================================================
   === sopt_printi ==========================================================
   ==========================================================================
    Print a SOPTINT vector.
   ========================================================================== */
void sopt_printi
(
    SOPTINT *i, /* array of SOPTINT's */
    SOPTINT  n, /* dimension of i */
    char *what  /* name of array */
)
{
    SOPTINT j ;
    printf ("%s =\n", what) ;
    for (j = 0; j < n; j++)
    {
        printf ("%ld %ld\n", (LONG) j, (LONG) i [j]) ;
    }
}

/* ==========================================================================
   === sopt_print_int =======================================================
   ==========================================================================
    Print an int vector.
   ========================================================================== */
void sopt_print_int 
(
    int     *i, /* array of int's */
    SOPTINT  n, /* dimension of i */
    char *what  /* name of array */
)
{
    SOPTINT j ;
    printf ("%s =\n", what) ;
    for (j = 0; j < n; j++)
    {
        printf ("%ld %i\n", (LONG) j, i [j]) ;
    }
}

/* ==========================================================================
   === sopt_print_iMATLAB ====================================================
   ==========================================================================
    Print text that can be fed to MATLAB to construct a dense int vector.
   ========================================================================== */
void sopt_printiMATLAB
(
    SOPTINT *i, /* array of SOPTINT's */
    SOPTINT  n, /* dimension of i */
    char *what  /* name of the array */
)
{
    SOPTINT j ;
    printf ("%s = [\n", what) ;
    for (j = 0; j < n; j++)
    {
        printf ("%ld\n", (LONG) i [j]) ;
    }
    printf ("] ;\n") ;
}

/* =========================================================================
   ===================== sopt_add ==========================================
   =========================================================================
   add scalar s to each component of array x
   ========================================================================= */

void sopt_add
(
    SOPTFLOAT       *x,  /* array to which s is added */
    SOPTFLOAT const  s,  /* scalar */
    SOPTINT   const  n   /* length of x */
)
{
    SOPTINT j, n5 ;
    SOPTFLOAT *xj ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] += s ;
    xj = x+j ;
    for (; j < n; j += 5)
    {
        *(xj++) += s ;
        *(xj++) += s ;
        *(xj++) += s ;
        *(xj++) += s ;
        *(xj++) += s ;
    }
}

/* =========================================================================
   ===================== sopt_addi =========================================
   =========================================================================
   add scalar s to each component of array x
   ========================================================================= */

void sopt_addi
(
    SOPTINT       *x,  /* array to which s is added */
    SOPTINT const  s,  /* scalar */
    SOPTINT const  n   /* length of x */
)
{
    SOPTINT j, n5 ;
    SOPTINT *xj ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] += s ;
    xj = x+j ;
    for (; j < n; j += 5)
    {
        *(xj++) += s ;
        *(xj++) += s ;
        *(xj++) += s ;
        *(xj++) += s ;
        *(xj++) += s ;
    }
}

/* =========================================================================
   ===================== sopt_scale ========================================
   =========================================================================
   Scale a SOPTFLOAT array x = s*y
   ========================================================================= */

void sopt_scale
(
    SOPTFLOAT       *x,  /* scaled array */
    SOPTFLOAT const *y,  /* array used for the scaling */
    SOPTFLOAT const  s,  /* scale */
    SOPTINT   const  n   /* length of x */
)
{
    if ( x == y )
    {
        SOPTFLOAT *xj ;
        if ( s == SOPTONE ) return ; /* nothing to do */
#ifndef NOBLAS
        if ( n >= DSCAL_START )
        {
            SOPTFLOAT S = s ;
            BLAS_INT int_one = 1 ;
            BLAS_INT N = (BLAS_INT) n ;
            SOPT_DSCAL (&N, &S, x, &int_one) ;
        }
        else
#endif
        {
            SOPTINT j, n5 ;
            n5 = n % 5 ;
            for (j = 0; j < n5; j++) x [j] *= s ;
            xj = x+j ;
            for (; j < n; j += 5)
            {
                *(xj++) *= s ;
                *(xj++) *= s ;
                *(xj++) *= s ;
                *(xj++) *= s ;
                *(xj++) *= s ;
            }
        }
    }
    else /* x != y */
    {
        if ( s == SOPTONE ) /* this is a copy */
        {
            sopt_copyx (x, y, n) ;
        }
        else if ( s == -SOPTONE ) /* s = -1 */
        {
            SOPTINT j, n5 ;
            n5 = n % 5 ;
            for (j = 0; j < n5; j++) x [j] = -y [j] ;
            for (; j < n; j += 5)
            {
                x [j]   = -y [j] ;
                x [j+1] = -y [j+1] ;
                x [j+2] = -y [j+2] ;
                x [j+3] = -y [j+3] ;
                x [j+4] = -y [j+4] ;
            }
        }
        else
        {
            SOPTINT j, n5 ;
            n5 = n % 5 ;
            for (j = 0; j < n5; j++) x [j] = s*y [j] ;
            for (; j < n; j += 5)
            {
                x [j]   = s*y [j] ;
                x [j+1] = s*y [j+1] ;
                x [j+2] = s*y [j+2] ;
                x [j+3] = s*y [j+3] ;
                x [j+4] = s*y [j+4] ;
            }
        }
    }
}

/* =========================================================================
   ===================== sopt_scale_max ====================================
   =========================================================================
   Scale a SOPTFLOAT array x = s*y and also evaluate ||y||_sup
   ========================================================================= */

SOPTFLOAT sopt_scale_max
(
    SOPTFLOAT       *x,  /* scaled array */
    SOPTFLOAT const *y,  /* array used for the scaling */
    SOPTFLOAT const  s,  /* scale */
    SOPTINT   const  n   /* length of x */
)
{
    SOPTFLOAT ymax, *xj ;
    ymax = SOPTZERO ;
    if ( x == y )
    {
        if ( s == SOPTONE ) /* nothing to do except compute the max */
        {
            ymax = sopt_sup_normx (y, n) ;
        }
#ifndef NOBLAS
        else if ( n >= DSCAL_START )
        {
            ymax = sopt_sup_normx (y, n) ;
            SOPTFLOAT S = s ;
            BLAS_INT int_one = 1 ;
            BLAS_INT N = (BLAS_INT) n ;
            SOPT_DSCAL (&N, &S, x, &int_one) ;
        }
#endif
        else
        {
            SOPTINT j, n5 ;
            n5 = n % 5 ;
            for (j = 0; j < n5; j++)
            {
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] *= s ;
            }
            xj = x+j ;
            for (; j < n; j += 5)
            {
                if ( ymax < fabs(*xj) ) ymax = fabs(*xj) ;
                *(xj++) *= s ;
                if ( ymax < fabs(*xj) ) ymax = fabs(*xj) ;
                *(xj++) *= s ;
                if ( ymax < fabs(*xj) ) ymax = fabs(*xj) ;
                *(xj++) *= s ;
                if ( ymax < fabs(*xj) ) ymax = fabs(*xj) ;
                *(xj++) *= s ;
                if ( ymax < fabs(*xj) ) ymax = fabs(*xj) ;
                *(xj++) *= s ;
            }
        }
    }
    else /* x != y */
    {
        if ( s == SOPTONE ) /* this is a copy */
        {
            ymax = sopt_sup_normx (y, n) ;
            sopt_copyx (x, y, n) ;
        }
        else if ( s == -SOPTONE ) /* s = -1 */
        {
            SOPTINT j, n5 ;
            n5 = n % 5 ;
            for (j = 0; j < n5; j++)
            {
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = -y [j] ;
            }
            for (; j < n; )
            {
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = -y [j] ; j++ ;
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = -y [j] ; j++ ;
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = -y [j] ; j++ ;
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = -y [j] ; j++ ;
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = -y [j] ; j++ ;
            }
        }
        else
        {
            SOPTINT j, n5 ;
            n5 = n % 5 ;
            for (j = 0; j < n5; j++)
            {
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = s*y [j] ;
            }
            for (; j < n; )
            {
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = s*y [j] ; j++ ;
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = s*y [j] ; j++ ;
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = s*y [j] ; j++ ;
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = s*y [j] ; j++ ;
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = s*y [j] ; j++ ;
            }
        }
    }
    return (ymax) ;
}

/* =========================================================================
   === sopt_step ============================================================
   =========================================================================
    Set xnew = x + alpha*d
   ========================================================================= */
void sopt_step
(
    SOPTFLOAT       *xnew, /* updated x vector */
    SOPTFLOAT const    *x, /* current x */
    SOPTFLOAT const    *d, /* search direction */
    SOPTFLOAT const alpha, /* stepsize */
    SOPTINT   const     n  /* dimension */
)
{
    SOPTINT j, n5 ;

    if ( x == NULL ) return ;
    n5 = n % 5 ;     /* n5 = n mod 5 */
    /* check if step size equals 1 */
    if ( alpha == SOPTONE )
    {
        for (j = 0; j < n5; j++)
        {
            xnew [j] = x [j] + d [j] ;
        }
        if ( xnew == x )
        {
            for (; j < n; )
            {
                 xnew [j] += d [j] ; j++ ;
                 xnew [j] += d [j] ; j++ ;
                 xnew [j] += d [j] ; j++ ;
                 xnew [j] += d [j] ; j++ ;
                 xnew [j] += d [j] ; j++ ;
            }
        }
        else
        {
            for (; j < n; )
            {
                 xnew [j] = x [j] + d [j] ; j++ ;
                 xnew [j] = x [j] + d [j] ; j++ ;
                 xnew [j] = x [j] + d [j] ; j++ ;
                 xnew [j] = x [j] + d [j] ; j++ ;
                 xnew [j] = x [j] + d [j] ; j++ ;
            }
        }
    }
    else if ( alpha == -SOPTONE )
    {
        for (j = 0; j < n5; j++)
        {
            xnew [j] = x [j] - d [j] ;
        }
        for (; j < n; )
        {
            xnew [j] = x [j] - d [j] ; j++ ;
            xnew [j] = x [j] - d [j] ; j++ ;
            xnew [j] = x [j] - d [j] ; j++ ;
            xnew [j] = x [j] - d [j] ; j++ ;
            xnew [j] = x [j] - d [j] ; j++ ;
        }
    }
    /* else step size is not 1 */
    else
    {
        for (j = 0; j < n5; j++)
        {
            xnew [j] = x [j] + alpha*d [j] ;
        }
        for (; j < n; )
        {
            xnew [j] = x [j] + alpha*d [j] ; j++ ;
            xnew [j] = x [j] + alpha*d [j] ; j++ ;
            xnew [j] = x [j] + alpha*d [j] ; j++ ;
            xnew [j] = x [j] + alpha*d [j] ; j++ ;
            xnew [j] = x [j] + alpha*d [j] ; j++ ;
        }
    }
}

/* =========================================================================
   === sopt_step_max =======================================================
   =========================================================================
    Set xnew = x + alpha*d, return dmax = ||d||_sup
   ========================================================================= */
SOPTFLOAT sopt_step_max
(
    SOPTFLOAT       *xnew, /* updated x vector */
    SOPTFLOAT const    *x, /* current x */
    SOPTFLOAT const    *d, /* search direction */
    SOPTFLOAT const alpha, /* stepsize */
    SOPTINT   const     n  /* dimension */
)
{
    SOPTFLOAT dmax ;
    SOPTINT j, n5 ;

    dmax = SOPTZERO ;
    n5 = n % 5 ;     /* n5 = n mod 5 */
    /* check if step size equals 1 */
    if ( alpha == SOPTONE )
    {
        for (j = 0; j < n5; j++)
        {
            xnew [j] = x [j] + d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
        }
        if ( xnew == x )
        {
            for (; j < n; )
            {
                xnew [j] += d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
                xnew [j] += d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
                xnew [j] += d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
                xnew [j] += d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
                xnew [j] += d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
            }
        }
        else
        {
            for (; j < n; )
            {
                xnew [j] = x [j] + d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
                xnew [j] = x [j] + d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
                xnew [j] = x [j] + d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
                xnew [j] = x [j] + d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
                xnew [j] = x [j] + d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
            }
        }
    }
    else if ( alpha == -SOPTONE )
    {
        for (j = 0; j < n5; j++)
        {
            xnew [j] = x [j] - d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
        }
        for (; j < n; )
        {
            xnew [j] = x [j] - d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
            xnew [j] = x [j] - d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
            xnew [j] = x [j] - d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
            xnew [j] = x [j] - d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
            xnew [j] = x [j] - d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
        }
    }
    /* else step size is not 1 */
    else
    {
        for (j = 0; j < n5; j++)
        {
            xnew [j] = x [j] + alpha*d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
        }
        for (; j < n; )
        {
            xnew [j] = x [j] + alpha*d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
            xnew [j] = x [j] + alpha*d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
            xnew [j] = x [j] + alpha*d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
            xnew [j] = x [j] + alpha*d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
            xnew [j] = x [j] + alpha*d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
        }
    }
    return (dmax) ;
}

/* =========================================================================
   ==== sopt_daxpy =========================================================
   =========================================================================
   Compute x = x + s * d
   ========================================================================= */
void sopt_daxpy
(
    SOPTFLOAT       *x, /* input and output vector */
    SOPTFLOAT const *d, /* direction vector */
    SOPTFLOAT const  s, /* stepsize */
    SOPTINT   const  n  /* length of the vectors */
)
{
#ifndef NOBLAS
    if ( n >= DAXPY_START )
    {
        SOPTFLOAT S = s ;
        BLAS_INT int_one = 1 ;
        BLAS_INT N = (BLAS_INT) n ;
        SOPT_DAXPY (&N, &S, (SOPTFLOAT *) d, &int_one, x, &int_one) ;
    }
    else
#endif
    {
        SOPTINT i, n5 ;
        n5 = n % 5 ;
        if ( s == -SOPTONE)
        {
            for (i = 0; i < n5; i++) x [i] -= d [i] ;
            for (; i < n; i += 5)
            {
                x [i]   -= d [i] ;
                x [i+1] -= d [i+1] ;
                x [i+2] -= d [i+2] ;
                x [i+3] -= d [i+3] ;
                x [i+4] -= d [i+4] ;
            }
        }
        else if ( s == SOPTONE )
        {
            for (i = 0; i < n5; i++) x [i] += d [i] ;
            for (; i < n; i += 5)
            {
                x [i]   += d [i] ;
                x [i+1] += d [i+1] ;
                x [i+2] += d [i+2] ;
                x [i+3] += d [i+3] ;
                x [i+4] += d [i+4] ;
            }
        }
        else
        {
            for (i = 0; i < n5; i++) x [i] += s*d [i] ;
            for (; i < n; i += 5)
            {
                x [i]   += s*d [i] ;
                x [i+1] += s*d [i+1] ;
                x [i+2] += s*d [i+2] ;
                x [i+3] += s*d [i+3] ;
                x [i+4] += s*d [i+4] ;
            }
        }
    }
}

/* =========================================================================
   === sopt_copyx ==========================================================
   =========================================================================
   Copy SOPTFLOAT vector y into vector x
   ========================================================================= */
void sopt_copyx
(
    SOPTFLOAT       *x, /* output of copy */
    SOPTFLOAT const *y, /* input of copy */
    SOPTINT   const  n  /* length of vectors */
)
{
    if ( (y == x) || (y == NULL) ) return ;
#ifndef NOBLAS
    if ( n >= DCOPY_START )
    {
        BLAS_INT int_one  = 1 ;
        BLAS_INT N = (BLAS_INT) n ;
        SOPT_DCOPY (&N, (SOPTFLOAT *) y, &int_one, x, &int_one) ;
    }
    else
#endif
    {
        SOPTINT i, n5 ;
        n5 = n % 5 ;
        for (i = 0; i < n5; i++) x [i] = y [i] ;
        for (; i < n; i += 5)
        {
            x [i]   = y [i] ;
            x [i+1] = y [i+1] ;
            x [i+2] = y [i+2] ;
            x [i+3] = y [i+3] ;
            x [i+4] = y [i+4] ;
        }
    }
}

/* =========================================================================
   === sopt_copyx_noblas ===================================================
   =========================================================================
   Copy SOPTFLOAT vector y into vector x
   NOTE: when y is a part of x, cannot use BLAS
   ========================================================================= */
void sopt_copyx_noblas
(
    SOPTFLOAT       *x, /* output of copy */
    SOPTFLOAT const *y, /* input of copy */
    SOPTINT   const  n  /* length of vectors */
)
{
    if ( (y == x) || (y == NULL) ) return ;
    SOPTINT i, n5 ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++) x [i] = y [i] ;
    for (; i < n; i += 5)
    {
        x [i]   = y [i] ;
        x [i+1] = y [i+1] ;
        x [i+2] = y [i+2] ;
        x [i+3] = y [i+3] ;
        x [i+4] = y [i+4] ;
    }
}

/* =========================================================================
   ======================== sopt_copyi =====================================
   =========================================================================
   Copy SOPTINT vector y into vector x
   ========================================================================= */
void sopt_copyi
(
    SOPTINT       *x, /* output of copy */
    SOPTINT const *y, /* input of copy */
    SOPTINT const  n  /* length of vectors */
)
{
    SOPTINT i, n5 ;
    if ( (x == NULL) || (y == x) ) return ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++) x [i] = y [i] ;
    for (; i < n; i += 5)
    {
        x [i]   = y [i] ;
        x [i+1] = y [i+1] ;
        x [i+2] = y [i+2] ;
        x [i+3] = y [i+3] ;
        x [i+4] = y [i+4] ;
    }
}

/* =========================================================================
   ======================== sopt_copy_int ==================================
   =========================================================================
   Copy int vector y into vector x
   ========================================================================= */
void sopt_copy_int
(
    int           *x, /* output of copy */
    int     const *y, /* input of copy */
    SOPTINT const  n  /* length of vectors */
)
{
    SOPTINT i, n5 ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++) x [i] = y [i] ;
    for (; i < n; i += 5)
    {
        x [i]   = y [i] ;
        x [i+1] = y [i+1] ;
        x [i+2] = y [i+2] ;
        x [i+3] = y [i+3] ;
        x [i+4] = y [i+4] ;
    }
}

/*  =========================================================================
    ==== sopt_dot ===========================================================
    =========================================================================
    Compute dot product of x and y
    ========================================================================= */
SOPTFLOAT sopt_dot
(
    SOPTFLOAT const *x, /* first vector */
    SOPTFLOAT const *y, /* second vector */
    SOPTINT   const  n  /* length of vectors */
)
{
#ifndef NOBLAS
    if ( n >= DDOT_START )
    {
        BLAS_INT int_one = 1 ;
        BLAS_INT N = (BLAS_INT) n ;
        return (SOPT_DDOT (&N, (SOPTFLOAT *) x, &int_one,
                               (SOPTFLOAT *) y, &int_one)) ;
    }
    else
#endif
    {
        SOPTINT i, n5 ;
        SOPTFLOAT t = SOPTZERO ;
        if ( n <= 0 ) return (t) ;
        n5 = n % 5 ;
        for (i = 0; i < n5; i++) t += x [i]*y [i] ;
        for (; i < n; i += 5)
        {
            t += x [i]*y [i] + x [i+1]*y [i+1] + x [i+2]*y [i+2] 
                             + x [i+3]*y [i+3] + x [i+4]*y [i+4] ;
        }
        return (t) ;
    }
}

/* =========================================================================
   ===================== sopt_initx ========================================
   =========================================================================
   Initialize a SOPTFLOAT array
   ========================================================================= */
void sopt_initx
(
    SOPTFLOAT      *x,  /* array to be initialized */
    SOPTFLOAT const s,  /* scalar */
    SOPTINT   const n   /* length of x */
)
{
    SOPTINT j, n5 ;
    SOPTFLOAT *xj ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] = s ;
    xj = x+j ;
    for (; j < n; j += 5)
    {
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
    }
}

/* =========================================================================
   ===================== sopt_initi ========================================
   =========================================================================
   Initialize a SOPTINT array
   ========================================================================= */
void sopt_initi
(
    SOPTINT      *x,  /* array to be initialized */
    SOPTINT const s,  /* scalar */
    SOPTINT const n   /* length of x */
)
{
    SOPTINT j, n5 ;
    SOPTINT *xj ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] = s ;
    xj = x+j ;
    for (; j < n; j += 5)
    {
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
    }
}

/* =========================================================================
   ===================== sopt_init_int =====================================
   =========================================================================
   Initialize an int array
   ========================================================================= */
void sopt_init_int
(
    int          *x,  /* array to be initialized */
    int     const s,  /* scalar */
    SOPTINT const n   /* length of x */
)
{
    SOPTINT j, n5 ;
    int     *xj ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] = s ;
    xj = x+j ;
    for (; j < n; j += 5)
    {
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
    }
}

/* =========================================================================
   ==== sopt_sup_normx =====================================================
   =========================================================================
   Compute sup norm of a SOPTFLOAT vector
   ========================================================================= */
SOPTFLOAT sopt_sup_normx
(
    SOPTFLOAT const *x, /* vector */
    SOPTINT   const  n  /* length of vector */
)
{
#ifndef NOBLAS
    if ( n >= IDAMAX_START )
    {
        BLAS_INT int_one = 1 ;
        BLAS_INT N = (BLAS_INT) n ;
        BLAS_INT i = SOPT_IDAMAX (&N, (SOPTFLOAT *) x, &int_one) ;
        return (fabs (x [i-1])) ; /* adjust for fortran indexing */
    }
    else
#endif
    {
        SOPTINT i, n5 ;
        SOPTFLOAT t ;
        t = SOPTZERO ;
        n5 = n % 5 ;

        for (i = 0; i < n5; i++) if ( t < fabs (x [i]) ) t = fabs (x [i]) ;
        for (; i < n; i += 5)
        {
            if ( t < fabs (x [i]  ) ) t = fabs (x [i]  ) ;
            if ( t < fabs (x [i+1]) ) t = fabs (x [i+1]) ;
            if ( t < fabs (x [i+2]) ) t = fabs (x [i+2]) ;
            if ( t < fabs (x [i+3]) ) t = fabs (x [i+3]) ;
            if ( t < fabs (x [i+4]) ) t = fabs (x [i+4]) ;
        }
        return (t) ;
    }
}

/* =========================================================================
   === sopt_supi ===========================================================
   =========================================================================
   Return largest entry of an SOPTINT vector
   ========================================================================= */
SOPTINT sopt_supi
(
    SOPTINT const *x, /* vector */
    SOPTINT const  n  /* length of vector */
)
{
    SOPTINT xsup ;
    SOPTINT j, n5 ;
    n5 = n % 5 ;             /* n5 = n mod 5 */
    xsup = -SuiteOPTinfint ; /* initializing xsup */
    for (j = 0; j < n5; j++)
    {
        if ( xsup < x [j] ) xsup = x [j] ;
    }
    for (; j < n; j += 5)
    {
        if ( xsup < x [j]   ) xsup = x [j] ;
        if ( xsup < x [j+1] ) xsup = x [j+1] ;
        if ( xsup < x [j+2] ) xsup = x [j+2] ;
        if ( xsup < x [j+3] ) xsup = x [j+3] ;
        if ( xsup < x [j+4] ) xsup = x [j+4] ;
    }
    return (xsup) ;
}

/* ========================================================================== */
/* === sopt_transpose ======================================================= */
/* ========================================================================== */
/*    Transpose a sparse matrix: B = A' */
/* ========================================================================== */
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
)
{
    SOPTINT i, j, p, q, pp ;

    /* ====================================================================== */
    /* === compute row counts of A ========================================== */
    /* ====================================================================== */

    for (i = 0; i < nrow; i++)
    {
        W [i] = 0 ;
    }
    p = 0 ;
    for (j = 1; j <= ncol; j++)
    {
        q = Ap [j] ;
        for (; p < q; p++)
        {
            W [Ai [p]]++ ;
        }
    }

    /* ====================================================================== */
    /* === compute column pointers of B given the counts ==================== */
    /* ====================================================================== */

    Bp [0] = 0 ;
    for (i = 0; i < nrow; i++)
    {
        Bp [i+1] = Bp [i] + W [i] ;
        W [i] = Bp [i] ;
    }

    /* ====================================================================== */
    /* === B = A' =========================================================== */
    /* ====================================================================== */

    p = 0 ;
    for (j = 0 ; j < ncol ; j++)
    {
        q = Ap [j+1] ;
        for (; p < q; p++)
        {
            pp = W [Ai [p]]++ ;
            Bi [pp] = j ;
            Bx [pp] = Ax [p] ;
        }
    }
}


/* ========================================================================== */
/* === sopt_sort_cols ======================================================= */
/* ========================================================================== */
/*    Perform double transpose to sort row indices in each column of matrix   */
/* ========================================================================== */
int sopt_sort_cols /* returned integer:
                         0 if conversion successful
                         SOPT_OUT_OF_MEMORY */
(
    SOPTINT    *Ap, /* column pointers */
    SOPTINT    *Ai, /* row indices */
    SOPTFLOAT  *Ax, /* numerical values */
    SOPTINT   *atp, /* row pointers for transpose */
    SOPTINT   *ati, /* column indices for transpose */
    SOPTFLOAT *atx, /* numerical values for transpose */
    SOPTINT   nrow, /* number of rows */
    SOPTINT   ncol  /* number of cols */
)
{
    int status ;
    SOPTINT *ATp, *ATi, *temp ;
    SOPTFLOAT *ATx ;
    status = 0 ;
    if ( atp == NULL )
    {
        ATp = (SOPTINT *)   sopt_malloc (&status, (nrow+1), sizeof (SOPTINT)) ;
        ATi = (SOPTINT *)   sopt_malloc (&status, Ap [ncol],sizeof (SOPTINT)) ;
        ATx = (SOPTFLOAT *) sopt_malloc (&status, Ap [ncol],sizeof (SOPTFLOAT));
    }
    else
    {
        ATp = atp ;
        ATi = ati ;
        ATx = atx ; 
    }
    temp = (SOPTINT *) sopt_malloc
                             (&status, SOPTMAX (nrow, ncol), sizeof (SOPTINT)) ;
    if ( status == SOPT_OUT_OF_MEMORY ) return (status) ;
    sopt_transpose (ATp, ATi, ATx, Ap, Ai, Ax, nrow, ncol, temp) ;
    sopt_transpose (Ap, Ai, Ax, ATp, ATi, ATx, ncol, nrow, temp) ;
    if ( atp == NULL )
    {
        sopt_free (ATp) ;
        sopt_free (ATi) ;
        sopt_free (ATx) ;
    }
    sopt_free (temp) ;
    return (0) ;
}

/* =========================================================================
   === sopt_convert_dense_to_sparse ========================================
   =========================================================================
   convert a dense matrix into a sparse matrix
   Sparse matrices are stored using 3 arrays:

        Ai - row indices of the nonzeros elements. Indices for each
             column should be in increasing order.
        Ax - the nonzeros numerical entries corresponding to elements of Ai
        Ap - the array of column pointers of size ncol + 1.  Ap [j] is
             location of first nonzero in Ax associated with column j and
             Ap [ncol] is the number of nonzeros in the matrix
   ========================================================================= */
int sopt_convert_dense_to_sparse /* returned integer:
                                  0 if conversion successful
                                  SOPT_OUT_OF_MEMORY
                                  SOPT_ERROR_IN_INPUT_MATRIX */
(
    SOPTINT            **Ap, /* return column pointers */
    SOPTINT            **Ai, /* return row indices */
    SOPTFLOAT          **Ax, /* return numerical entries */
    SOPTFLOAT const      *A, /* matrix entries in dense format */
    SOPTINT   const    nrow, /* number of rows in A */
    SOPTINT   const    ncol, /* number of cols in A */
    int       const by_rows  /* T => matrix stored by rows, F => by columns */
)
{
    int status ;
    SOPTINT i, j, k, nnz, *ap, *ap_shift, *ai ;
    SOPTFLOAT *ax ;
    if ( (nrow == EMPTY) || (ncol == EMPTY) )
    {
        printf ("\nWhen converting a dense matrix to a sparse matrix, it\n"
                "was found that either the number of rows or the number\n"
                "of columns was not provided.\n") ;
        return (SOPT_ERROR_IN_INPUT_MATRIX) ;
    }
    /* count the number of nonzeros in each column */
    *Ap = ap = (SOPTINT *) sopt_malloc (&status, (ncol+1), sizeof (SOPTINT)) ;
    if ( status == SOPT_OUT_OF_MEMORY ) return (status) ;
    sopt_initi (ap, (SOPTINT) 0, ncol+1) ;
    k = 0 ;
    if ( by_rows )
    {
        for (i = 0; i < nrow; i++)
        {
            SOPTINT const l = k + ncol ;
            ap_shift = ap-k ;
            for (; k < l; k++)
            {
                if ( A [k] != SOPTZERO )
                {
                    ap_shift [k]++ ;
                }
            }
        }
        /* ap now contains the number of nonzeros in each column of A
           change ap to pointers to the start of each column */
        k = 0 ;
        for (j = 0; j < ncol; j++)
        {
            SOPTINT const l = k + ap [j] ;
            ap [j] = k ;
            k = l ;
        }
        ap [ncol] = k ;
        if ( k == 0 )
        {
            *Ap = *Ai = NULL ; *Ax = NULL ;
            return (0) ; /* the matrix is all zero */
        }
        ai = (SOPTINT *)   sopt_malloc (&status, k, sizeof (SOPTINT)) ;
        ax = (SOPTFLOAT *) sopt_malloc (&status, k, sizeof (SOPTFLOAT)) ;
        if ( status == SOPT_OUT_OF_MEMORY ) return (status) ;
        k = 0 ;
        for (i = 0; i < nrow; i++)
        {
            SOPTINT const l = k + ncol ;
            ap_shift = ap-k ;
            for (; k < l; k++)
            {
                if ( A [k] != SOPTZERO )
                {
                    j = ap_shift [k]++ ;
                    ai [j] = i ;
                    ax [j] = A [k] ;
                }
            }
        }
        for (j = ncol-1; j > 0; j--)
        {
            ap [j] = ap [j-1] ;
        }
        ap [0] = 0 ;
    }
    else /* matrix is stored by columns */
    {
        /* When a dense matrix stored by columns is converted to a sparse
           matrix by the program sopt_convert_to_sparse, the leading
           dimension must be >= the number of rows in the matrix. */
        ap [0] = 0 ;
        nnz = 0 ;
        for (j = 0; j < ncol; j++)
        {
            SOPTINT const l = k + nrow ;
            for (; k < l; k++)
            {
                if ( A [k] != SOPTZERO ) nnz++ ;
            }
            ap [j+1] = nnz ;
        }
        ai = (SOPTINT *)   sopt_malloc (&status, nnz, sizeof (SOPTINT)) ;
        ax = (SOPTFLOAT *) sopt_malloc (&status, nnz, sizeof (SOPTFLOAT)) ;
        if ( status == SOPT_OUT_OF_MEMORY ) return (status) ;
        nnz = k = 0 ;
        for (j = 0; j < ncol; j++)
        {
            SOPTINT const k0 = k ;
            SOPTINT const l = k + nrow ;
            for (; k < l; k++)
            {
                if ( A [k] != SOPTZERO )
                {
                    ai [nnz] = k - k0 ;
                    ax [nnz] = A [k] ;
                    nnz++ ;
                }
            }
        }
    }
    *Ai = ai ;
    *Ax = ax ;
    return (0) ;
}
/* =========================================================================
   === sopt_convert_triple_to_sparse =======================================
   =========================================================================
   convert triples Ti (row indices), Tj (column indices), Tx (numerical values)
   to              Ai (row indices), Ap (column pointer), Ax (numerical values)

   If either the nrow or the ncol argument are equal to EMPTY, then there
   values are set to the 1 + max(Ti) or max(Tj) respectively.

   Sparse matrices are stored using 3 arrays:

        Ai - row indices of the nonzeros elements. Indices for each
             column should be in increasing order.
        Ax - the nonzeros numerical entries corresponding to elements of Ai
        Ap - the array of column pointers of size ncol + 1.  Ap [j] is
             location of first nonzero in Ax associated with column j and
             Ap [ncol] is the number of nonzeros in the matrix
   ========================================================================= */
int sopt_convert_triple_to_sparse /* returned integer:
                                   0 if conversion successful
                                   SOPT_OUT_OF_MEMORY
                                   SOPT_ERROR_IN_INPUT_MATRIX
                                   SOPT_MATRIX_ELEMENT_WAS_ZERO */
(
    SOPTINT         **Ap, /* column pointers */
    SOPTINT         **Ai, /* row indices (increasing order in each column) */
    SOPTFLOAT       **Ax, /* numerical entries */
    SOPTINT        *nrow, /* return 1 + largest element of Ti if EMPTY */
    SOPTINT        *ncol, /* return 1 + largest element of Tj if EMPTY */
    SOPTINT   const  *Ti, /* row    indices   of nonzero entries */
    SOPTINT   const  *Tj, /* column indices   of nonzero entries */
    SOPTFLOAT const  *Tx, /* nonzero numerical values */
    SOPTINT   const  nnz, /* number of nonzeros in the matrix */
    int const order_cols, /* TRUE if row indices in each column of the sparse
                             matrix format should be put increasing order */
    int const        sym  /* TRUE => matrix symmetric and only elements on the
                                     diagonal and one side of the diagonal given
                             FALSE => all the matrix elements are given */
                                    
)
{
    int status ;
    SOPTINT imax, i, j, k, l, Ncol, *ap, *ai ;
    SOPTFLOAT *ax ;

    /* set default status */
    status = 0 ;

    /* immediately return if no matrix is input */
    if ( (*nrow == 0) || (*ncol == 0) || (nnz == 0) ||
         (Ti == NULL) || (Tj == NULL) || (Tx == NULL) )
    {
        *Ap = *Ai = NULL ; *Ax = NULL ;
        return (status) ;
    }

    if ( nnz < 0 )
    {
        printf ("\nWhen converting a matrix in triples format to a sparse\n"
                "matrix, it was found that the number of nonzeros Tnz for\n"
                "the triples was negative.\n") ;
        return (SOPT_ERROR_IN_INPUT_MATRIX) ;
    }
    Ncol = *ncol ;
    if ( Ncol < 0 )
    {
        Ncol = 1 + sopt_supi (Tj, nnz) ;
        if ( sym )
        {
            i = 1 + sopt_supi (Ti, nnz) ;
            if ( i > Ncol ) Ncol = i ;
        }
    }

    /* count the number of nonzeros in each column */
    *Ap = ap = (SOPTINT *) sopt_malloc (&status, Ncol+1, sizeof (SOPTINT)) ;
    if ( status == SOPT_OUT_OF_MEMORY ) return (status) ;

    /* count number of nonzeros in each column */
    k = 0 ;
    imax = -1 ;
    sopt_initi (ap, (SOPTINT) 0, Ncol) ;
    for (l = 0; l < nnz; l++)
    {
        j = Tj [l] ;
        if ( j < 0 )
        {
            printf ("\nWhen converting a matrix from triples format to sparse\n"
                    "format, a column index was negative.\n") ;
            return (SOPT_ERROR_IN_INPUT_MATRIX) ;
        }
        else if ( j >= Ncol )
        {
            printf ("\nWhen converting a matrix from triples format to sparse\n"
                    "format, a column index %ld exceeded the number of\n"
                    "columns of the input matrix. NOTE: The software employs\n"
                    "C indexing where the first row and column is 0\n",
                     (LONG) j) ;
            return (SOPT_ERROR_IN_INPUT_MATRIX) ;
        }
        ap [j]++ ;
        i = Ti [l] ;
        if ( i < 0 )
        {
            printf ("\nWhen converting a matrix from triples format to sparse\n"
                    "format, a row index was negative.\n") ;
            return (SOPT_ERROR_IN_INPUT_MATRIX) ;
        }
        if ( i > imax ) imax = i ;
        if ( sym && (i != j) )
        {
            ap [i]++ ;
            if ( j > imax ) imax = j ;
        }
        if ( Tx [l] == SOPTZERO )
        {
            ap [j]-- ;
            if ( sym && (i != j) ) ap [i]-- ;
            status = SOPT_MATRIX_ELEMENT_WAS_ZERO ;
        }
    }
    if ( status == SOPT_MATRIX_ELEMENT_WAS_ZERO )
    {
        printf ("When converting a matrix from triples format to sparse\n"
                "matrix format, it was found that one or more matrix\n"
                "elements in triples format were zero. These zeros were\n"
                "ignored\n") ;
    }
    imax++ ;
    if ( *nrow < 0 ) *nrow = imax ;
    else
    {
        if ( *nrow < imax )
        {
            printf ("When converting a matrix from triples format to a sparse\n"
                    "matrix, it was found that the number of rows = %ld was\n"
                    "less than the number of rows = %ld of the input matrix.\n",
                    (LONG) *nrow, (LONG) imax) ;
            return (SOPT_ERROR_IN_INPUT_MATRIX) ;
        }
    }

    /* ap now contains the number of nonzeros in each column of A;
       change ap to pointers to start of each column */
    k = 0 ;
    for (j = 0; j < Ncol; j++)
    {
        l = k + ap [j] ;
    /* flag column containing matrix elements: ap [j] <= 0 means that
       column j contains matrix elements, but none have been put there yet */
        if ( l > k ) ap [j] = -(k+1) ;
        else         ap [j] =  k ;
        k = l ;
    }
    ap [Ncol] = k ;
    if ( k == 0 ) return (0) ; /* the matrix is all zero */
    *Ai = ai = (SOPTINT *)   sopt_malloc (&status, k, sizeof (SOPTINT)) ;
    *Ax = ax = (SOPTFLOAT *) sopt_malloc (&status, k, sizeof (SOPTFLOAT)) ;
    if ( status == SOPT_OUT_OF_MEMORY ) return (status) ;

    int col_out_of_order = 0 ; /* changes to 1 when column found with rows
                                  out-of-order */
    /* store entries in Ai and Ax */
    for (k = 0; k < nnz; k++)
    {
        i = Ti [k] ;
        j = Tj [k] ;
        l = ap [j] ;
        SOPTFLOAT Txk = Tx [k] ;
        if ( !Txk ) continue ; /* skip zero matrix elements */
        if ( l < 0 ) /* no elements have been put in column so far */
        {
            l = -(l+1) ;
            ai [l] = i ;
            ap [j] = l + 1 ;
        }
        else
        {
            ap [j]++ ;
            ai [l] = i ;
            if ( order_cols )
            {
                if ( i <= ai [l-1] )
                {
                    col_out_of_order = 1 ;
                }
            }
        }
        ax [l] = Txk ;
        if ( sym && (i != j) ) /* if symmetric matrix, swap i and j if i != j */
        {
            l = ap [i] ;
            if ( l < 0 ) /* no elements have been put in column so far */
            {
                l = -(l+1) ;
                ai [l] = j ;
                ap [i] = l + 1 ;
            }
            else
            {
                ap [i]++ ;
                ai [l] = j ;
                if ( order_cols )
                {
                    if ( j <= ai [l-1] )
                    {
                        col_out_of_order = 1 ;
                    }
                }
            }
            ax [l] = Txk ;
        }
    }

    /* shift ap array down */
    for (j = Ncol-1; j > 0; j--)
    {
        ap [j] = ap [j-1] ;
    }
    ap [0] = 0 ;

    /* perform double transpose to sort the columns so that row indices
       are in increasing order when the columns are not yet sorted */
    if ( order_cols && col_out_of_order )
    {
        status = sopt_sort_cols (ap, ai, ax, NULL, NULL, NULL, *nrow, Ncol) ;
    }
    *ncol = Ncol ;
    return (status) ;
}

/* ==========================================================================
   === sopt_check_matrix ====================================================
   ==========================================================================
    Check a sparse matrix for consistency. The following checks are
    implemented: that the row indices in each column are strictly increasing,
    that Ap [0] = 0 and Ap [j] <= Ap [j+1], and that the elements
    in Ax are all nonzero.
   ========================================================================== */
int sopt_check_matrix /* returned integer:
                         0 if no error was detected
                         SOPT_ERROR_IN_INPUT_MATRIX */
(
    SOPTINT   const  *Ap, /* column pointers */
    SOPTINT   const  *Ai, /* row indices */
    SOPTFLOAT const  *Ax, /* numerical entries */
    SOPTINT   const ncol  /* number of columns in matrix */
)
{
    SOPTINT j, p ;

    /* if there is no matrix, then return immediately */
    if ( (Ap == NULL) || (Ai == NULL) || (Ax == NULL) || (ncol == 0) )
    {
        return (0) ;
    }
    /* check that row indices in each column are in increasing order */
    if ( Ap [0] != 0 )
    {
        printf ("in check_matrix: Ap [0] != 0\n") ;
        return (SOPT_ERROR_IN_INPUT_MATRIX) ;
    }
    p = 0 ;
    for (j = 0; j < ncol; j++)
    {
        SOPTINT const q = Ap [j+1] ;
        if ( p > q )
        {
            printf ("\nIn check_matrix: Ap [%ld] > Ap [%ld]\n",
                     (LONG) j, (LONG) j+1) ;
        }
        if ( q > p )
        {
            SOPTINT oldAi = Ai [p] ;
            for (p++; p < q; p++)
            {
                if ( Ai [p] <= oldAi )
                {
                    printf ("\nIn check_matrix, the row indices for\n"
                    "column %ld are not strictly increasing. Either a\n"
                    "row index repeats or the row indices for the column\n"
                    "are not sorted in increasing order.\n\n", (LONG) j) ;
                    return (SOPT_ERROR_IN_INPUT_MATRIX) ;
                }
                oldAi = Ai [p] ;
            }
        }
    }
    SOPTINT const nnz = Ap [ncol] ;
    for (j = 0; j < nnz; j++)
    {
        if ( Ax [j] != SOPTZERO )
        {
            printf ("\nIn check_matrix, a matrix element in column %ld"
                    "was zero.\n", (LONG) j) ;
            return (SOPT_ERROR_IN_INPUT_MATRIX) ;
        }
    }
    return (0) ;
}

/* ====================================================================== */
/* === sopt_minsortx ==================================================== */
/* ======================================================================
       ________________________________________________________
      |                                                        |
      |       sort a SOPTFLOAT array in increasing order       |
      |                                                        |
      |    input:                                              |
      |                                                        |
      |         x     --SOPTFLOAT array of numbers of length n |
      |         w     --SOPTINT working array of length n      |
      |         n     --number of array elements to sort       |
      |                                                        |
      |    output:                                             |
      |                                                        |
      |         x     --original array (SOPTFLOAT *) length n  |
      |         y     --indices of x giving increasing order   |
      |                 (SOPTINT *) of length n                |
      |________________________________________________________| */

void sopt_minsortx
(
    SOPTINT         *y, /* n-by-1 (output) */
    SOPTFLOAT const *x, /* n-by-1 (input not modified) */
    SOPTINT         *w, /* n-by-1, (input, working array) */
    SOPTINT          n  /* number of elements to sort */
)
{
    SOPTINT *yi, *wi, i, j, k, l, m, p, q ;
    SOPTFLOAT s, t ;

    y [0] = 0 ;
    if ( n < 2 ) return ;
    if ( n < 3 )
    {
        if ( x [0] > x [1] )
        {
            y [0] = 1 ;
            y [1] = 0 ;
        }
        else y [1] = 1 ;
        return ;
    }

    j = k = 0 ;
    for (i = 1; i < n; i++)
    {
        if ( x [i] < x [j] )
        {
            w [k] = i ;
            k = i ;
        }
        y [i] = j = i ;
    }

    w [k] = n ;
    while ( k > 0 )
    {
        l = m = 0 ;
        while ( l < n )
        {
            i = l ;
            p = y [i] ;
            s = x [p] ;
            j = w [i] ;
            k = j ;
            if ( j == n )
            {
                y [i] = j ;
                l = j ;
                w [m] = p ;
                k += m - i ;
                yi = y+(i-m) ;
                for (m++; m < k; m++) w [m] = yi [m] ;
            }
            else
            {
                q = y [j] ;
                t = x [q] ;
                l = w [j] ;
                y [i] = l ;
                while ( 1 )
                {
                    if ( s > t )
                    {
                        w [m] = q ;
                        m++ ;
                        j++ ;
                        if ( j == l )
                        {
                            w [m] = p ;
                            k += m - i ;
                            yi = y+(i-m) ;
                            for (m++; m < k; m++) w [m] = yi [m] ;
                            break ;
                        }
                        q = y [j] ;
                        t = x [q] ;
                    }
                    else
                    {
                        w [m] = p ;
                        m++ ;
                        i++ ;
                        if ( i == k )
                        {
                            w [m] = q ;
                            k = m + l - j ;
                            yi = y+(j-m) ;
                            for (m++; m < k; m++) w [m] = yi [m] ;
                            break ;
                        }
                        p = y [i] ;
                        s = x [p] ;
                    }
                }
            }
        }
        if ( y [0] == n )
        {
            for (i = 0; i < n; i++) y [i] = w [i] ;
            return ;
        }

        l = m = 0 ;
        while ( l < n )
        {
            i = l ;
            p = w [i] ;
            s = x [p] ;
            j = y [i] ;
            k = j ;
            if ( j == n )
            {
                w [i] = j ;
                l = j ;
                y [m] = p ;
                k += m - i ;
                wi = w+(i-m) ;
                for (m++; m < k; m++) y [m] = wi [m] ;
            }
            else
            {
                q = w [j] ;
                t = x [q] ;
                l = y [j] ;
                w [i] = l ;
                while ( 1 )
                {
                    if ( s > t )
                    {
                        y [m] = q ;
                        m++ ;
                        j++ ;
                        if ( j == l )
                        {
                            y [m] = p ;
                            k += m - i ;
                            wi = w+(i-m) ;
                            for (m++; m < k; m++) y [m] = wi [m] ;
                            break ;
                        }
                        q = w [j] ;
                        t = x [q] ;
                    }
                    else
                    {
                        y [m] = p ;
                        m++ ;
                        i++ ;
                        if ( i == k )
                        {
                            y [m] = q ;
                            k = m + l - j ;
                            wi = w+(j-m) ;
                            for (m++; m < k; m++) y [m] = wi [m] ;
                            break ;
                        }
                        p = w [i] ;
                        s = x [p] ;
                    }
                }
            }
        }
        if ( y [0] == n ) return ;
    }
}

/* ====================================================================== */
/* === sopt_minsorti ==================================================== */
/* ======================================================================
       ________________________________________________________
      |                                                        |
      |       sort an SOPTINT array in increasing order        |
      |                                                        |
      |    input:                                              |
      |                                                        |
      |         x  --array of numbers (SOPTINT *) of length n  |
      |         w  --working array (SOPTINT *) of length n     |
      |         n  --number of array elements to sort          |
      |                                                        |
      |    output:                                             |
      |                                                        |
      |         x  --original array (SOPTINT *) of length n    |
      |         y  --indices of x giving increasing order      |
      |              (SOPTINT *) of length n                   |
      |________________________________________________________| */

void sopt_minsorti
(
    SOPTINT        *y, /* n-by-1 (output) */
    SOPTINT  const *x, /* n-by-1 (input not modified) */
    SOPTINT        *w, /* n-by-1, (input, working array) */
    SOPTINT         n  /* number of elements to sort */
)
{
    SOPTINT *yi, *wi, i, j, k, l, m, p, q ;
    SOPTINT s, t ;

    y [0] = 0 ;
    if ( n < 2 ) return ;

    j = k = 0 ;
    for (i = 1; i < n; i++)
    {
        if ( x [i] < x [j] )
        {
            w [k] = i ;
            k = i ;
        }
        y [i] = j = i ;
    }

    w [k] = n ;
    while ( k > 0 )
    {
        l = m = 0 ;
        while ( l < n )
        {
            i = l ;
            p = y [i] ;
            s = x [p] ;
            j = w [i] ;
            k = j ;
            if ( j == n )
            {
                y [i] = j ;
                l = j ;
                w [m] = p ;
                k += m - i ;
                yi = y+(i-m) ;
                for (m++; m < k; m++)
                {
                    w [m] = yi [m] ;
                }
            }
            else
            {
                q = y [j] ;
                t = x [q] ;
                l = w [j] ;
                y [i] = l ;
                while ( 1 )
                {
                    if ( s > t )
                    {
                        w [m] = q ;
                        m++ ;
                        j++ ;
                        if ( j == l )
                        {
                            w [m] = p ;
                            k += m - i ;
                            yi = y+(i-m) ;
                            for (m++; m < k; m++)
                            {
                                w [m] = yi [m] ;
                            }
                            break ;
                        }
                        q = y [j] ;
                        t = x [q] ;
                    }
                    else
                    {
                        w [m] = p ;
                        m++ ;
                        i++ ;
                        if ( i == k )
                        {
                            w [m] = q ;
                            k = m + l - j ;
                            yi = y+(j-m) ;
                            for (m++; m < k; m++)
                            {
                                w [m] = yi [m] ;
                            }
                            break ;
                        }
                        p = y [i] ;
                        s = x [p] ;
                    }
                }
            }
        }
        if ( y [0] == n )
        {
            for (i = 0; i < n; i++)
            {
                y [i] = w [i] ;
            }
            return ;
        }

        l = m = 0 ;
        while ( l < n )
        {
            i = l ;
            p = w [i] ;
            s = x [p] ;
            j = y [i] ;
            k = j ;
            if ( j == n )
            {
                w [i] = j ;
                l = j ;
                y [m] = p ;
                k += m - i ;
                wi = w+(i-m) ;
                for (m++; m < k; m++)
                {
                    y [m] = wi [m] ;
                }
            }
            else
            {
                q = w [j] ;
                t = x [q] ;
                l = y [j] ;
                w [i] = l ;
                while ( 1 )
                {
                    if ( s > t )
                    {
                        y [m] = q ;
                        m++ ;
                        j++ ;
                        if ( j == l )
                        {
                            y [m] = p ;
                            k += m - i ;
                            wi = w+(i-m) ;
                            for (m++; m < k; m++)
                            {
                                y [m] = wi [m] ;
                            }
                            break ;
                        }
                        q = w [j] ;
                        t = x [q] ;
                    }
                    else
                    {
                        y [m] = p ;
                        m++ ;
                        i++ ;
                        if ( i == k )
                        {
                            y [m] = q ;
                            k = m + l - j ;
                            wi = w+(j-m) ;
                            for (m++; m < k; m++)
                            {
                                y [m] = wi [m] ;
                            }
                            break ;
                        }
                        p = w [i] ;
                        s = x [p] ;
                    }
                }
            }
        }
        if ( y [0] == n ) return ;
    }
}
