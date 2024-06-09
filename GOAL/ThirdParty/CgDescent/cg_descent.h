#ifndef _CG_DESCENT_H_
#define _CG_DESCENT_H_

#include "sopt.h"

#define CGZERO ((CGFLOAT) 0)
#define CGONE  ((CGFLOAT) 1)

/* define functions, routines, & constants in CG using sopt.h definitions */
#define CGFLOAT SOPTFLOAT
#define CGINT SOPTINT
/* infinite float */
#define CGINF SuiteOPTinf
/* infinite integer */
#define CGINFINT SuiteOPTinfint
#define CGMAX SOPTMAX
#define CGMIN SOPTMIN
#define FALSE SuiteOPTfalse
#define TRUE SuiteOPTtrue
#define cg_error sopt_error
#define cg_print_TF sopt_print_TF
#define cg_scale0 sopt_scale
#define cg_scale sopt_scale
#define cg_step sopt_step
#define cg_daxpy0 sopt_daxpy
#define cg_daxpy sopt_daxpy
#define cg_initx sopt_initx
#define cg_initi sopt_initi
#define cg_sup_normx sopt_sup_normx
#define cg_dot0 sopt_dot
#define cg_dot sopt_dot
#define cg_copy0 sopt_copyx
#define cg_copy sopt_copyx
#define cg_sort_cols sopt_sort_cols
#define cg_timer sopt_timer
#define pasa_cg_timer sopt_timer
#define cg_convert_triple_to_sparse sopt_convert_triple_to_sparse
#define cg_convert_dense_to_sparse sopt_convert_dense_to_sparse
#define cg_malloc sopt_malloc
#define cg_free sopt_free

#ifdef PASA
#define PASA_CG_COM pasacom,cgcom
#define XXCG(name) pasa_cg_ ## name
#else
#define PASA_CG_COM cgcom
#define XXCG(name) cg_ ## name
#endif

/* ==========================================================================
 * =========== status returned by CG_DESCENT ================================
 * ========================================================================== */

#define CG_ERROR_TOLERANCE_SATISFIED                             (0)
#define CG_ITERATIONS_EXCEED_MAXITS                              (200)
#define CG_SLOPE_ALWAYS_NEGATIVE                                 (201)
#define CG_LINE_SEARCH_STEPS_EXCEED_MAXSTEPS                     (202)
#define CG_SEARCH_DIRECTION_NOT_DESCENT_DIRECTION                (203)
#define CG_WOLFE_CONDITIONS_NOT_SATISFIED                        (204)
#define CG_DEBUGGER_IS_ON_AND_FUNCTION_VALUE_INCREASES           (205)
#define CG_NO_COST_OR_GRADIENT_IMPROVEMENT                       (206)
#define CG_OUT_OF_MEMORY                                         (207)
#define CG_QUADRATIC_OBJECTIVE_NO_LOWER_BOUND                    (208)
#define CG_STARTING_FUNCTION_VALUE_INFINITE_OR_NAN               (209)
#define CG_EXCESSIVE_UPDATING_OF_PERT_EPS                        (210)
#define CG_FUNCTION_NAN_OR_INF                                   (211)
#define CG_QP_LINEAR_TERM_GIVEN_BUT_HPROD_MISSING                (212)
#define CG_N_IS_EMPTY                                            (213)
#define CG_ERROR_IN_INPUT_MATRIX                                 (214)
#define CG_HESSIAN_ELEMENT_REPEATS                               (215)

#define CG_START_MESSAGES                                        (200)
#define CG_END_MESSAGES                                          (299)

/* not seen by the caller, used for flow constrol */
#define CG_OK                                                    (-1)
#define CG_WOLFE_OK                                              (-2)
#define CG_WOLFE_NOT_OK                                          (-3)
#define CG_NEW_PERT                                              (-4)
#define CG_INTERVAL_OK                                           (-5)
#define CG_ERROR_TOLERANCE_DOES_NOT_HOLD                         (-6)
#define CG_RESTART                                               (-7)
#define CG_CONTINUE                                              (-8)
#define CG_HITS_BOUNDARY                                         (-9)

/* -------------------------------------------------------------------------- */
/* cg_descent version information */
/* -------------------------------------------------------------------------- */
#define CG_DATE "February 15, 2022"
#define CG_MAIN_VERSION 8
#define CG_SUB_VERSION 0
#define CG_SUBSUB_VERSION 0


#pragma region SzxExtension
#include <functional>
using CgEval = std::function<void(CGFLOAT* fg, CGFLOAT* x, CGINT n)>;
using CgEval2 = std::function<void(CGFLOAT* f, CGFLOAT* g, CGFLOAT* x, CGINT n)>;
#pragma endregion SzxExtension



/* ==========================================================================
 * =============== cg structures ============================================
 * ========================================================================== */

/* ============================================================================
   CGparm is a structure containing parameters used in cg_descent.
   CGdefault gives default values for these parameters.
============================================================================ */
typedef struct CGparm_struct
{
    /* T => print status of run
       F => do not print status of run */
    int PrintStatus ;

    /* T => print pasa statistics
       F => do not print statistics */
    int PrintStat ;

    /* T => print parameter values 
       F => do not print parameter values */
    int PrintParm ;

    /* Level 0  = no printing, ... , Level 2 = maximum printing */
    int PrintLevel ;

    /* There are multiple ways to evaluate some of the parameters in
       cg_descent, some are faster but less stable than others. If
       FastLA = TRUE, then the fast options are always used. */
    int FastLA ;

    /* replace the Hessian Q of a quadratic by Q + QPshift*I */
    CGFLOAT QPshift ;

    /* nominal stopping criterion for cg_descent is ||grad||_infty <= grad_tol*/
    CGFLOAT grad_tol ;

    /* CG_DESCENT no longer includes the stopping criterion

         ||proj_grad||_infty <= grad_tol*(1 + |f_k|).

      If the optimization problem is unconstrained, then the stopping
      criterion is

         ||proj_grad||_infty <= testtol

      where testtol = max(grad_tol,initial ||grad||_infty*StopFact) and
      the default value of StopFact is zero. If the optimization problem
      contains constraints, testtol is the pasa stopping criterion
      switchfactor*global_error. This value for testtol is the PASA
      criterion to stop solving the unconstrained problem and return
      to the gradient project algorithm. */

    CGFLOAT StopFac ;

    /* T => check that f_k+1 - f_k <= debugtol*C_k
       F => no checking of function values */
    int debug ;
    CGFLOAT debugtol ;

    /* If the sparse Hessian for a QP is provided using the triples format,
       then T => check that there are no duplications of row indices in a
       column, while F => skip this check */
    int CheckMatrix ;

    /* if step is nonzero, it is the initial step of the initial line search
       NOTE: if the user specifies the parameter bbk in PASA for the initial
             bb step, then bbk replaces the value of step in cg */
    CGFLOAT step ;

    /* LBFGS = 0 => use cg_descent
               1 => use L-BFGS
               2 => use L-BFGS when LBFGSmemory >= n, use cg_descent if memory<n
               3 => use L-BFGS when LBFGSmemory >= n, else use limited memory CG
    */
    int LBFGS ;

    /* if LBFGS is used, then LBFGSmemory is the number of vectors in memory */
    int LBFGSmemory ;

    /* abort cg after maxit iterations */
    CGINT maxit ;

    /* conjugate gradient method restarts after (n*restart_fac) iterations */
    CGFLOAT restart_fac ;

    /* factor in [0, 1] used to compute average cost magnitude C_k as follows:
       Q_k = 1 + (Qdecay)Q_k-1, Q_0 = 0,  C_k = C_k-1 + (|f_k| - C_k-1)/Q_k */
    CGFLOAT Qdecay ;

    /* terminate after nslow iterations without strict improvement in
       either function value or gradient */
    int nslow ;

    /* factor by which eps grows when line search fails during contraction */
    CGFLOAT egrow ;

    /* T => attempt quadratic interpolation in line search when
                |f_k+1 - f_k|/f_k <= QuadCutoff
       F => no quadratic interpolation step */
    int    QuadStep ;
    CGFLOAT QuadCutOff ;

    /* maximum factor by which a quad step can reduce the step size */
    CGFLOAT QuadSafe ;

    CGFLOAT psi_lo ; /* in performing a QuadStep, we evaluate at point
                        betweeen [psi_lo, psi_hi]*psi2*previous step */
    CGFLOAT psi_hi ;
    CGFLOAT   psi1 ; /* for approximate quadratic, use gradient at
                        psi1*psi2*previous step for initial stepsize */

    CGFLOAT   qeps ; /* parameter in cost error for quadratic restart
                        criterion */
    CGFLOAT  qrule ; /* parameter used to decide if cost is quadratic */
    int   qrestart ; /* number of iterations the function should be
                        nearly quadratic before a restart */

    /* T => when possible, use a cubic step in the line search */
    int UseCubic ;

    /* use cubic step when |f_k+1 - f_k|/|f_k| > CubicCutOff */
    CGFLOAT CubicCutOff ;

    /* |f| < SmallCost*starting cost => skip QuadStep and set PertRule = FALSE*/
    CGFLOAT SmallCost ;

    /* maximum factor secant step increases stepsize in expansion phase */
    CGFLOAT ExpandSafe ;

    /* factor by which secant step is amplified during expansion phase
       where minimizer is bracketed */
    CGFLOAT SecantAmp ;

    /* factor by which rho grows during expansion phase where minimizer is
       bracketed */
    CGFLOAT RhoGrow ;

   /*  If the currect derivative <= BigDfactor * derivative at starting point,
       then we evaluate the function at the current point even though
       AvoidFeval is TRUE. AvoidFeval is TRUE when we think that
       the function values have converged; in this case, line search
       is mostly based on the derivative. */
    CGFLOAT BigDfactor ;

    /* If approxstep is TRUE, use approximate Wolfe line search.
       If approxstep is FALSE, then use ordinary line search and
       switch to the approximate step when |f - fnew| < ApproxSwitchFactor*|fR|,
       where fR is an estimate of the function size.  In the conjugate
       gradient code, an averaging technique is used to estimate function
       size. */
    int approxstep ;
    CGFLOAT ApproxSwitchFactor ;

    /* As the cg_descent converges, the function values typically
       approach a constant value. When this happens, the cubic interpolation
       step in the line search loses its accuracy and it is better to
       use a secant step based on the derivative of the objective function.
       The cost has converged when the relative change in the objective
       function <= CostConverge */
    CGFLOAT CostConverge ;

    CGFLOAT FuncGradSwitchFactor ;

    /* When performing an approximate Wolfe line search, we require
       that the new function value <= perturbation of the prior
       function value where the perturbation tries to take into
       rounding errors associated with the function value.
       There are two different perturbations:
       PertRule = 1 => fpert is f + Parm->pert_eps*|f| (relative change)
       PertRule = 0 => fpert is f + Parm->pert_eps     (absolute change) */
    int PertRule ;
    CGFLOAT pert_eps ;

    /* Maximum number of contractions in cg_contract. If it cannot find a
       step that either satisfies the Wolfe conditions or which has derivative
       >= 0 within ncontract attempts, then it is felt that fpert is too
       small and it will be increased so that the current function value
       is less than fpert. To increase fpert, we increase the value of
       pert_eps which is used to compute fpert. */
    int ncontract ;

    /* When pert_eps is increased, we multiply the new value by the growth
       factor eps_grow to ensure a healthy growth in pert_eps. */
    CGFLOAT eps_grow ;

    /* Maximum number of times that pert_eps is recomputed before a line
       search error is declared. */
    int neps ;

    CGFLOAT    cgdelta ; /* Wolfe line search parameter */
    CGFLOAT    cgsigma ; /* Wolfe line search parameter */
    int       maxsteps ; /* max number of tries to find acceptable step */
    CGFLOAT  stepdecay ; /* decay factor for bracket interval width */
    CGFLOAT        rho ; /* growth factor when searching for initial
                            bracketing interval */
    CGFLOAT       psi0 ; /* factor used in starting guess for iteration 1 */
    CGFLOAT       psi2 ; /* when starting a new cg iteration, our initial
                            guess for the line search stepsize is
                            psi2*previous step */
    CGFLOAT  BetaLower ; /* parameter connected with lower bound for beta */
    CGFLOAT      theta ; /* parameter describing the cg_descent family */
    int   AdaptiveTheta ; /* T => choose theta adaptively, F => use theta */

    /* When an infinite or nan objective value is encountered, the
       stepsize is reduced in an effort to find a finite objective
       value. infdecay is the initial decay factor that is used when an
       infinite or nan objective value is encountered, ninf_tries is
       the number of attempts we make to find a finite objective value,
       and infdecay_rate is a factor by which infdecay is multiplied
       after each attampt to find a finite objective value. */
    CGFLOAT cg_infdecay ;
    CGFLOAT cg_infdecay_rate ;
    int  cg_ninf_tries ;

/* ============ LIMITED MEMORY CG PARAMETERS ================================ */
    /* SubCheck and SubSkip control the frequency with which the subspace
       condition is checked. It it checked for SubCheck*mem iterations and
       if it is not activated, then it is skipped for Subskip*mem iterations
       and Subskip is doubled. Whenever the subspace condition is satisfied,
       SubSkip is returned to its original value. */
    int SubCheck ;
    int SubSkip ;

    /* when relative distance from current gradient to subspace <= eta0,
       enter subspace if subspace dimension = mem (eta0 = 0 means gradient
        inside subspace) */

    CGFLOAT   eta0 ; /* corresponds to eta0*eta0 in the paper */

    /* when relative distance from current gradient to subspace >= eta1,
       leave subspace (eta1 = 1 means gradient orthogonal to subspace) */
    CGFLOAT   eta1 ; /* corresponds to eta1*eta1 in the paper */

    /* when relative distance from current gradient to subspace <= eta2,
       always enter subspace (invariant space) */
    CGFLOAT   eta2 ;

    /* number of expansion steps in cg when solving bound constrained problem
       through pasa */
    int     nexpand ;

    /* factor by which the stepsize grows in each expansion step */
    CGFLOAT step_growth ;
} CGparm ;

/* --------------------------------------------------------------------------
    CGstat is a structure containing statistics which are returned to the
    user when cg_descent terminates.  They can be print by user
    using cg_print_stat or cg_print_status, or by setting the parameters
    PrintStat or PrintStatus to TRUE.
   -------------------------------------------------------------------------- */
typedef struct CGstat_struct
{
    int         status ; /* returned status from cg_descent */
    CGFLOAT          f ; /* function value at solution */
    CGFLOAT        err ; /* sup norm of the gradient */
    CGFLOAT   grad_tol ; /* error tolerance */
    CGFLOAT        tol ; /* computing tolerance: gnorm <= tol */
    CGINT        maxit ; /* maximum number of iterations */
    int  cg_ninf_tries ; /* number of tries to find finite objective value */
    CGFLOAT       oldf ; /* old function value when debugger fails */
    CGFLOAT       newf ; /* new function value when debugger fails */
    int       maxsteps ; /* max number of attempts in the line search */
    int        NegDiag ; /* T => negative diagonal encountered in QP factor */
    CGINT         iter ; /* number of iterations in cg */
    CGINT        nfunc ; /* function evaluations in cg */
    CGINT        ngrad ; /* gradient evaluations in cg */
    int        IterSub ; /* number subspace iterations in limited memory cg */
    int         NumSub ; /* number of subspaces in limited memory cg */
} CGstat ;

/* --------------------------------------------------------------------------
    CGhess is a structure that is passed to the user when cg_descent wishes
    to evaluate the Hessian of the objective.
   -------------------------------------------------------------------------- */
typedef struct CGhess_struct
{
    CGINT     ncol ; /* number of cols in the Hessian */
    CGFLOAT     *x ; /* the point where the Hessian is evaluated */
    CGINT      *Hp ; /* size ncol+1, the column pointers */
    CGINT      *Hi ; /* size Hp [ncol], the row indices */
    CGFLOAT    *Hx ; /* size Hp [ncol], numerical values in the Hessian */
} CGhess ;

typedef struct CGdata_struct
{
    /* -------- cg input data -------- */
    CGFLOAT     *x ; /* size n, points to a starting guess
                        If NULL, then cg_descent allocates memory of size n
                        for x and sets the value of its entries to zero.
                        Any allocated memory is freed by cg_terminate. */
    CGINT        n ; /* dimension of x */
    CGparm   *Parm ; /* parameter structure */
    CGstat   *Stat ; /* statistics structure */

    CGFLOAT  *Work ; /* NULL => cg_descent should allocate real work
                          space. Otherwise, see allocation section of
                          cg_descent to determine memory requirements. */

    /* for a general nonlinear function, value (f, x, n) is the function value*/
    CgEval value ;
    /* for a general nonlinear function, grad (g, x, n) is the gradient */
    CgEval grad ;
    /* for a general nonlinear function, valgrad (f, g, x, n) gives both
       function value f and gradient g at x (size n), the valgrad argument
       is optional. */
    CgEval2 valgrad ;

    /* hess evaluates the Hessian at a given point, not currently used */
    void    (*hess) (CGhess *) ;

    CGFLOAT *c ; /* linear coefficients in a QP */

    /* If the problem is a QP, then instead of providing routines
       to evaluate the objective function and its gradient, the user
       could provide the Hessian of the objective
       and any linear term in the objective (argument c above).
       If the linear term is not given, then it is taken to be zero.
       Whenever the objective value or its gradient is needed by pasa,
       it will be computed internally using the provided Hessian.
       There are three ways to input the Hessian matrix:

       1. A dense packed matrix containing the matrix entries.

       2. The nonzero elements in the matrix can be stored as a triple:
              HTi  (row    indices of nonzero matrix elements)
              HTj  (column indices of nonzero matrix elements)
              HTx  (numerical values of nonzero matrix elements)
              Hnz  (number of nonzeros in the matrix)

       3. The matrix can be stored in standard sparse matrix format:
              Hp   (column pointers, the number of entries is n + 1, Hp [j] is
                    location of first nonzero in Hx associated with column j)
              Hi   (the row indices of the nonzero matrix elements)
              Hx   (numerical values of the nonzero elements)
       When H does not exist, these elements of the data structure are NULL */

    /* dense matrix */
    CGFLOAT *Hdense ; /* size n by n, numerical values of all elements in H
                         by rows or by columns (the matrix is symmetric) */

    /* triples format */
    CGINT      *HTj ; /* size Hnz (column indices) */
    CGINT      *HTi ; /* size Hnz (row indices) */
    CGFLOAT    *HTx ; /* size Hnz (nonzero matrix elements) */
    CGINT       Hnz ; /* number of entries in HTx */
    int        Hsym ; /* TRUE  => matrix is symmetric and only matrix elements
                                  on the main diagonal and on one side are given
                         FALSE => all the matrix elements are given
                                  (the matrix could still be symmetric) */

    /* sparse matrix format */
    CGINT       *Hp ; /* column pointers for Hessian, size ncol + 1 */
    CGINT       *Hi ; /* row indices for Hessian */
    CGFLOAT     *Hx ; /* nonzero numerical values in Hessian */

    /* Alternatively, for a quadratic objective, the user could provide
       a routine for computing the product between the Hessian and a
       vector.  hprod (Hd, d, n) computes Hd = H*d where d
       has length n and the Hessian H is n by n. */
    void   (*hprod) (CGFLOAT *, CGFLOAT *, CGINT) ;

    /* the following are used internally when the code allocates memory */
    CGFLOAT *x_created ; /* pointer to memory created for x */
    int      H_created ; /* = TRUE if cg_descent creates the sparse Hessian
                            = FALSE otherwise */
} CGdata ;

typedef struct CGcom_struct /* common variables */
{
    CGdata        *cgdata ; /* pointer to cgdata structure */
    /* first two arguments below not used in unconstrained cg */
    CGFLOAT          *AAd ; /* (A'A)*d (Hessian of penalty) * d in CG */
    CGFLOAT         *gpen ; /* gradient of penalty term in CG */

    /* cg_value (f, x, n) */
    CgEval value ;
    /* cg_grad (g, x, n) */
    CgEval grad ;
    /* cg_valgrad (f, g, x, n)*/
    CgEval2 valgrad ;

    int          QuadCost ; /* TRUE if objective is quadratic */
    int      hprod_status ; /* 0 => user provided, 1 => builtin, 2 => none */
    CGFLOAT            *c ; /* linear term when objective quadratic */

    /* sparse matrix format for a QP */
    CGINT             *Hp ; /* column pointers for Hessian, size ncol + 1 */
    CGINT             *Hi ; /* row indices for Hessian */
    CGFLOAT           *Hx ; /* nonzero numerical values in Hessian */

    /* Hessian times vector for a quadratic objective. For a QP either
       the sparse matrix form or the product form is used, not both */
    void   (*hprod) (CGFLOAT *, CGFLOAT *, CGINT) ;

    CGFLOAT *work_created ; /* work array created if none provided */
    CGFLOAT          f_at ; /* stepsize where function is evaluated */
    CGFLOAT         df_at ; /* stepsize where derivative is evaluated */
    CGFLOAT           *Qd ; /* Hessian * search direction d for quadratics */
    CGFLOAT       QPshift ; /* replace Hessian Q of quadratic by Q+QPshift*I */
    CGINT               n ; /* problem dimension, saved for reference */
    CGINT       FirstIter ; /* first iteration of pasa */
    int           NegDiag ; /* F => no negative diagonal element in QR factor */
    int            QuadOK ; /* T (quadratic step successful) */
    int          UseCubic ; /* T (use cubic step) F (use secant step) */
    int          PertRule ; /* 1 => estimated error in function value is eps*f,
                               0 => estimated error in function value is eps */
    int             QuadF ; /* T => function appears to be quadratic */
    int              neps ; /* number of time eps updated */
    CGFLOAT         fpert ; /* perturbation is pert_eps*|f| if PertRule = 1 */
    CGFLOAT      pert_eps ; /* current value of pert_eps */
    CGFLOAT       maxstep ; /* used in pasa for storing maximum allowed step */
    CGFLOAT     SmallCost ; /* |f| <= SmallCost => set PertRule = F */
    CGFLOAT         alpha ; /* stepsize along search direction */
    CGFLOAT             f ; /* function value for step alpha */
    CGFLOAT            df ; /* directional derivative for step alpha */
    CGFLOAT            f0 ; /* old function value */
    CGFLOAT           df0 ; /* old derivative */
    CGFLOAT            Ck ; /* average cost as given by the rule:
                               Qk = Qdecay*Qk + 1, Ck += (fabs (f) - Ck)/Qk */
    CGFLOAT      wolfe_hi ; /* upper bound for slope in Wolfe test */
    CGFLOAT      wolfe_lo ; /* lower bound for slope in Wolfe test */
    CGFLOAT     awolfe_hi ; /* upper bound for slope, approximate Wolfe test */
    int        approxstep ; /* T (use approximate Wolfe line search)
                               F (use             Wolfe line search) */
    int        AvoidFeval ; /* T when function values converges */
    int             Wolfe ; /* becomes T when Wolfe line search performed */
    CGFLOAT      alphaold ; /* previous value for stepsize alpha */
    CGFLOAT            *x ; /* current iterate */
    CGFLOAT            *d ; /* current search direction */
    CGFLOAT            *g ; /* gradient at x */
    CGFLOAT         *xnew ; /* x + alpha*d */
    CGFLOAT         *gnew ; /* gradient at x + alpha*d */
    /* if Wolfe line search fails and we try approxstep, need the following */
    CGFLOAT        savedf ;
    CGFLOAT     savealpha ;
    CGFLOAT        savefb ;
    int            saveqb ;
    CGFLOAT      maxeqerr ; /* ||Ax-b||_inf for active constraints */
    CGparm          *Parm ; /* user parameters */
    CGstat          *Stat ; /* CG statistics */
} CGcom ;

/* ==========================================================================
   ================ prototypes ==============================================
   ========================================================================== */
int cg_descent
(
    CGdata  *cgdata /* CG data structure */
) ;

CGdata * cg_setup (void) ;

void cg_terminate
(
    CGdata **DataHandle
) ;

int cg_wrapup
(
    int       status,
    CGcom     *cgcom
) ;

int cg_evaluate
(
    CGFLOAT alpha_good, /* a value of alpha for which function is finite */
    CGFLOAT     *Alpha, /* stepsize along the search direction */
    const char   *what, /* fg = eval func & grad, g = grad only,f = func only */
    CGcom         *Com
) ;

int cg_line
(
    int       repeat, /* TRUE => Wolfe search failed, retry using approxstep */
    CGcom     *cgcom
) ;

int cg_contract
(
    CGFLOAT       *A, /* left side of bracketing interval */
    CGFLOAT      *fA, /* function value at a */
    CGFLOAT      *dA, /* derivative at a */
    CGFLOAT       *B, /* right side of bracketing interval */
    CGFLOAT      *fB, /* function value at b */
    CGFLOAT      *dB, /* derivative at b */
    CGcom     *cgcom
) ;

/* ==========================================================================
   ================ protypes for pure cg_descent routines ===================
   ========================================================================== */

void cg_default
(
    CGparm *Parm /* pointer to parameter structure */
) ;

void cg_print_status
(
    CGdata *Data /* pointer to cgdata structure */
) ;

void cg_print_stat
(
    CGdata *Data /* pointer to cgdata structure */
) ;

void cg_print_parm
(
    CGdata *Data /* pointer to cgdata structure */
) ;

int cg_Wolfe
(
    CGFLOAT  alpha, /* stepsize */
    CGFLOAT      f, /* function value associated with stepsize alpha */
    CGFLOAT   dphi, /* derivative value associated with stepsize alpha */
    CGcom   *cgcom
) ;

CGFLOAT cg_cubic
(
    CGFLOAT  a,
    CGFLOAT fa, /* function value at a */
    CGFLOAT da, /* derivative at a */
    CGFLOAT  b,
    CGFLOAT fb, /* function value at b */
    CGFLOAT db  /* derivative at b */
) ;

void cg_builtin_hprod
(
    CGFLOAT  *Hd, /* Hd = H*d */
    CGFLOAT   *d,
    CGINT      n, /* number of entries in d, H is n by n */
    CGINT    *Hp, /* column pointers for Hessian */
    CGINT    *Hi, /* row indices for Hessian */
    CGFLOAT  *Hx  /* nonzero numerical entries in Hessian */
) ;

CGFLOAT cg_update_2
(
    CGFLOAT *gold, /* old g */
    CGFLOAT *gnew, /* new g */
    CGFLOAT    *d, /* d */
    CGINT       n /* length of vectors */
) ;

void cg_update_beta
(
    CGFLOAT *oldproj,
    CGFLOAT *newproj,
    CGFLOAT   *GkPyk,
    CGFLOAT   *YkPyk,
    CGINT          n  /* length of vectors */
) ;

void cg_update_d
(
    CGFLOAT      *d,
    CGFLOAT  *gproj,
    CGFLOAT    beta,
    CGINT         n  /* length of vectors */
) ;

/* limited memory CG routines */

void cg_matvec
(
    CGFLOAT *y, /* product vector */
    CGFLOAT *A, /* dense matrix */
    CGFLOAT *x, /* input vector */
    CGINT    n, /* number of columns of A */
    CGINT    m, /* number of rows of A */
    int      w  /* T => y = A*x, F => y = A'*x */
) ;

void cg_trisolve
(
    CGFLOAT *x, /* right side on input, solution on output */
    CGFLOAT *R, /* dense matrix */
    int      m, /* leading dimension of R */
    int      n, /* dimension of triangular system */
    int      w  /* T => Rx = y, F => R'x = y */
) ;

void cg_Yk
(
    CGFLOAT    *y, /*output vector */
    CGFLOAT *gold, /* initial vector */
    CGFLOAT *gnew, /* search direction */
    CGFLOAT  *yty, /* y'y */
    CGINT       n  /* length of the vectors */
) ;

CGFLOAT cg_update_inf2
(
    CGFLOAT   *gold, /* old g */
    CGFLOAT   *gnew, /* new g */
    CGFLOAT      *d, /* d */
    CGFLOAT *gnorm2, /* 2-norm of g */
    CGINT         n  /* length of vectors */
) ;

#endif
