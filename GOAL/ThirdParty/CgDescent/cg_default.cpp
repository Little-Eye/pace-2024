#include "cg_descent.h"

/* ==========================================================================
   === cg_default =========================================================
   ==========================================================================
    Set CG default parameter values
   ========================================================================== */
void cg_default
(
    CGparm *Parm /* pointer to parameter structure */
)
{

    /* T => print status of run
       F => do not print status of run */
    Parm->PrintStatus = TRUE ;

    /* T => print pasa statistics
       F => do not print statistics */
    Parm->PrintStat = FALSE ;

    /* T => print parameter values 
       F => do not print parameter values */
    Parm->PrintParm = FALSE ;

    /* Level 0 (no printing), 1 (key branching),
             2 (more printing outside loops), 3 (also print inside loops) */
    Parm->PrintLevel = 0 ;

    /* There are multiple ways to evaluate some of the parameters in
       cg_descent, some are faster but less stable than others. If
       FastLA = TRUE, then the fast options are always used. */
    Parm->FastLA = FALSE ;

    /* replace the Hessian Q of a quadratic by Q + QPshift*I */
    Parm->QPshift = CGZERO ;

    /* nominal stopping criterion for cg_descent is ||grad||_infty <= grad_tol*/
    Parm->grad_tol = 1.e-6 ;

    /* CG_DESCENT no longer includes the stopping criterion

         ||grad||_infty <= grad_tol*(1 + |f_k|).

      The stopping criterion is now

         ||grad||_infty <= testtol

      where testtol = max(grad_tol,initial ||grad||_infty*StopFact) and
      the default value of StopFact is zero. If the optimization problem
      contains constraints and cg_descent is called by PASA, then
      testtol is the pasa stopping criterion switchfactor*global_error.
      This value for testtol is the PASA criterion to stop solving the
      unconstrained problem and return to the gradient project algorithm. */

    Parm->StopFac = CGZERO ;

    /* debug = T => check that f_k+1 - f_k <= debugtol*fR.
       debug = F => no checking of function values */
    Parm->debug = FALSE ;
    Parm->debugtol = 1.e-4 ;


    /* If the sparse Hessian for a QP is provided using the triples format,
       then T => check that there are no duplications of row indices in a
       column of the matrix, while F => skip this check */
    Parm->CheckMatrix = FALSE ;

    /* if step is nonzero, it is the initial step of the initial line search */
    Parm->step = CGZERO ;

    /* 0 => use cg_descent
       1 => use L-BFGS
       2 => use L-BFGS when LBFGSmemory >= n, use cg_descent otherwise
       3 => use L-BFGS when LBFGSmemory >= n, use limited memory CG otherwise */
    Parm->LBFGS = 3 ;

    /* if LBFGS is used, then LBFGSmemory is the number of vectors in memory */
    Parm->LBFGSmemory = 11 ;

    Parm->maxit = CGINFINT ; /* max # conjugate gradient iterations */

    /* conjugate gradient method restarts after (n*restart_fac) iterations */
    Parm->restart_fac = 6.0 ;

    /* Factor in [0, 1] that is used to estimate the average
       cost magnitude C_k.  Q_k = 1 + (Qdecay)Q_k-1, Q_0 = 0,
       C_k = C_k-1 + (|f_k| - C_k-1)/Q_k */
    Parm->Qdecay = 0.7 ;

    /* terminate after 2*n + nslow iterations without strict improvement in
       either function value or gradient */
    Parm->nslow = 1000 ;

    /* T => attempt quadratic interpolation in line search when
                |f_k+1 - f_k|/f_k <= QuadCutOff
       F => no quadratic interpolation step */
    Parm->QuadStep = TRUE ;
    Parm->QuadCutOff = 1.e-12 ;

    /* maximum factor by which a quad step can reduce the step size */
    Parm->QuadSafe = 1.e-10 ;

    /* for a QuadStep, function evaluated on interval
       [psi_lo, phi_hi]*psi2*previous step */
    Parm->psi_lo = 0.1 ;
    Parm->psi_hi = 10. ;

    /* when the function is approximately quadratic, use gradient at
       psi1*psi2*previous step for estimating initial stepsize */
    Parm->psi1 = 1.0 ;

    /* parameter used in cost error estimate for quadratic restart criterion */
    Parm->qeps = 1.e-12 ;

    /* treat cost as quadratic if
       |1 - (cost change)/(quadratic cost change)| <= qrule */
    Parm->qrule = 1.e-8 ;

    /* number of iterations the function is nearly quadratic before
       it is treated as quadratic */
    Parm->qrestart = 6 ;

    /* T => when possible, use a cubic step in the line search */
    Parm->UseCubic = TRUE ;

    /* use cubic step when |f_k+1 - f_k|/|f_k| > CubicCutOff */
    Parm->CubicCutOff = 1.e-12 ;

    /* |f| < SmallCost*starting cost => skip QuadStep and set PertRule = FALSE*/
    Parm->SmallCost = 1.e-30 ;

    /* maximum factor secant step increases stepsize in expansion phase */
    Parm->ExpandSafe = 200. ;

    /* factor by which secant step is amplified during expansion phase
       where minimizer is bracketed */
    Parm->SecantAmp = 1.05 ;

/* =================== LINE SEARCH ========================================== */
    /* T => use approximate Wolfe line search
       F => use ordinary Wolfe, switch to approximate when
              |fR - f| < ApproxSwitchFactor*|fR|
       ApproxSwitchFactor = 0 => ordinary Wolfe condition if approxstep = F */
    Parm->approxstep = FALSE ;
    Parm->ApproxSwitchFactor = 1.e-3 ;

    /* As the cg_descent converges, the function values typically
       approach a constant value. When this happens, the cubic interpolation
       step in the line search loses its accuracy and it is better to
       use a secant step based on the derivative of the objective function.
       The cost has converged when the relative change in the objective
       function <= CostConverge */
    Parm->CostConverge = 1.e-10 ;

    /* Wolfe line search parameter between 0 and 0.5
       phi (a) - phi (0) <= cgdelta phi'(0) */
    Parm->cgdelta = 0.1 ;

    /* In a Wolfe line search, it is also required that gnew'*dk >=
       Wolfe_sigma*gk'*dk where sigma is between cgdelta and 1. */
    Parm->cgsigma = 0.9 ;

    /* maximum number of attempts to find an acceptable stepsize */
    Parm->maxsteps = (int) 99 ;

    /* decay factor for bracket interval, if the interval width does
       not decay by this factor when high-order methods are used,
       then employ a bisection step */
    Parm->stepdecay = 0.66 ;

    /* When an infinite or nan objective value is encountered, the
       stepsize is reduced in an effort to find a finite objective
       value. infdecay is the initial decay factor that is used when an
       infinite or nan objective value is encountered, ninf_tries is
       the number of attempts we make to find a finite objective value,
       and infdecay_rate is a factor by which indecay is multiplied
       after each attempt to find a finite objective value. */

    /* initial contraction factor when hitting an infinite objective value */
    Parm->cg_infdecay = 0.5 ;

    /* infdecay is multiplied by infdecay_rate after each attempt to
       find a finite objective value */
    Parm->cg_infdecay_rate = 0.9 ;

    /* number of times we try to find a finite objective value before
       declaring an error */
    Parm->cg_ninf_tries = 20 ;

    /* growth factor in search for initial bracket interval */
    Parm->rho = 5.0 ;

    /* factor by which rho grows during expansion phase where minimizer is
       bracketed */
    Parm->RhoGrow = 2.0 ;

    /* if the currect derivative <= BigDfactor * derivative at starting point,
       then we evaluate the function at the current point even though
       AvoidFeval is TRUE. When AvoidFeval is TRUE when we think that
       the function values have converged; in this case, line search
       is mostly based on the derivative. */
    Parm->BigDfactor = 1000. ;

    /* When performing an approximate Wolfe line search, we require
       that the new function value <= perturbation of the prior
       function value where the perturbation tries to take into
       rounding errors associated with the function value.
       There are two different perturbations:
       PertRule = 1 => fpert is f + Parm->pert_eps*|f| (relative change)
       PertRule = 0 => fpert is f + Parm->pert_eps     (absolute change) */
    Parm->PertRule = 1 ;
    Parm->pert_eps = 1.e-6 ;

    /* Maximum number of contractions in cg_contract. If it cannot find a
       step that either satisfies the Wolfe conditions or which has derivative
       >= 0 within ncontract attempts, then it is felt that fpert is too
       small and it will be increased so that the current function value
       is less than fpert. To increase fpert, we increase the value of
       pert_eps which is used to compute fpert. */
    Parm->ncontract = (int) 5 ;

    /* When pert_eps is increased, we multiply the new value by the growth
       factor eps_grow to ensure a healthy growth in pert_eps. */
    Parm->eps_grow = 20.0 ;

    /* Maximum number of times that pert_eps is recomputed before a line
       search error is declared. */
    Parm->neps = (int) 5 ;

    /* starting guess for line search =
         psi0 ||x_0||_infty over ||g_0||_infty if x_0 != 0
         psi0 |f(x_0)|/||g_0||_2               otherwise */
    Parm->psi0 = .01 ;      /* factor used in starting guess for iteration 1 */

    /* when starting a new cg iteration, our initial guess for the line
       search stepsize is psi2*previous step */
    Parm->psi2 = 2.0 ;

    /* lower bound for cg update parameter beta is BetaLower*d_k'g_k/ ||d_k||^2
       (lower bound for beta needed to ensure global convergence of cg) */
    Parm->BetaLower = 0.4 ;

    /* value of the parameter theta in the cg_descent update formula:
       W. W. Hager and H. Zhang, A survey of nonlinear conjugate gradient
       methods, Pacific Journal of Optimization, 2 (2006), pp. 35-58. */
    Parm->theta = 1.0 ;

    /* choose theta adaptively if AdaptiveTheta = T */
    Parm->AdaptiveTheta = FALSE ;

/* ============ LIMITED MEMORY CG PARAMETERS ================================ */
    /* SubCheck and SubSkip control the frequency with which the subspace
       condition is checked. It it checked for SubCheck*mem iterations and
       if it is not activated, then it is skipped for Subskip*mem iterations
       and Subskip is doubled. Whenever the subspace condition is satisfied,
       SubSkip is returned to its original value. */
    Parm->SubCheck = 8 ;
    Parm->SubSkip = 4 ;

    /* when relative distance from current gradient to subspace <= eta0,
       enter subspace if subspace dimension = mem (eta0 = 0 means gradient
       inside subspace) */
    Parm ->eta0 = 0.001 ; /* corresponds to eta0*eta0 in the paper */

    /* when relative distance from current gradient to subspace >= eta1,
       leave subspace (eta1 = 1 means gradient orthogonal to subspace) */
    Parm->eta1 = 0.900 ; /* corresponds to eta1*eta1 in the paper */

    /* when relative distance from current gradient to subspace <= eta2,
       always enter subspace (invariant space) */
    Parm->eta2 = 1.e-10 ;

    /* number of expansion steps in cg when solving bound constrained problem
       through pasa */
    Parm->nexpand = 10 ;

    /* factor by which the stepsize grows in each expansion step */
    Parm->step_growth = 4. ;
}
