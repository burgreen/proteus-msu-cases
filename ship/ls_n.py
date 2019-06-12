from proteus     import *
from ls_p        import *
from user_param  import *

timeIntegration = BackwardEuler_cfl
stepController  = Min_dt_controller

femSpaces = {0:basis}

massLumping       = False
conservativeFlux  = None
numericalFluxType = NCLS.NumericalFlux
subgridError      = NCLS.SubgridError(coefficients,nd)
shockCapturing    = NCLS.ShockCapturing(coefficients,nd,shockCapturingFactor=ls_shockCapturingFactor,lag=ls_lag_shockCapturing)

fullNewtonFlag            = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton

nonlinearSmoother = None
linearSmoother    = None

matrix = SparseMatrix

multilevelLinearSolver                 = KSP_petsc4py
levelLinearSolver                      = KSP_petsc4py
if useOldPETSc: multilevelLinearSolver = PETSc
if useOldPETSc: levelLinearSolver      = PETSc
if useSuperlu:  multilevelLinearSolver = LU
if useSuperlu:  levelLinearSolver      = LU

linear_solver_options_prefix        = 'ncls_'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest         = 'r-true'

tolFac             = 0.0
linTolFac          = 0.01
l_atol_res         = 0.01*ls_nl_atol_res
nl_atol_res        = 1.00*ls_nl_atol_res
l_atol_res         = 1.e-0;
nl_atol_res        = 1.e-0;

if 'nl_atol_ls' in user_param.tols: 
  nl_atol_res  = user_param.tols['nl_atol_ls']

useEisenstatWalker = False

maxNonlinearIts = 80
maxLineSearches = 0
