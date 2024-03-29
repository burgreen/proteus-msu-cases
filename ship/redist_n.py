from proteus    import *
from redist_p   import *
from main_param import *

tolFac             = 0.0
linTolFac          = 0.0
l_atol_res         = 0.01*rd_nl_atol_res
nl_atol_res        = rd_nl_atol_res

if 'nl_atol_redist' in user_param.tols: 
  nl_atol_res  = user_param.tols['nl_atol_redist']

useEisenstatWalker = False

if redist_Newton:
    timeIntegration                     = NoIntegration
    stepController                      = Newton_controller
    maxNonlinearIts                     = 25
    maxLineSearches                     = 0
    nonlinearSolverConvergenceTest      = 'rits'
    levelNonlinearSolverConvergenceTest = 'rits'
    linearSolverConvergenceTest         = 'r-true'
else:
    timeIntegration                     = BackwardEuler_cfl
    stepController                      = RDLS.PsiTC
    runCFL                              = 0.5
    psitc['nStepsForce']                = 6
    psitc['nStepsMax']                  = 25
    psitc['reduceRatio']                = 3.0
    psitc['startRatio']                 = 1.0
    rtol_res[0]                         = 0.0
    atol_res[0]                         = nl_atol_res
    useEisenstatWalker                  = False
    maxNonlinearIts                     = 1
    maxLineSearches                     = 0
    nonlinearSolverConvergenceTest      = 'rits'
    levelNonlinearSolverConvergenceTest = 'rits'
    linearSolverConvergenceTest         = 'r-true'

femSpaces = {0:basis}
       
massLumping       = False
numericalFluxType = DoNothing    
conservativeFlux  = None
subgridError      = RDLS.SubgridError(coefficients,nd)
shockCapturing    = RDLS.ShockCapturing(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=rd_lag_shockCapturing)

fullNewtonFlag             = True
multilevelNonlinearSolver  = Newton
levelNonlinearSolver       = Newton

nonlinearSmoother = NLGaussSeidel
linearSmoother    = None

matrix = SparseMatrix

multilevelLinearSolver                 = KSP_petsc4py
levelLinearSolver                      = KSP_petsc4py
if useOldPETSc: multilevelLinearSolver = PETSc
if useOldPETSc: levelLinearSolver      = PETSc
if useSuperlu:  multilevelLinearSolver = LU
if useSuperlu:  levelLinearSolver      = LU

linear_solver_options_prefix = 'rdls_'
