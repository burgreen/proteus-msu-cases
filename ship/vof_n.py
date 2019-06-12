from proteus    import *
from vof_p      import *
from main_param import *
from user_param import *

timeIntegration = BackwardEuler_cfl
stepController  = Min_dt_controller

femSpaces = {0:basis}

massLumping       = False
numericalFluxType = VOF.NumericalFlux
conservativeFlux  = None
subgridError      = VOF.SubgridError(coefficients=coefficients,nd=nd)
shockCapturing    = VOF.ShockCapturing(coefficients,nd,shockCapturingFactor=vof_shockCapturingFactor,lag=vof_lag_shockCapturing)

fullNewtonFlag = True
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

linear_solver_options_prefix        = 'vof_'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest         = 'r-true'

tolFac             = 0.0
linTolFac          = 0.01
l_atol_res         = 0.01*vof_nl_atol_res
nl_atol_res        = vof_nl_atol_res

if 'abs_tol_vof' in user_param.tols: 
  nl_atol_res  = user_param.tols['abs_tol_vof']

useEisenstatWalker = False

maxNonlinearIts = 50
maxLineSearches = 0

#----- gwb mods
#linTolFac          = 1.e-5*0
#l_atol_res         = 2.e-6
#nl_atol_res        = 3.e-4
#print( 'vof_nl_atol_res', vof_nl_atol_res )

