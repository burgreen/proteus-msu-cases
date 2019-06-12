# turb_k_n.py

import turb_k_p as physics

from proteus.default_n  import *
from proteus            import StepControl
from proteus            import TimeIntegration
from proteus            import NonlinearSolvers
from proteus            import LinearSolvers
from proteus            import LinearAlgebraTools
from proteus.mprans     import Kappa
from main_param         import *
from user_param         import * 

timeIntegration = TimeIntegration.BackwardEuler_cfl
stepController  = StepControl.Min_dt_controller

if timeIntegration == "VBDF": timeIntegration = TimeIntegration.VBDF
if timeIntegration == "VBDF": timeOrder = 2

femSpaces = {0: basis}

massLumping       = False
numericalFluxType = Kappa.NumericalFlux
conservativeFlux  = None
subgridError      = Kappa.SubgridError(coefficients=physics.coefficients,
                                       nd=nd)
shockCapturing    = Kappa.ShockCapturing(coefficients=physics.coefficients,
                                         nd=nd,
                                         #shockCapturingFactor=user_param.kappa_shockCapturingFactor,
                                         #lag=user_param.kappa_lag_shockCapturing)
                                         shockCapturingFactor=kappa_shockCapturingFactor,
                                         lag=kappa_lag_shockCapturing)

fullNewtonFlag            = True
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver      = NonlinearSolvers.Newton

nonlinearSmoother = None
linearSmoother    = None
#printNonlinearSolverInfo = True

matrix = LinearAlgebraTools.SparseMatrix

if not useOldPETSc and not useSuperlu:
    multilevelLinearSolver = LinearSolvers.KSP_petsc4py
    levelLinearSolver      = LinearSolvers.KSP_petsc4py
else:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver      = LinearSolvers.LU

linear_solver_options_prefix        = 'kappa_'
levelNonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest         = 'rits'

tolFac             = 0.
linTolFac          = 0.
l_atol_res         = 0.001*kappa_nl_atol_res
nl_atol_res        = kappa_nl_atol_res

if 'nl_atol_turb_k' in user_param.tols: 
  nl_atol_res  = user_param.tols['nl_atol_turb_k']


useEisenstatWalker = False

maxNonlinearIts = 50
maxLineSearches = 0
