# turb_e_n.py

import turb_e_p as physics

from proteus.default_n import *
from proteus           import StepControl
from proteus           import TimeIntegration
from proteus           import NonlinearSolvers
from proteus           import LinearSolvers
from proteus           import LinearAlgebraTools
from proteus.mprans    import Dissipation
from proteus           import Context
from main_param        import *
from user_param        import *

timeIntegration = TimeIntegration.BackwardEuler_cfl
stepController  = StepControl.Min_dt_controller

if timeIntegration == "VBDF": timeIntegration = TimeIntegration.VBDF
if timeIntegration == "VBDF": timeOrder = 2

femSpaces                 = {0: basis}
elementQuadrature         = elementQuadrature
elementBoundaryQuadrature = elementBoundaryQuadrature

massLumping       = False
numericalFluxType = Dissipation.NumericalFlux
conservativeFlux  = None
subgridError      = Dissipation.SubgridError(coefficients=physics.coefficients,
                                             nd=nd)
shockCapturing    = Dissipation.ShockCapturing(coefficients=physics.coefficients,
                                               nd=nd,
                                               shockCapturingFactor=dissipation_shockCapturingFactor,
                                               lag=dissipation_lag_shockCapturing)
fullNewtonFlag            = True
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver      = NonlinearSolvers.Newton

nonlinearSmoother = None
linearSmoother    = None

matrix = LinearAlgebraTools.SparseMatrix

multilevelLinearSolver = LinearSolvers.LU
levelLinearSolver      = LinearSolvers.LU

if not useSuperlu: multilevelLinearSolver = LinearSolvers.KSP_petsc4py
if not useSuperlu: levelLinearSolver      = LinearSolvers.KSP_petsc4py

linear_solver_options_prefix        = 'dissipation_'
levelNonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest         = 'rits'

tolFac             = 0.0
linTolFac          = 0.0
l_atol_res         = 0.001*dissipation_nl_atol_res
nl_atol_res        = dissipation_nl_atol_res

if 'abs_tol_turb_e' in user_param.tols: 
  nl_atol_res  = user_param.tols['abs_tol_turb_e']

useEisenstatWalker = False

maxNonlinearIts = 50
maxLineSearches = 0
