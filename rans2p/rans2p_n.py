# rans2p_n.py

from proteus    import *
from rans2p_p   import *
from main_param import *
from user_param import *

timeIntegration = BackwardEuler_cfl
stepController  = Min_dt_controller

femSpaces = {0:basis, 1:basis, 2:basis, 3:basis}

massLumping       = False
numericalFluxType = None
conservativeFlux  = None

numericalFluxType = RANS2P.NumericalFlux
subgridError      = RANS2P.SubgridError(coefficients,nd,lag=ns_lag_subgridError,hFactor=hFactor)
shockCapturing    = RANS2P.ShockCapturing(coefficients,nd,ns_shockCapturingFactor,lag=ns_lag_shockCapturing)

fullNewtonFlag            = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton

nonlinearSmoother = None
linearSmoother    = SimpleNavierStokes3D

matrix = SparseMatrix

multilevelLinearSolver                 = KSP_petsc4py
levelLinearSolver                      = KSP_petsc4py
if useOldPETSc: multilevelLinearSolver = PETSc
if useOldPETSc: levelLinearSolver      = PETSc
if useSuperlu:  multilevelLinearSolver = LU
if useSuperlu:  levelLinearSolver      = LU

linear_solver_options_prefix        = 'rans2p_'
levelNonlinearSolverConvergenceTest = user_param.knob['ns']['nl_test'] 
linearSolverConvergenceTest         = 'r-true'

tolFac             = 0.0
linTolFac          = 0.0
l_atol_res         = 0.01*ns_nl_atol_res
nl_atol_res        = ns_nl_atol_res
nl_atol_res        = max( ns_nl_atol_res, user_param.knob['ns']['nl_atol'] )

useEisenstatWalker = False
maxNonlinearIts    = user_param.knob['ns']['nl_its']
maxLineSearches    = 0
conservativeFlux   = {0:'pwl-bdm-opt'}
