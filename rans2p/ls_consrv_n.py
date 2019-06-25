from proteus      import *
from ls_consrv_p  import *
from main_param   import *
from user_param   import *

timeIntegrator  = ForwardIntegrator
timeIntegration = NoIntegration

femSpaces = {0:basis}

subgridError      = None
massLumping       = False
numericalFluxType = DoNothing
conservativeFlux  = None
shockCapturing    = None

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

linear_solver_options_prefix = 'mcorr_'
levelNonlinearSolverConvergenceTest = user_param.knob['mcorr']['nl_test']
linearSolverConvergenceTest  = 'r-true'

tolFac             = 0.0
linTolFac          = 0.0
l_atol_res         = 0.01*mcorr_nl_atol_res
nl_atol_res        = min( mcorr_nl_atol_res, user_param.knob['mcorr']['nl_atol'] )

useEisenstatWalker = False

maxNonlinearIts = user_param.knob['mcorr']['nl_its']
maxLineSearches = 0
