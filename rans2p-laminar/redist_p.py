# The redistancing equation in the sloshbox test problem.

from math              import *
from proteus           import *
from proteus.default_p import *
from proteus.mprans    import RDLS
from main_param        import *
from user_param        import *

LevelModelType = RDLS.LevelModel

coefficients = RDLS.Coefficients(
    applyRedistancing = True,
    epsFact           = epsFact_redistance,
    nModelId          = 2,
    rdModelId         = 3,
    useMetrics        = useMetrics 
)

def getDBC_rd(x,flag):
    pass
    
dirichletConditions = {
  0: getDBC_rd
  #0: lambda x,flag: None
}

weakDirichletConditions = {
  0: RDLS.setZeroLSweakDirichletBCsSimple
}

advectiveFluxBoundaryConditions =  {
}

diffusiveFluxBoundaryConditions = {
  0: {}
}

initialConditions = {
  0: user_param.IC_field_phi()
}
