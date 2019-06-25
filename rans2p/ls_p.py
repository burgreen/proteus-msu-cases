from proteus           import *
from proteus.default_p import *
from proteus.mprans    import NCLS
from main_param        import *
from user_param        import *

LevelModelType = NCLS.LevelModel

coefficients = NCLS.Coefficients(
    V_model    = 0,
    RD_model   = 3,
    ME_model   = 2,
    checkMass  = False, 
    useMetrics = useMetrics,
    epsFact    = epsFact_consrv_heaviside,
    sc_uref    = ls_sc_uref,
    sc_beta    = ls_sc_beta
)


dirichletConditions = {
  0: lambda x,flag: domain.bc[domain.meshtag_bcIdx[flag]].clsvof_dirichlet.init_cython()
}

advectiveFluxBoundaryConditions = {
}

diffusiveFluxBoundaryConditions = {
  0: {}
}

initialConditions = {
  0: user_param.IC_field_phi()
}


