from proteus                        import *
from proteus.default_p              import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.mprans                 import VOF
from main_param                     import *
from user_param                     import *

LevelModelType = VOF.LevelModel

RD_model                = 3
LS_model                = 2
if modeVF < 2: RD_model = None
if modeVF < 2: LS_model = None

coefficients = VOF.Coefficients(
    LS_model   = LS_model,
    V_model    = 0,
    RD_model   = RD_model,
    ME_model   = 1,
    checkMass  = False,
    useMetrics = useMetrics,
    epsFact    = epsFact_vof,
    sc_uref    = vof_sc_uref,
    sc_beta    = vof_sc_beta 
)
 
def getDBC_vof(x,flag):
    return None

def getAFBC_vof(x,flag):
    return lambda x,t: 0.0
    #return domain.bc[domain.meshtag_bcIdx[flag]].vof_advective.init_cython()

dirichletConditions = {
  #0: getDBC_vof
  0: lambda x,flag: domain.bc[domain.meshtag_bcIdx[flag]].vof_dirichlet.init_cython()
}

advectiveFluxBoundaryConditions = {
  #0: getAFBC_vof
  0: lambda x,flag: domain.bc[domain.meshtag_bcIdx[flag]].vof_advective.init_cython()
}

diffusiveFluxBoundaryConditions = {
  0: {}
}

class sdf_smoothed:
    def uOfXT(self,x,t):
        sd = user_param.IC_signed_distance(x)
        fac = epsFact_consrv_heaviside * user_param.mesh_nominal_spacing
        return smoothedHeaviside( fac, sd )
	    
initialConditions  = {
  0: sdf_smoothed() 
}
