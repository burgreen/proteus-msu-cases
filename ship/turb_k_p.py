# turb_k_p.py

from proteus.default_p  import *
from proteus.mprans     import Kappa
from user_param         import *
from main_param         import *

LevelModelType = Kappa.LevelModel

dissipation_model_flag = 1
if user_param.useRANS == 2: dissipation_model_flag = 2
if user_param.useRANS == 3: dissipation_model_flag = 3

RD_model          = None
LS_model          = None
ME_model          = 1
dissipation_model = 2

coefficients = Kappa.Coefficients(
  V_model                = 0 + int(movingDomain),
  ME_model               = ME_model,
  LS_model               = LS_model,
  RD_model               = RD_model,
  dissipation_model      = dissipation_model,
  dissipation_model_flag = dissipation_model_flag, # 1=K-epsilon, 2=K-omega 1998, 3=K-omega 1988
  useMetrics             = useMetrics,             # main_param.useMetrics,
  rho_0                  = rho_0, # main_param
  nu_0                   = nu_0, # main_param
  rho_1                  = rho_1, # main_param
  nu_1                   = nu_1, # main_param
  #g                     = gravity, # main_param
  g                      = numpy.array([user_param.gravity[0],
                                        user_param.gravity[1],
                                        user_param.gravity[2]], dtype='d'),
  nd                     = nd, #main_param.nd,
  c_mu                   = 0.09,
  sigma_k                = 1.0,
  sc_uref                = kappa_sc_uref, # main_param
  sc_beta                = kappa_sc_beta # main_param
)

'''
dirichletConditions = {
  0: lambda x, flag: domain.bc[flag].k_dirichlet.init_cython()
}

advectiveFluxBoundaryConditions = {
  0: lambda x, flag: domain.bc[flag].k_advective.init_cython()
}

diffusiveFluxBoundaryConditions = {
  0: {},
  1: {1: lambda x, flag: domain.bc[flag].k_diffusive.init_cython()}
}
'''

dirichletConditions = {
  0: lambda x,flag: domain.bc[domain.meshtag_bcIdx[flag]].k_dirichlet.init_cython()
}

advectiveFluxBoundaryConditions = {
  0: lambda x,flag: domain.bc[domain.meshtag_bcIdx[flag]].k_advective.init_cython()
}

diffusiveFluxBoundaryConditions = {
  0: { 0: lambda x,flag: domain.bc[domain.meshtag_bcIdx[flag]].k_diffusive.init_cython() }
}

turb_k = user_param.IC_field_value['turb_k']

initialConditions = {
  0: user_param.IC_field_constant(turb_k)
}

'''
class ConstantIC:
    def __init__(self, cval=0.):
        self.cval = cval
    def uOfXT(self, x, t):
        return self.cval

initialConditions = {0: ConstantIC(cval=0)}
'''
