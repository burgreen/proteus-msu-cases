from proteus            import *
from proteus.default_p  import *
from proteus.mprans     import Dissipation
from main_param         import *
from user_param  import *

LevelModelType = Dissipation.LevelModel
RD_model       = None
LS_model       = None
ME_model       = 2

#kappa_model   = 1
#should not include this dissipation_model = 3

dissipation_model_flag = 1
if user_param.useRANS == 2: dissipation_model_flag = 2
if user_param.useRANS == 3: dissipation_model_flag = 3

print('dissipation_model_flag =', dissipation_model_flag)

gravity  = numpy.array([user_param.gravity[0],
                        user_param.gravity[1],
                        user_param.gravity[2]], dtype='d'),

coefficients = Dissipation.Coefficients(
    V_model     = 0 + int(movingDomain),
    ME_model    = ME_model,
    LS_model    = LS_model,
    RD_model    = RD_model,
    #kappa_model = kappa_model,
    dissipation_model_flag = dissipation_model_flag, #1 -- K-epsilon, 2 -- K-omega 1998, 3 -- K-omega 1988
    useMetrics  = useMetrics, #user_param
    rho_0       = rho_0,
    rho_1       = rho_1,
    nu_0        = nu_0,
    nu_1        = nu_1,
    g           = gravity
    c_mu        = 0.09,
    sigma_e     = 1.0,
    sc_uref     = dissipation_sc_uref,
    sc_beta     = dissipation_sc_beta
)

dirichletConditions = {
  0: lambda x, flag: domain.bc[flag].dissipation_dirichlet.init_cython()
}

advectiveFluxBoundaryConditions = {
  0: lambda x, flag: domain.bc[flag].dissipation_advective.init_cython()
}

diffusiveFluxBoundaryConditions = {
  0: {},
  1: { 1: lambda x, flag: domain.bc[flag].dissipation_diffusive.init_cython() }
}

turb_e = user_param.IC_field_value['turb_e']

initialConditions = {
  0: user_param.IC_field_constant(turb_e)
}

'''
class ConstantIC:
    def __init__(self, cval=0.):
        self.cval = cval
    def uOfXT(self, x, t):
        return self.cval

kInflow = 0.
#xw hardwired
dissipationInflow =0.0
#dissipationInflow = coefficients.c_mu*kInflow**(1.5)/(0.03*tank_dim[nd-1])
if useRANS >= 2:
    dissipationInflow = dissipationInflow/(kInflow+1.0e-12)

initialConditions = {0: ConstantIC(cval=dissipationInflow*0.001)}
'''
