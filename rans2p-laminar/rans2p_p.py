# rans2p_p.py

from proteus            import *
from proteus.default_p  import *
from proteus.mprans     import RANS2P
from main_param         import *
from user_param         import *

LevelModelType = RANS2P.LevelModel

LS_model                = 2
if modeVF < 2: LS_model = None

# 2-phase definition
coefficients_2p = RANS2P.Coefficients(
  epsFact                = epsFact_viscosity,
  sigma                  = 0.0,
  rho_0                  = rho_0,
  nu_0                   = nu_0,
  rho_1                  = rho_1,
  nu_1                   = nu_1,
  g                      = g,
  nd                     = nd,
  VF_model               = 1,
  LS_model               = LS_model,
  epsFact_density        = epsFact_density,
  stokes                 = False,
  useVF                  = useVF,
	useRBLES               = useRBLES,
	useMetrics             = useMetrics,
  eb_adjoint_sigma       = 1.0,
  forceStrongDirichlet   = 0,
  turbulenceClosureModel = ns_closure
)

# 1-phase definition
coefficients_1p = RANS2P.Coefficients(
  epsFact              = epsFact_viscosity,
  rho_0                = rho_0,
  nu_0                 = nu_0,
  rho_1                = rho_1,
  nu_1                 = nu_1,
  g                    = g,
  nd                   = nd,
  LS_model             = None,
  epsFact_density      = epsFact_density,
  stokes               = False,#useStokes,
  #gwb forceStrongDirichlet = True,
  forceStrongDirichlet = False,
  eb_adjoint_sigma     = 1.0,
  eb_penalty_constant  = 10.0,
  useRBLES             = 0.0,
  useMetrics           = 1.0
)

coefficients                           = coefficients_1p 
if user_param.nphase > 1: coefficients = coefficients_2p 

dirichletConditions = {
  0: lambda x,flag: domain.bc[domain.meshtag_bcIdx[flag]].p_dirichlet.init_cython(),
  1: lambda x,flag: domain.bc[domain.meshtag_bcIdx[flag]].u_dirichlet.init_cython(),
  2: lambda x,flag: domain.bc[domain.meshtag_bcIdx[flag]].v_dirichlet.init_cython(),
  3: lambda x,flag: domain.bc[domain.meshtag_bcIdx[flag]].w_dirichlet.init_cython()
}
advectiveFluxBoundaryConditions = {
  0: lambda x,flag: domain.bc[domain.meshtag_bcIdx[flag]].p_advective.init_cython(),
  1: lambda x,flag: domain.bc[domain.meshtag_bcIdx[flag]].u_advective.init_cython(),
  2: lambda x,flag: domain.bc[domain.meshtag_bcIdx[flag]].v_advective.init_cython(),
  3: lambda x,flag: domain.bc[domain.meshtag_bcIdx[flag]].w_advective.init_cython() 
}
diffusiveFluxBoundaryConditions = {
  0: {},
  1: {1: lambda x,flag: domain.bc[domain.meshtag_bcIdx[flag]].u_diffusive.init_cython()},
  2: {2: lambda x,flag: domain.bc[domain.meshtag_bcIdx[flag]].v_diffusive.init_cython()},
  3: {3: lambda x,flag: domain.bc[domain.meshtag_bcIdx[flag]].w_diffusive.init_cython()} 
}
fluxBoundaryConditions = {
  0: 'mixedFlow',
  1: 'mixedFlow',
  2: 'mixedFlow',
  3: 'mixedFlow'
}

p0 = user_param.IC_field_value['p']
u0 = user_param.IC_field_value['u']
v0 = user_param.IC_field_value['v']
w0 = user_param.IC_field_value['w']

IC_p                           = user_param.IC_field_constant(p0)
if user_param.nphase > 1: IC_p = user_param.IC_field_p()

initialConditions = {
  0: IC_p,
  1: user_param.IC_field_constant(u0),
  2: user_param.IC_field_constant(v0),
  3: user_param.IC_field_constant(w0)
}
