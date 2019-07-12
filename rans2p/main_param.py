# main_param.py

from math import *

import proteus.MeshTools
from   proteus            import Domain
from   proteus.msu        import MeshFileDomain
from   proteus.default_n  import *   
from   proteus.Profiling  import logEvent
from   case               import std_knob as std_knob
from   user_param         import *

# Discretization -- input options    

useOldPETSc     = False
useSuperlu      = False
spaceOrder      = user_param.spaceOrder
useHex          = False
useRBLES        = 0.0
useMetrics      = 1.0
applyCorrection = True
useVF           = 1.0
modeVF          = 0
#redist_Newton   = True # assigned via useMetrics 
useRANS         = user_param.useRANS      # 0 -- None # 1 -- K-Epsilon # 2 -- K-Omega
he              = 0.1    # mesh size
weak_bc_penalty_constant = 100.0
ns_forceStrongDirichlet = False

runCFL = std_knob.knob['misc']['runCFL']

if user_param.nphase > 1: modeVF = 2

# Input checks

if spaceOrder not in [1,2]:
    print( "INVALID: spaceOrder" + spaceOrder )
    sys.exit()    
    
if useRBLES not in [0.0, 1.0]:
    print( "INVALID: useRBLES" + useRBLES  )
    sys.exit()

if useMetrics not in [0.0, 1.0]:
    print( "INVALID: useMetrics" )
    sys.exit()
    
#  Discretization   

nd = 3

if spaceOrder == 1:
   hFactor = 1.0
   if useHex:
      basis = C0_AffineLinearOnCubeWithNodalBasis
      elementQuadrature         = CubeGaussQuadrature(nd,2)
      elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,2)     	 
   else:
      basis = C0_AffineLinearOnSimplexWithNodalBasis
      elementQuadrature         = SimplexGaussQuadrature(nd,3)
      elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3) 	    
elif spaceOrder == 2:
    hFactor = 0.5
    if useHex:    
      basis = C0_AffineLagrangeOnCubeWithNodalBasis
      elementQuadrature         = CubeGaussQuadrature(nd,4)
      elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,4)    
    else:    
      basis = C0_AffineQuadraticOnSimplexWithNodalBasis	
      elementQuadrature         = SimplexGaussQuadrature(nd,4)
      elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

# Domain and mesh

nLevels = 1
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
use_petsc4py=True

print( '----------------', user_param.filename )

domain = MeshFileDomain.MeshFileDomain(user_param.filename,3) 

print( '----------------', domain )

for key in user_param.bc_zone.keys():
  bc = user_param.bc_zone[key]
  domain.bc_zone_define( name=key, meshtag=bc['meshtag'],  condition=bc['condition'] )

for key in user_param.ele_zone.keys():
  bc = user_param.ele_zone[key]
  domain.ele_zone_define( name=key, meshtag=bc['meshtag'],  condition=bc['condition'] )

# Time stepping

dt_init  = user_param.dt_init
dt_fixed = user_param.dt_fixed
dt_fixed_steps = user_param.dt_fixed_steps

# Numerical parameters

print('---------------------------------------')
print(' case-specific knob values')
for e in user_param.knob.keys():
  for v in user_param.knob[e].keys():
    if user_param.knob[e][v] != std_knob.knob[e][v]: 
      print('  knob['+e+']['+v+'] = ', user_param.knob[e][v] )
print('---------------------------------------')

ns_forceStrongDirichlet = False

if useMetrics:
    ns_shockCapturingFactor = user_param.knob['ns']['sc_fac']
    ns_lag_shockCapturing   = user_param.knob['ns']['sc_lag']
    ns_lag_subgridError     = user_param.knob['ns']['sge_lag']

    ls_shockCapturingFactor = user_param.knob['ls']['sc_fac']
    ls_lag_shockCapturing   = user_param.knob['ls']['sc_lag']
    ls_sc_uref              = user_param.knob['ls']['sc_uref']
    ls_sc_beta              = user_param.knob['ls']['sc_beta']

    vof_shockCapturingFactor = user_param.knob['vof']['sc_fac']
    vof_lag_shockCapturing   = user_param.knob['vof']['sc_lag']
    vof_sc_uref              = user_param.knob['vof']['sc_uref']
    vof_sc_beta              = user_param.knob['vof']['sc_beta']

    rd_shockCapturingFactor  = user_param.knob['rd']['sc_fac']
    rd_lag_shockCapturing    = user_param.knob['rd']['sc_lag']
    redist_Newton            = user_param.knob['rd']['newton']

    kappa_shockCapturingFactor = user_param.knob['turb_k']['sc_fac']
    kappa_lag_shockCapturing   = user_param.knob['turb_k']['sc_lag']
    kappa_sc_uref              = user_param.knob['turb_k']['sc_uref']
    kappa_sc_beta              = user_param.knob['turb_k']['sc_beta']

    dissipation_shockCapturingFactor = user_param.knob['turb_e']['sc_fac']
    dissipation_lag_shockCapturing   = user_param.knob['turb_e']['sc_lag']
    dissipation_sc_uref              = user_param.knob['turb_e']['sc_uref']
    dissipation_sc_beta              = user_param.knob['turb_e']['sc_beta']

    epsFact_density    = user_param.knob['epsFact']['std']
    epsFact_viscosity  = epsFact_curvature = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = user_param.knob['epsFact']['rd']
    epsFact_consrv_diffusion = user_param.knob['epsFact']['consrv_diffusion']

else:
    ns_shockCapturingFactor  = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref  = 1.0
    vof_sc_beta  = 1.0
    rd_shockCapturingFactor  = 0.9
    rd_lag_shockCapturing = False
    epsFact_density    = 1.5
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0

# nl_atol_res assignments

tol_std = max( 1.0e-8, 0.1*he**2/2.0 )
tol_rd  = max( 1.0e-8, 0.1*he )
tol_std = 1.0e-5
tol_rd  = 1.0e-5

ns_nl_atol_res          = tol_std
vof_nl_atol_res         = tol_std
ls_nl_atol_res          = tol_std
rd_nl_atol_res          = tol_rd
mcorr_nl_atol_res       = tol_std
kappa_nl_atol_res       = tol_std
dissipation_nl_atol_res = tol_std

#turbulence: 1-classic-smagorinsky, 2-dynamic-smagorinsky, 3-k-epsilon, 4-k-omega

ns_closure = 0 
if useRANS == 1: ns_closure = 3
if useRANS == 2: ns_closure = 4

# fluid phases

phase_0 = user_param.phase[0]
phase_1                               = user_param.phase[0]
if len(user_param.phase) > 1: phase_1 = user_param.phase[1] 

rho_0 = user_param.fluid[phase_0]['rho']
mu_0  = user_param.fluid[phase_0]['mu']
nu_0  = mu_0/rho_0

rho_1 = user_param.fluid[phase_1]['rho']
mu_1  = user_param.fluid[phase_1]['mu']
nu_1  = mu_1/rho_1

print( 'nphase, phases =', user_param.nphase, phase_0, phase_1 )

# Surface tension

sigma_01 = 0.0

# Gravity

g = user_param.gravity
