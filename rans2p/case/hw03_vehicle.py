# user_param.py  case: hw03 vehicle situated midstream

import math 
from proteus import default_p
from . import std_knob as std_knob
from . import std_bc as std_bc

# time stepping

dt_init  = 0.0001
dt_fixed = 0.0001
dt_fixed_steps = 2 

# case/geometry properties

mesh = {}
mesh['filename'] = 'mesh/hw03.tetgen'
mesh['nominal_spacing'] = 0.5

case = {}
case['num_phases']  = 2
case['gravity_dir'] = '-z' # choose from: '+x','-x','+y','-y','+z','-z'
case['waterline']   = 2.0

IC_field_value = {} # initial conditions
IC_field_value['u'] = 1.
IC_field_value['v'] = 0.
IC_field_value['w'] = 0.
IC_field_value['p'] = 0.
IC_field_value['turb_k'] = 0.001
IC_field_value['turb_e'] = 0.001

turbulence = {}
turbulence['model'] = 'laminar' # choose from: 'laminar','k-omega'
turbulence['omega_wall'] = 0.001

# fluid properties

fluid = {}
fluid['water'] = { 'rho':998.2, 'mu':0.001, 'phi':0. }
fluid['air']   = { 'rho':1.205, 'mu':1.807e-5, 'phi':1. }

phase = ['water','air']
if case['num_phases'] == 1: phase = ['water']

nphase = len(phase)

# Gravity

gx = 0.0; gy = 0.0; gz = 0.0;

if case['num_phases'] > 1: 
  if case['gravity_dir'] == '+x': gx =  9.81
  if case['gravity_dir'] == '-x': gx = -9.81
  if case['gravity_dir'] == '+y': gy =  9.81
  if case['gravity_dir'] == '-y': gy = -9.81
  if case['gravity_dir'] == '+z': gz =  9.81
  if case['gravity_dir'] == '-z': gz = -9.81

gravity = [gx,gy,gz]

if 'x' in case['gravity_dir']: gravity_axis = 0
if 'y' in case['gravity_dir']: gravity_axis = 1
if 'z' in case['gravity_dir']: gravity_axis = 2

# input options    

spaceOrder = 1
useRANS = 0      # 0-None  1-K-Epsilon  2-K-Omega
if turbulence['model'] == 'laminar': useRANS = 0
if turbulence['model'] == 'k-omega': useRANS = 2

# initial conditions

waterLine = case['waterline']

def signedDistance(x):
    sd = x[gravity_axis] - waterLine 
    return sd
    
def IC_signed_distance(x):
    phi = x[gravity_axis] - waterLine 
    return phi

class IC_field_constant:
    def __init__(self,value): self.value = value
    def uOfXT(self,x,t): return self.value
    def uOfX(self,x):    return self.value

class IC_field_p:
    def __init__(self): pass
    def uOfXT(self,x,t):
        rho_air     = fluid['air']['rho']
        rho_water   = fluid['water']['rho']
        g_component = gravity[gravity_axis]
        waterLevel  = waterLine
        p_air   = rho_air    * g_component * ( default_p.L[gravity_axis] - waterLevel )
        p_water = rho_water  * g_component * ( waterLevel - x[gravity_axis] )
        p_hydrostatic = p_air
        if IC_signed_distance(x) < 0: p_hydrostatic = p_water
        return -p_hydrostatic

class IC_field_phi:       
    def uOfXT(self,x,t): return IC_signed_distance(x)

# mesh and BCs 

filename = mesh['filename']
mesh_nominal_spacing = mesh['nominal_spacing']

u = IC_field_value['u']
v = IC_field_value['v']
w = IC_field_value['w']
Vmag = math.sqrt( u*u + v*v + w*w )

ele_fluid  = { 'type':'fluid' }
ele_zone = {}
ele_zone['fluid'] = { 'meshtag': 1, 'condition': ele_fluid }

bc_wall   = { 'method': std_bc.noSlip }
bc_slip   = { 'method': std_bc.freeSlip }
bc_inlet  = { 'method': std_bc.velocityInlet_rans2p, 
              'Vmag': Vmag,
              'sdf': IC_signed_distance,
              'fluid': fluid,
              'smoothing': mesh_nominal_spacing,
              'epsFact': 1.5,
}
bc_outlet = { 'method': std_bc.outflow_rans2p, 
              'sdf': IC_signed_distance,
              'fluid': fluid,
              'gravity': gravity,
              'gravity_axis': gravity_axis,
              'waterLevel': waterLine,
              'smoothing': mesh_nominal_spacing,
              'epsFact': 1.5,
}
bc_open   = { 'method': std_bc.open,
              'sdf': IC_signed_distance,
              'fluid': fluid,
              'gravity': gravity,
              'gravity_axis': gravity_axis,
              'waterLevel': waterLine,
              'smoothing': mesh_nominal_spacing,
              'epsFact': 1.5,
}
bc_interior = { 'method': std_bc.interior }

bc_zone = {}
bc_zone['int']  = { 'meshtag': 0,  'condition': bc_interior }
bc_zone['io-i'] = { 'meshtag': 4,  'condition': bc_inlet }
bc_zone['io-o'] = { 'meshtag': 6,  'condition': bc_outlet }
bc_zone['io-l'] = { 'meshtag': 5,  'condition': bc_wall }
bc_zone['io-r'] = { 'meshtag': 7,  'condition': bc_wall }
bc_zone['io-t'] = { 'meshtag': 8,  'condition': bc_open }
bc_zone['sw-r'] = { 'meshtag': 9,  'condition': bc_wall }
bc_zone['sw-b'] = { 'meshtag': 12, 'condition': bc_wall }
bc_zone['sw-w'] = { 'meshtag': 13, 'condition': bc_wall } 

# knobs that control system convergence
# change from default values via: knob['vof']['sc_fac'] = 0.1

knob = std_knob.set_defaults()  
knob['vof']['sc_fac']     = 0.1
knob['rd']['sc_fac']      = 0.1
knob['ns']['nl_atol']     = 1.e-4
knob['vof']['nl_atol']    = 1.e-4
knob['ls']['nl_atol']     = 1.e-4
knob['rd']['nl_atol']     = 1.e-3
knob['mcorr']['nl_atol']  = 1.e-4
knob['turb_k']['nl_atol'] = 1.e-4
knob['turb_e']['nl_atol'] = 1.e-4
