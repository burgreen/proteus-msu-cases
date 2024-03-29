# ip02.py srf-piercing-cylinder

import math
from proteus import default_p
from . import std_knob as std_knob
from . import std_bc as std_bc

# time stepping

dt_init  = 0.0001
dt_fixed = 0.001
dt_fixed_steps = 3 

# case/geometry properties

mesh = {}
mesh['filename'] = 'mesh/ip02.tetgen'
mesh['nominal_spacing'] = 0.05

case = {}
case['num_phases']  = 2
case['gravity_dir'] = '-z' # choose from: '+x','-x','+y','-y','+z','-z'
case['waterline']   = 0.0

IC_field_value = {} # initial conditions
IC_field_value['u'] = -1.
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

def IC_signed_distance(x):
    phi = x[gravity_axis] - waterLine 
    return phi

def signedDistance(x):
    sd = x[gravity_axis] - waterLine 
    return sd
    ''' 
    https://github.com/erdc/air-water-vv/blob/master/2d/hydraulicStructures/broad_crested_weir/broad_crested_weir.py
    phi_x = x[0] - waterLine_x
    phi_z = x[1] - waterLine_z
    phi_z_outflow = x[1] - outflow_level
    if phi_x <= 0.0:
        if phi_z < 0.0:
            return max(phi_x, phi_z)
        else:
            return phi_z
    else:
        if phi_z_outflow < 0.0:
            return phi_z_outflow
        else:
            if phi_z < 0.0:
                return min(phi_x, phi_z_outflow)
            else:
    return min(sqrt(phi_x ** 2 + phi_z ** 2), phi_z_outflow)
    '''

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
    def __init__(self): pass
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
bc_zone['interior'] = { 'meshtag': 0,  'condition': bc_interior }
bc_zone['inlet']    = { 'meshtag': 4,  'condition': bc_inlet }
bc_zone['outlet']   = { 'meshtag': 5,  'condition': bc_outlet }
bc_zone['cyl']      = { 'meshtag': 6,  'condition': bc_wall }
bc_zone['y0']       = { 'meshtag': 7,  'condition': bc_wall }
bc_zone['y1']       = { 'meshtag': 8,  'condition': bc_wall }
bc_zone['z0']       = { 'meshtag': 9,  'condition': bc_wall }
bc_zone['z1']       = { 'meshtag': 10, 'condition': bc_wall }

# knobs that control system convergence
# change from default values via: knob['vof']['sc_fac'] = 0.1

knob = std_knob.set_defaults()  

knob['vof']['sc_fac']     = 0.1
knob['rd']['sc_fac']      = 0.1
knob['ls']['sc_fac']      = 0.1

knob['ns']['nl_atol']     = 1.e-4
knob['vof']['nl_atol']    = 1.e-4
knob['ls']['nl_atol']     = 1.e-4
knob['rd']['nl_atol']     = 1.e-3
knob['mcorr']['nl_atol']  = 1.e-4
knob['turb_k']['nl_atol'] = 1.e-4
knob['turb_e']['nl_atol'] = 1.e-4
