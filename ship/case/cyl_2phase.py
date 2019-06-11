# user_param.py  case: marin_2_phase treated as open duct

import math
from proteus import default_p
from . import cyl_2phase_bc as custom

# fluid properties

fluid = {}
fluid['water'] = { 'rho':998.2, 'mu':0.001, 'phi':0. }
fluid['air']   = { 'rho':1.205, 'mu':1.807e-5, 'phi':1. }

phase = ['water','air']

nphase = len(phase)

# Gravity

gravity = [0.,0.,0.]

if nphase > 1: gravity = [0.,-9.8,0.] 

# input options    

spaceOrder = 1
useRANS    = 0      # 0 -- None # 1 -- K-Epsilon # 2 -- K-Omega

# time stepping

dt_init  = 0.0001
dt_fixed = 0.0001
dt_fixed_steps = 0

# initial conditions

waterLine_y = 1.75

def IC_signed_distance(x):
    phi_y = x[1]-waterLine_y 
    return phi_y

class IC_field_constant:
    def __init__(self,value): self.value = value
    def uOfXT(self,x,t): return self.value
    def uOfX(self,x):    return self.value

class IC_field_p:
    def __init__(self): pass
    def uOfXT(self,x,t):
        rho_air     = fluid['air']['rho']
        rho_water   = fluid['water']['rho']
        g_component = gravity[1]
        waterLevel  = waterLine_y
        p_air   = rho_air    * g_component * ( default_p.L[1] - waterLevel )
        p_water = rho_water  * g_component * ( waterLevel - x[1] )
        p_hydrostatic = p_air
        if IC_signed_distance(x) < 0: p_hydrostatic = p_water
        return -p_hydrostatic

class IC_field_phi:       
    def __init__(self): pass
    def uOfXT(self,x,t): return IC_signed_distance(x)

IC_field_value = {}
IC_field_value['p'] = 0.
IC_field_value['u'] = 1.
IC_field_value['v'] = 0.
IC_field_value['w'] = 0.

# mesh and BCs 

filename = 'mesh/cyl_01.tetgen'

mesh_nominal_spacing = 0.05

u = IC_field_value['u']
v = IC_field_value['v']
w = IC_field_value['w']
Vmag = math.sqrt( u*u + v*v + w*w )

bc_wall   = { 'type':'noSlip' }
bc_slip   = { 'type':'freeSlip' }
bc_inlet  = { 'type':'velocityInlet_rans2p', 
              'Vmag': Vmag,
              'sdf': IC_signed_distance,
              'fluid': fluid,
              'smoothing': mesh_nominal_spacing,
}
bc_outlet = { 'type':'outflow_rans2p', 
              'sdf': IC_signed_distance,
              'fluid': fluid,
              'gravity': gravity,
              'waterLevel': waterLine_y,
              'smoothing': mesh_nominal_spacing,
}
bc_open = { 'type':'open' }
bc_interior = { 'type':'interior' }

bc_zone = {}
bc_zone['interior'] = { 'meshtag': 0,  'condition': bc_interior }
bc_zone['cyl']      = { 'meshtag': 7,  'condition': bc_wall }
bc_zone['x0']       = { 'meshtag': 1,  'condition': bc_inlet }
bc_zone['x1']       = { 'meshtag': 2,  'condition': bc_outlet, 'custom': custom.outflow_rans2p }
bc_zone['y0']       = { 'meshtag': 3,  'condition': bc_slip }
bc_zone['y1']       = { 'meshtag': 4,  'condition': bc_slip }
bc_zone['z0']       = { 'meshtag': 5,  'condition': bc_slip }
bc_zone['z1']       = { 'meshtag': 6,  'condition': bc_slip }
  
ele_fluid  = { 'type':'fluid' }
  
ele_zone = {}
ele_zone['fluid'] = { 'meshtag': 1, 'condition': ele_fluid }

tols = {}
