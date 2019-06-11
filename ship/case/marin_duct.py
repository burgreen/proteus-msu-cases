# user_param.py  case: marin_2_phase treated as open duct

from math    import *
from proteus import default_p

# fluid properties

fluid = {}
fluid['water'] = { 'rho':998.2, 'mu':0.001, 'phi':0. }
fluid['air']   = { 'rho':1.205, 'mu':1.807e-5, 'phi':1. }

phase = ['water','air']

nphase = len(phase)

# Gravity

gravity = [0.,0.,0.]

if nphase > 1: gravity = [0.,0.,-9.8] 

# input options    

spaceOrder = 1
useRANS    = 0  # 0=None # 1=K-Epsilon # 2=K-Omega

# time stepping

dt_init  = 0.01
dt_fixed = 0.01
dt_fixed_steps = 20 

# initial conditions

waterLine_z = 0.3

def IC_signed_distance(x):
    phi_z = x[2]-waterLine_z 
    return phi_z

class IC_field_constant:
    def __init__(self,value): self.value = value
    def uOfXT(self,x,t): return self.value
    def uOfX(self,x):    return self.value

class IC_field_p:
    def __init__(self): pass
    def uOfXT(self,x,t):
        rho_air     = fluid['air']['rho']
        rho_water   = fluid['water']['rho']
        g_component = gravity[2]
        waterLevel  = waterLine_z
        p_air   = rho_air    * g_component * ( default_p.L[2] - waterLevel )
        p_water = rho_water  * g_component * ( waterLevel - x[2] )
        p_hydrostatic = p_air
        if IC_signed_distance(x) < 0: p_hydrostatic = p_water
        return -p_hydrostatic

class IC_field_phi:       
    def __init__(self): pass
    def uOfXT(self,x,t): return IC_signed_distance(x)

IC_field_value = {}
IC_field_value['p'] = 0.
IC_field_value['u'] = -1.
IC_field_value['v'] = 0.
IC_field_value['w'] = 0.

# mesh and BCs 

filename = 'mesh/marin_duct.tetgen'

mesh_nominal_spacing = 0.05

bc_wall   = { 'type':'NoSlip' }
bc_inlet  = { 'type':'velocityInlet_rans2p', 
              'Vmag':1.0,
              'sdf': IC_signed_distance,
              'fluid': fluid,
              'smoothing': mesh_nominal_spacing,
}
bc_outlet = { 'type':'outflow_rans2p', 
              'sdf': IC_signed_distance,
              'fluid': fluid,
              'gravity': gravity,
              'waterLevel': waterLine_z,
              'smoothing': mesh_nominal_spacing,
}
bc_open = { 'type':'open' }
bc_interior = { 'type':'interior' }

bc_zone = {}
bc_zone['x+']     = { 'meshtag': 2,  'condition': bc_inlet }
bc_zone['x-']     = { 'meshtag': 1,  'condition': bc_outlet }
bc_zone['z-']     = { 'meshtag': 3,  'condition': bc_wall }
bc_zone['z+']     = { 'meshtag': 4,  'condition': bc_open }
bc_zone['y-']     = { 'meshtag': 5,  'condition': bc_wall }
bc_zone['y+']     = { 'meshtag': 6,  'condition': bc_wall }
bc_zone['box_x-'] = { 'meshtag': 7,  'condition': bc_wall }
bc_zone['box_x+'] = { 'meshtag': 8,  'condition': bc_wall }
bc_zone['box_z+'] = { 'meshtag': 9,  'condition': bc_wall }
bc_zone['box_y-'] = { 'meshtag': 10, 'condition': bc_wall }
bc_zone['box_y+'] = { 'meshtag': 11, 'condition': bc_wall }
bc_zone['interior'] = { 'meshtag': 0, 'condition': bc_interior }
  
ele_fluid  = { 'type':'fluid' }
  
ele_zone = {}
ele_zone['fluid'] = { 'meshtag': 1, 'condition': ele_fluid }

tols = {}