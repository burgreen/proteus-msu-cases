# user_param.py  case: hw03 vehicle situated midstream

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
useRANS    = 0      # 0 -- None # 1 -- K-Epsilon # 2 -- K-Omega

# time stepping

dt_fixed = 0.0001
dt_fixed_steps = 1 

# initial conditions

waterLine_z = 4.0
waterLine_z = 9.0
waterLine_z = 2.0

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
    def uOfXT(self,x,t): return IC_signed_distance(x)

IC_field_value = {}
IC_field_value['p'] = 0.
IC_field_value['u'] = 1.
IC_field_value['v'] = 0.
IC_field_value['w'] = 0.

# mesh and BCs 

filename = 'mesh/pw-hw03.tetgen'

mesh_nominal_spacing = 0.05

bc_wall   = { 'type':'FreeSlip' }
bc_open   = { 'type':'open' }
bc_inlet  = { 'type':'velocityInlet_rans2p', 
              'Vmag':1.0,
              'sdf': IC_signed_distance,
              'fluid': fluid
}
bc_outlet = { 'type':'outflow_rans2p', 
              'sdf': IC_signed_distance,
              'fluid': fluid,
              'gravity': gravity,
              'waterLevel': waterLine_z
}

bc_zone = {}
bc_zone['io-i'] = { 'meshtag': 4,  'condition': bc_inlet }
bc_zone['io-o'] = { 'meshtag': 6,  'condition': bc_outlet }
bc_zone['io-l'] = { 'meshtag': 5,  'condition': bc_wall }
bc_zone['io-r'] = { 'meshtag': 7,  'condition': bc_wall }
bc_zone['io-t'] = { 'meshtag': 8,  'condition': bc_open }
bc_zone['sw-r'] = { 'meshtag': 9,  'condition': bc_wall }
bc_zone['sw-b'] = { 'meshtag': 12, 'condition': bc_wall }
bc_zone['sw-w'] = { 'meshtag': 13, 'condition': bc_wall } 

ele_fluid  = { 'type':'fluid' }
ele_zone = {}
ele_zone['fluid'] = { 'meshtag': 1, 'condition': ele_fluid }
