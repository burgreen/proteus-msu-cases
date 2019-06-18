# user_param.py  case: marin_2_phase treated as open duct

import math
from proteus import default_p

# fluid properties

fluid = {}
fluid['water'] = { 'rho':998.2, 'mu':0.001, 'phi':0. }
fluid['air']   = { 'rho':1.205, 'mu':1.807e-5, 'phi':1. }

phase = ['water','air']

nphase = len(phase)

# Gravity

gravity = [0.,0.,0.]

if nphase > 1: gravity = [0.,0.,-9.81] 

# input options    

spaceOrder = 1
useRANS    = 0      # 0 -- None # 1 -- k-e # 2 -- K-Omega

# time stepping

dt_init  = 0.001
dt_fixed = 0.001
dt_fixed_steps = 0

# initial conditions

waterLine_z = 0.5

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

# mesh and BCs 

filename = 'mesh/mesh_wigley.tetgen'

mesh_nominal_spacing = 0.05

hull_length = 1.0
hull_beam   = hull_length/10.0
hull_draft  = hull_length/16.0
hull_center = (0.0, 0.0, 0.5*hull_length)
waterLevel  = 0.5*hull_length

Fr = 0.25 # 0.51  0.0
Um = Fr*math.sqrt(math.fabs(gravity[2])*hull_length)
Re = hull_length*Um*fluid['water']['rho']/fluid['water']['mu']

residence_time = 0.0
if Um > 0.0: residence_time = hull_length/Um

IC_field_value = {}
IC_field_value['p'] = 0.
IC_field_value['u'] = Um
IC_field_value['v'] = 0.
IC_field_value['w'] = 0.
IC_field_value['turb_k'] = 0.001
IC_field_value['turb_e'] = 0.001

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
              'turb_k': IC_field_value['turb_k'],
              'turb_e': IC_field_value['turb_e'],
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
bc_zone['interior'] = { 'meshtag': 0,  'condition': bc_interior }
bc_zone['ship']     = { 'meshtag': 7,  'condition': bc_wall }
bc_zone['inlet']    = { 'meshtag': 5,  'condition': bc_inlet }
bc_zone['outlet']   = { 'meshtag': 3,  'condition': bc_outlet }
bc_zone['y0']       = { 'meshtag': 2,  'condition': bc_slip }
bc_zone['y1']       = { 'meshtag': 4,  'condition': bc_slip }
bc_zone['bot']      = { 'meshtag': 1,  'condition': bc_slip }
bc_zone['air']      = { 'meshtag': 6,  'condition': bc_open }
  
ele_fluid  = { 'type':'fluid' }
  
ele_zone = {}
ele_zone['fluid'] = { 'meshtag': 1, 'condition': ele_fluid }

tols = {}