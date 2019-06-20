# cyl_1p.py

import math
from proteus import default_p
from . import std_bc as std_bc

# fluid properties

fluid = {}
fluid['water'] = { 'rho':998.2, 'mu':0.001, 'phi':0. }
fluid['air']   = { 'rho':1.205, 'mu':1.807e-5, 'phi':1. }

phase = ['water']

nphase = len(phase)

# Gravity

gravity = [0.,0.,0.]

if nphase > 1: gravity = [0.,0.,-9.8] 

# input options    

spaceOrder = 1
useRANS    = 0      # 0 -- None # 1 -- K-Epsilon # 2 -- K-Omega

# time stepping

dt_init  = 0.01
dt_fixed = 0.01
dt_fixed_steps = 2

# initial conditions

waterLine_z = 0.0

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
IC_field_value['u'] = 1.
IC_field_value['v'] = 0.
IC_field_value['w'] = 0.

# mesh and BCs 

filename = 'mesh/cyl_02.tetgen'

mesh_nominal_spacing = 0.05

u = IC_field_value['u']
v = IC_field_value['v']
w = IC_field_value['w']
Vmag = math.sqrt( u*u + v*v + w*w )

bc_wall   = { 'method': std_bc.noSlip }
bc_slip   = { 'method': std_bc.freeSlip }
bc_inlet  = { 'method': std_bc.velocityInlet, 
              'Vmag': Vmag,
              'sdf': IC_signed_distance,
              'fluid': fluid,
              'smoothing': mesh_nominal_spacing,
}
bc_outlet = { 'method': std_bc.outflow, 
              'sdf': IC_signed_distance,
              'fluid': fluid,
              'gravity': gravity,
              'waterLevel': waterLine_z,
              'smoothing': mesh_nominal_spacing,
}
bc_open = { 'method': std_bc.open }
bc_interior = { 'method': std_bc.interior }

bc_zone = {}
#bc_zone['interior'] = { 'meshtag': 0,  'condition': bc_interior }
bc_zone['cyl']      = { 'meshtag': 0,  'condition': bc_wall }
bc_zone['x0']       = { 'meshtag': 1,  'condition': bc_inlet }
bc_zone['x1']       = { 'meshtag': 2,  'condition': bc_outlet }
bc_zone['y0']       = { 'meshtag': 3,  'condition': bc_slip }
bc_zone['y1']       = { 'meshtag': 4,  'condition': bc_slip }
bc_zone['z0']       = { 'meshtag': 5,  'condition': bc_slip }
bc_zone['z1']       = { 'meshtag': 6,  'condition': bc_slip }
  
ele_fluid  = { 'type':'fluid' }
  
ele_zone = {}
ele_zone['fluid'] = { 'meshtag': 1, 'condition': ele_fluid }

tols = {}