# visc-channel.py

import math
from . import std_bc as std_bc

# mesh and BCs

# fluid properties

fluid = {}
fluid['water'] = { 'rho':998.2, 'mu':0.001 }
fluid['air']   = { 'rho':1.205, 'mu':1.807e-5 }

# phase = ['water','air']
phase = ['water']

nphase = len(phase)

# Gravity

gravity = [0.,0.,0.]

if len(phase) > 1: gravity = [0.,0.,-9.8] 

# input options    

spaceOrder = 1
useRANS    = 0      # 0 -- None # 1 -- K-Epsilon # 2 -- K-Omega

# time stepping

dt_init  = 0.01
dt_fixed = 0.01
dt_fixed_steps = 1

# initial conditions

class IC_field_constant:
    def __init__(self,value): self.value = value
    def uOfXT(self,x,t): return self.value
    def uOfX(self,x):    return self.value

class IC_field_zero:
    def __init__(self): pass
    def uOfXT(self,x,t): return 0.0

IC_field_value = {}
IC_field_value['p'] = 0.
IC_field_value['u'] = 1.
IC_field_value['v'] = 0.
IC_field_value['w'] = 0.

u = IC_field_value['u']
v = IC_field_value['v']
w = IC_field_value['w']
Vmag = math.sqrt( u*u + v*v + w*w )

# mesh and bcs

filename = 'mesh/visc-0_01.tetgen'

mesh_nominal_spacing = 0.05

bc_noslip   = { 'method': std_bc.noSlip }
bc_slip     = { 'method': std_bc.freeSlip }
bc_inlet    = { 'method': std_bc.velocityInlet, 'Vmag':Vmag }
bc_outlet   = { 'method': std_bc.outflow, 'p':IC_field_value['p'] }
bc_interior = { 'method': std_bc.interior }
  
bc_zone = {}
bc_zone['x+']  = { 'meshtag': 1,  'condition': bc_inlet  }
bc_zone['x-']  = { 'meshtag': 2,  'condition': bc_outlet }
bc_zone['y-']  = { 'meshtag': 3,  'condition': bc_noslip }
bc_zone['y+']  = { 'meshtag': 4,  'condition': bc_slip }
bc_zone['z-']  = { 'meshtag': 5,  'condition': bc_slip }
bc_zone['z+']  = { 'meshtag': 6,  'condition': bc_slip }

ele_fluid  = { 'type':'fluid' }

ele_zone = {}
ele_zone['fluid'] = { 'meshtag': 1, 'condition': ele_fluid }

tols = {}