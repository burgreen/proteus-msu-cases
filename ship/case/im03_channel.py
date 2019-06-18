# im03_channel.py

import math
from proteus import default_p
from . import std_bc as std

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

dt_init  = -0.0001
dt_fixed = 0.0001
dt_fixed_steps = 2 

# initial conditions

waterLine_z = 0.0
#outlet_waterLine_z = 0.0

def IC_signed_distance(x):
    phi_z = x[2]-waterLine_z 
    #phi_y = x[1]-(-180.)
    return phi_z

def signedDistance(x):
    phi_z = x[1] - waterLine_z
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
IC_field_value['u'] = 0.
IC_field_value['v'] = -10.0
IC_field_value['w'] = 0.

# mesh and BCs 

mesh_nominal_spacing = 0.05

u = IC_field_value['u']
v = IC_field_value['v']
w = IC_field_value['w']
Vmag = math.sqrt( u*u + v*v + w*w )

bc_wall   = { 'type':'noSlip' }
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
              'waterLevel': waterLine_z,
              'smoothing': mesh_nominal_spacing,
}
bc_open = { 'type':'open' }
bc_interior = { 'type':'interior' }

bc_zone = {}

filename = 'mesh/im03a.tetgen'

bc_zone['interior'] = { 'meshtag': 0,  'condition': bc_interior, 'custom':std.interior }
bc_zone['channel']  = { 'meshtag': 5,  'condition': bc_wall,     'custom':std.noSlip }
bc_zone['inlet']    = { 'meshtag': 1,  'condition': bc_inlet }
bc_zone['outlet']   = { 'meshtag': 2,  'condition': bc_outlet,   'custom':std.outflow_rans2p }
bc_zone['side1']    = { 'meshtag': 6,  'condition': bc_wall,     'custom':std.noSlip }
bc_zone['side2']    = { 'meshtag': 7,  'condition': bc_wall,     'custom':std.noSlip }
bc_zone['ground']   = { 'meshtag': 4,  'condition': bc_wall,     'custom':std.noSlip }
bc_zone['air']      = { 'meshtag': 3,  'condition': bc_open,     'custom':std.open }

ele_fluid  = { 'type':'fluid' }
  
ele_zone = {}
ele_zone['fluid'] = { 'meshtag': 1, 'condition': ele_fluid }

tols = {}
'''
tols['nl_atol_rans2p'] = 1.e-5
tols['nl_atol_vof'   ] = 1.e-5
tols['nl_atol_ls'    ] = 1.e-0
tols['nl_atol_redist'] = 1.e-2
tols['nl_atol_mcorr' ] = 1.e-5
tols['nl_atol_turb_k'] = 1.e-5
tols['nl_atol_turb_e'] = 1.e-5
'''
tols['nl_atol_rans2p'] = 1.e-4
tols['nl_atol_vof'   ] = 1.e-4
tols['nl_atol_ls'    ] = 1.e-4
tols['nl_atol_redist'] = 1.e-3 # a lot looser
tols['nl_atol_mcorr' ] = 1.e-4
tols['nl_atol_turb_k'] = 1.e-4
tols['nl_atol_turb_e'] = 1.e-4
