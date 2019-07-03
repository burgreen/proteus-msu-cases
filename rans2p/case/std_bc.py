# case/std_bc.py

import sys
#import numpy as np
#from proteus.Profiling import logEvent
#from proteus import MeshTools
#from proteus.mprans import BoundaryConditions as bc
#from proteus import Domain
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus import default_p

def velocityInlet_rans2p( bc, condition ):
    """
    Sets velocity inlet bc for twoPhase flows
    """
    bc.reset()

    normal  = condition['_normal']
    epsFact = condition['epsFact']
    he      = condition['smoothing']

    Vmag = 0.
    if 'Vmag' in condition: Vmag = condition['Vmag']
    Vmag_air   = Vmag*1
    Vmag_water = Vmag

    turb_k = 0.
    turb_e = 0.
    if 'turb_k' in condition: turb_k = condition['turb_k']
    if 'turb_e' in condition: turb_k = condition['turb_e']

    sdf   = condition['sdf']
    fluid = condition['fluid']

    phi_air   = fluid['air']['phi']
    phi_water = fluid['water']['phi']

    def vel_dirichlet(i):
            def dirichlet(x,t):
                phi = sdf(x)
                #if phi <= 0.:
                #    H = 0.0
                #elif 0 < phi <= smoothing:
                #    H = smoothedHeaviside(smoothing / 2., phi - smoothing / 2.)
                #else:
                #    H = 1.0
                H = 1.
                if phi <= 0.: H = 0.
                vel = H * Vmag_air*normal[i] + (1-H) * Vmag_water*normal[i]
                return vel
            return dirichlet

    def vof_dirichlet(x,t):
            phi = sdf(x)
            #if phi >= smoothing:
            #    H = 1.
            #elif smoothing > 0 and -smoothing < phi < smoothing:
            #    H = smoothedHeaviside(smoothing, phi)
            #elif phi <= -smoothing:
            #    H = 0.
            #H = 1.
            #if phi <= 0.: H = 0.
            #vof = H * phi_air + (1-H) * phi_water
            vof = smoothedHeaviside( epsFact*he, phi )
            return vof

    def p_advective(x,t):
            phi = sdf(x)
            #if phi <= 0.:
            #    H = 0.0
            #elif 0 < phi <= smoothing:
            #    H = smoothedHeaviside(smoothing / 2., phi - smoothing / 2.)
            #else:
            #    H = 1.0
            H = 1.
            if phi <= 0.: H = 0.
            u = H * Vmag_air + (1-H) * Vmag_water
            # This is the normal velocity, based on the inwards boundary
            # orientation -b_or
            #u_p = np.sum(u * np.abs(b_or))
            #return -u_p
            return -u

    def ls_dirichlet(x,t):
            phi = sdf(x)
            return phi

    def constant(c):
        def function(x,t): return c
        return function

    bc.u_dirichlet.uOfXT = vel_dirichlet(0)
    bc.u_advective.uOfXT = None
    bc.u_diffusive.uOfXT = None

    bc.v_dirichlet.uOfXT = vel_dirichlet(1)
    bc.w_advective.uOfXT = None
    bc.w_diffusive.uOfXT = None

    bc.w_dirichlet.uOfXT = vel_dirichlet(2)
    bc.w_advective.uOfXT = None
    bc.w_diffusive.uOfXT = None

    if 'u_xt' in condition: bc.u_dirichlet.uOfXT = condition['u_xt']
    if 'v_xt' in condition: bc.v_dirichlet.uOfXT = condition['v_xt']
    if 'w_xt' in condition: bc.w_dirichlet.uOfXT = condition['w_xt']

    bc.clsvof_dirichlet.uOfXT = ls_dirichlet

    bc.k_dirichlet.uOfXT = constant(turb_k)
    #bc.k_advective.
    #bc.k_diffusive.

    bc.dissipation_dirichlet.uOfXT = constant(turb_e)
    #bc.dissipation_advective.
    #bc.dissipation_diffusive.

    bc.vof_dirichlet.uOfXT = vof_dirichlet
    bc.vof_advective.uOfXT = None
    ##.vof_diffusive.uOfXT = does not exist

    bc.p_dirichlet.uOfXT = None
    bc.p_advective.uOfXT = p_advective
    ##.p_diffusive.uOfXT = does not exist

    if 'flux_xt' in condition: bc.p_advective.uOfXT = condition['flux_xt']
    
def outflow_rans2p( bc, condition ):
    """
    Sets outflow bc for twoPhase flows
    """
    bc.reset()

    sdf        = condition['sdf']
    fluid      = condition['fluid']
    gravity    = condition['gravity']
    waterLevel = condition['waterLevel']
    axis       = condition['gravity_axis']
    epsFact    = condition['epsFact']
    he         = condition['smoothing']
    fac        = 1.0 

    rho_air   = fluid['air']['rho']
    rho_water = fluid['water']['rho']
    mu_air    = fluid['air']['mu']
    mu_water  = fluid['water']['mu']
    phi_air   = fluid['air']['phi']
    phi_water = fluid['water']['phi']

    def old_p_dirichlet(x,t):  # assumes a constant fixed fluid height at the outlet
        g_component = gravity[2]
        p_air   = rho_air    * g_component * ( default_p.L[2] - waterLevel )
        p_water = rho_water  * g_component * ( waterLevel - x[2] )
        p_hydrostatic = p_air
        #if sdf(x) < 0: p_hydrostatic = p_water
        if x[2] < waterLevel: p_hydrostatic = p_water
        return -p_hydrostatic

    def p_dirichlet(x,t):
      p_L = default_p.L[axis]*rho_air*gravity[axis]
      phi_L = default_p.L[axis] - waterLevel
      phi = x[axis] - waterLevel
      return p_L - gravity[axis]*(  rho_water*(phi_L - phi)
                                  + (rho_air -rho_water)*(smoothedHeaviside_integral(epsFact*he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact*he,phi)))

    def u_diffusive(x,t):
        g = p_dirichlet(x,t)
        phi = sdf(x)
        if phi >= smoothing: H = 1.
        elif -smoothing < phi < smoothing: H = smoothedHeaviside(smoothing, phi)
        elif phi <= -smoothing: H = 0.
        #H = 1.
        #if phi <= 0.: H = 0.
        diff =  H*(mu_air)*g + (1-H)*(mu_water)*g
        return fac*diff

    def vof_dirichlet(x,t):
        phi = sdf(x)
        #if phi >= smoothing: H = 1.
        #elif -smoothing < phi < smoothing: H = smoothedHeaviside(smoothing, phi)
        #elif phi <= -smoothing: H = 0.
        #H = 1.
        #if phi <= 0.: H = 0.
        #return H * phi_air + (1-H) * phi_water
        vof = smoothedHeaviside( epsFact*he, phi )
        return vof

    def constant(c):
        def function(x,t): return c
        return function

    bc.u_dirichlet.uOfXT = None # constant(0.)
    bc.u_advective.uOfXT = None
    bc.u_diffusive.uOfXT = constant(0.)

    bc.v_dirichlet.uOfXT = None # constant(0.)
    bc.v_advective.uOfXT = None
    bc.v_diffusive.uOfXT = constant(0.)

    bc.w_dirichlet.uOfXT = None # constant(0.)
    bc.w_advective.uOfXT = None
    bc.w_diffusive.uOfXT = constant(0.)

    bc.k_dirichlet.uOfXT = constant(0.)
    bc.k_advective.uOfXT = None
    bc.k_diffusive.uOfXT = constant(0.)

    bc.clsvof_dirichlet.uOfXT = None

    bc.dissipation_dirichlet.uOfXT = constant(0.)
    bc.dissipation_advective.uOfXT = None
    bc.dissipation_diffusive.uOfXT = constant(0.)

    bc.vof_dirichlet.uOfXT = None # vof_dirichlet
    bc.vof_advective.uOfXT = None
    ##.vof_diffusive.uOfXT = does not exist

    bc.p_dirichlet.uOfXT = p_dirichlet
    bc.p_advective.uOfXT = None
    ##.p_diffusive.uOfXT = does not exist

    '''
    if U is not None:
            def get_inlet_ux_dirichlet(i):
                def ux_dirichlet(x, t):
                    phi = x[vert_axis] - seaLevel
                    if phi <= 0.:
                        H = 0.0
                    elif 0 < phi <= smoothing:
                        H = smoothedHeaviside(smoothing / 2., phi - smoothing / 2.)
                    else:
                        H = 1.0
                    return H * Uwind[i] + (1 - H) * U[i]
                return ux_dirichlet

            if Uwind is None:
                Uwind = np.zeros(3)
            U = np.array(U)
            Uwind = np.array(Uwind)
            self.u_dirichlet.uOfXT = get_inlet_ux_dirichlet(0)
            self.v_dirichlet.uOfXT = get_inlet_ux_dirichlet(1)
            self.w_dirichlet.uOfXT = get_inlet_ux_dirichlet(2)
            self.u_diffusive.resetBC()
    '''

def freeSlip( bc, condition ):
    '''
    Sets free slip conditions at the boundary
    '''
    bc.reset()

    def constant(c):
        def function(x,t): return c
        return function

    bc.u_dirichlet.uOfXT = None
    bc.u_advective.uOfXT = constant(0.)
    bc.u_diffusive.uOfXT = constant(0.)

    bc.v_dirichlet.uOfXT = None
    bc.v_advective.uOfXT = constant(0.)
    bc.v_diffusive.uOfXT = constant(0.)

    bc.w_dirichlet.uOfXT = None
    bc.w_advective.uOfXT = constant(0.)
    bc.w_diffusive.uOfXT = constant(0.)

    bc.clsvof_dirichlet.uOfXT = None

    bc.vof_dirichlet.uOfXT = None
    bc.vof_advective.setConstantBC(0.)
    ##.vof_diffusive.uOfXT = does not exist

    bc.k_dirichlet.setConstantBC(0.) 
    #bc.k_advective.
    bc.k_diffusive.setConstantBC(0.)

    #bc.dissipation_dirichlet.
    #bc.dissipation_advective.
    bc.dissipation_diffusive.setConstantBC(0.)  
    
    bc.p_dirichlet.uOfXT = None
    bc.p_advective.setConstantBC(0.)
    #bc.p_diffusive.uOfXT = None

    bc.pInc_dirichlet.uOfXT = None
    bc.pInc_advective.setConstantBC(0.)  
    bc.pInc_diffusive.setConstantBC(0.)

    bc.pInit_dirichlet.uOfXT = None
    bc.pInit_advective.setConstantBC(0.)
    bc.pInit_diffusive.setConstantBC(0.)

    # dirichlet
    # advective        
    bc.us_advective.setConstantBC(0.)
    bc.vs_advective.setConstantBC(0.)
    bc.ws_advective.setConstantBC(0.)
    bc.vos_advective.setConstantBC(0.)
    # diffusive
    bc.us_diffusive.setConstantBC(0.)
    bc.vs_diffusive.setConstantBC(0.)
    bc.ws_diffusive.setConstantBC(0.)

def interior( bc, condition ):
    """
    Sets interior bc for rans3p flows
    """
    bc.reset()

    def constant(c):
        def function(x,t): return c
        return function

    bc.p_dirichlet.uOfXT = None

    bc.u_diffusive.setConstantBC(0.)
    bc.v_diffusive.setConstantBC(0.)
    bc.w_diffusive.setConstantBC(0.)
    bc.pInc_diffusive.setConstantBC(0.)
    bc.us_diffusive.setConstantBC(0.)
    bc.vs_diffusive.setConstantBC(0.)
    bc.ws_diffusive.setConstantBC(0.)
    bc.k_diffusive.setConstantBC(0.)
    bc.dissipation_diffusive.setConstantBC(0.)

    bc.clsvof_dirichlet.uOfXT = None

    bc.vof_advective.setConstantBC(0.)
    bc.vos_advective.setConstantBC(0.)

def open( bc, condition ):
    """
    Sets open air bc 
    """
    bc.reset()

    def constant(c):
        def function(x,t): return c
        return function

    fluid      = condition['fluid']
    sdf        = condition['sdf']
    gravity    = condition['gravity']
    waterLevel = condition['waterLevel']
    axis       = condition['gravity_axis']
    epsFact    = condition['epsFact']
    he         = condition['smoothing']
    fac        = 1.0 

    rho_air   = fluid['air']['rho']
    rho_water = fluid['water']['rho']
    mu_air    = fluid['air']['mu']
    mu_water  = fluid['water']['mu']
    phi_air   = fluid['air']['phi']
    phi_water = fluid['water']['phi']

    def p_dirichlet(x,t):
      p_L = default_p.L[axis]*rho_air*gravity[axis]
      phi_L = default_p.L[axis] - waterLevel
      phi = x[axis] - waterLevel
      return p_L - gravity[axis]*(  rho_water*(phi_L - phi)
                                  + (rho_air -rho_water)*(smoothedHeaviside_integral(epsFact*he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact*he,phi)))
    def vof_dirichlet(x,t):
        phi = sdf(x)
        vof = smoothedHeaviside( epsFact*he, phi )
        return vof

    bc.u_dirichlet.uOfXT = None
    bc.u_advective.uOfXT = None
    bc.u_diffusive.uOfXT = constant(0.)

    bc.v_dirichlet.uOfXT = None
    bc.v_advective.uOfXT = None
    bc.v_diffusive.uOfXT = constant(0.)

    bc.w_dirichlet.uOfXT = None
    bc.w_advective.uOfXT = None
    bc.w_diffusive.uOfXT = constant(0.)

    bc.clsvof_dirichlet.uOfXT = None

    bc.vof_dirichlet.uOfXT = vof_dirichlet
    bc.vof_advective.uOfXT = None
    ##.vof_diffusive.uOfXT = does not exist

    bc.p_dirichlet.uOfXT = p_dirichlet
    bc.p_advective.uOfXT = None
    ##.p_diffusive.uOfXT = does not exist

    bc.pInc_dirichlet.uOfXT = constant(0.)
    bc.pInc_advective.uOfXT = None
    bc.pInc_diffusive.uOfXT = None
    
    bc.pInit_dirichlet.uOfXT = constant(0.)
    bc.pInit_advective.uOfXT = None
    bc.pInit_diffusive.uOfXT = None

def noSlip( bc, condition ):
    """
    Sets no slip conditions at the boundary
    """
    bc.reset()

    def constant(c):
        def function(x,t): return c
        return function

    bc.u_dirichlet.uOfXT = constant(0.)
    bc.u_advective.uOfXT = constant(0.)
    bc.u_diffusive.uOfXT = None

    bc.v_dirichlet.uOfXT = constant(0.)
    bc.v_advective.uOfXT = constant(0.)
    bc.v_diffusive.uOfXT = None

    bc.w_dirichlet.uOfXT = constant(0.)
    bc.w_advective.uOfXT = constant(0.)
    bc.w_diffusive.uOfXT = None

    bc.clsvof_dirichlet.uOfXT = None

    bc.vof_dirichlet.uOfXT = None
    bc.vof_advective.setConstantBC(0.)
    ##.vof_diffusive.uOfXT = does not exist

    bc.k_dirichlet.setConstantBC(0.)
    #bc.k_advective.
    bc.k_diffusive.setConstantBC(0.)

    #bc.dissipation_dirichlet.
    #bc.dissipation_advective.
    bc.dissipation_diffusive.setConstantBC(0.)   

    bc.p_dirichlet.uOfXT = None
    bc.p_advective.setConstantBC(0.)
    #bc.p_diffusive.uOfXT = None

    bc.pInc_dirichlet.uOfXT = None
    bc.pInc_advective.setConstantBC(0.)  
    bc.pInc_diffusive.setConstantBC(0.)

    bc.pInit_dirichlet.uOfXT = None
    bc.pInit_advective.setConstantBC(0.)
    bc.pInit_diffusive.setConstantBC(0.)

    bc.vos_advective.setConstantBC(0.)
    bc.us_dirichlet.setConstantBC(0.)
    bc.vs_dirichlet.setConstantBC(0.)
    bc.ws_dirichlet.setConstantBC(0.)  

