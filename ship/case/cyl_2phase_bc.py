# cyl_2phase_bc.py

import sys
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus import default_p

def velocityInlet_rans2p( bc, condition, normal ):
    """
    Sets velocity inlet bc for twoPhase flows
    """
    bc.reset()

    Vmag = 0.
    if 'Vmag' in condition: Vmag = condition['Vmag']
    Vmag_air   = Vmag*0
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
            H = 1.
            if phi <= 0.: H = 0.
            vof = H * phi_air + (1-H) * phi_water
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

    bc.k_dirichlet.uOfXT = constant(turb_k)
    #bc.k_advective.
    #bc.k_diffusive.

    bc.dissipation_dirichlet.uOfXT = constant(turb_e)
    #bc.dissipation_advective.
    #bc.dissipation_diffusive.

    bc.vof_dirichlet.uOfXT = vof_dirichlet
    bc.vof_advective.uOfXT = constant(0.)
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
    smoothing  = 0.05 
    fac        = 0.1 

    rho_air   = fluid['air']['rho']
    rho_water = fluid['water']['rho']
    mu_air    = fluid['air']['mu']
    mu_water  = fluid['water']['mu']
    phi_air   = fluid['air']['phi']
    phi_water = fluid['water']['phi']

    def p_dirichlet(x,t):  # assumes a constant fixed fluid height at the outlet
        g_component = gravity[1]
        p_air   = rho_air    * g_component * ( default_p.L[1] - waterLevel )
        p_water = rho_water  * g_component * ( waterLevel - x[1] )
        p_hydrostatic = p_air
        #if sdf(x) < 0: p_hydrostatic = p_water
        if x[1] < waterLevel: p_hydrostatic = p_water
        return -p_hydrostatic

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
        if phi >= smoothing: H = 1.
        elif -smoothing < phi < smoothing: H = smoothedHeaviside(smoothing, phi)
        elif phi <= -smoothing: H = 0.
        #H = 1.
        #if phi <= 0.: H = 0.
        return H * phi_air + (1-H) * phi_water

    def constant(c):
        def function(x,t): return c
        return function

    bc.u_dirichlet.uOfXT = constant(0.)
    bc.u_advective.uOfXT = None
    bc.u_diffusive.uOfXT = u_diffusive # some level is needed for stability

    bc.v_dirichlet.uOfXT = constant(0.)
    bc.v_advective.uOfXT = None
    bc.v_diffusive.uOfXT = constant(0.)

    bc.w_dirichlet.uOfXT = constant(0.)
    bc.w_advective.uOfXT = None
    bc.w_diffusive.uOfXT = constant(0.)

    bc.k_dirichlet.uOfXT = constant(0.)
    bc.k_advective.uOfXT = None
    bc.k_diffusive.uOfXT = constant(0.)

    bc.dissipation_dirichlet.uOfXT = constant(0.)
    bc.dissipation_advective.uOfXT = None
    bc.dissipation_diffusive.uOfXT = constant(0.)

    bc.vof_dirichlet.uOfXT = vof_dirichlet
    bc.vof_advective.uOfXT = None
    ##.vof_diffusive.uOfXT = does not exist

    bc.p_dirichlet.uOfXT = p_dirichlet
    bc.p_advective.uOfXT = None
    ##.p_diffusive.uOfXT = does not exist

