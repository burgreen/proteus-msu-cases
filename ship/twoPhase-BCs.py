def set_bc_velocityInlet_twoPhase( bc, condition, normal ):
    """
    Sets velocity inlet bc for twoPhase flows
    """
    bc.reset()

    Vmag  = condition['Vmag']
    sdf   = condition['sdf']
    fluid = condition['fluid']
    Vmag_wind = Vmag

    def vel_dirichlet(i):
            def dirichlet(x,t):
                phi = sdf(x)
                H = 1.
                if phi <= 0.: H = 0.
                v = H * Vmag_wind*normal[i] + (1-H) * Vmag*normal[i]
                return v
            return dirichlet

    def vof_dirichlet(x,t):
            phi = sdf(x)
            H = 1.
            if phi <= 0.: H = 0.
            phi_air     = fluid['air']['phi']
            phi_water   = fluid['water']['phi']
            vof = H * phi_air + (1-H) * phi_water
            return vof

    def p_advective(x,t):
            phi = sdf(x)
            H = 1.
            if phi <= 0.: H = 0.
            u = H * Vmag_wind + (1 - H) * Vmag
            return -u

    bc.u_dirichlet.uOfXT = vel_dirichlet(0)
    bc.v_dirichlet.uOfXT = vel_dirichlet(1)
    bc.w_dirichlet.uOfXT = vel_dirichlet(2)

    bc.p_advective.uOfXT = p_advective

    bc.vof_dirichlet.uOfXT = vof_dirichlet

def set_bc_outflow_twoPhase( bc, condition ):
    """
    Sets outflow bc for twoPhase flows
    """
    bc.reset()

    sdf        = condition['sdf']
    fluid      = condition['fluid']
    gravity    = condition['gravity']
    waterLevel = condition['waterLevel']

    def p_dirichlet(x,t):  # assumes a constant fixed fluid height at the outlet
        rho_air     = fluid['air']['rho']
        rho_water   = fluid['water']['rho']
        g_component = gravity[2]
        p_air   = rho_air    * g_component * ( default_p.L[2] - waterLevel )
        p_water = rho_water  * g_component * ( waterLevel - x[2] )
        p_hydrostatic = p_air
        if sdf(x) < 0: p_hydrostatic += p_water
        return -p_hydrostatic

    def u_diffusive(x,t):
        g = p_dirichlet(x,t)
        phi = sdf(x)
        H = 1.
        if phi <= 0.: H = 0.
        mu_air     = fluid['air']['mu']
        mu_water   = fluid['water']['mu']
        return H*(mu_air)*g + (1-H)*(mu_water)*g

    def vof_dirichlet(x,t):
        phi = sdf(x)
        H = 1.
        if phi <= 0.: H = 0.
        phi_air     = fluid['air']['phi']
        phi_water   = fluid['water']['phi']
        return H * phi_air + (1-H) * phi_water

    bc.u_dirichlet.setConstantBC(0.)
    bc.v_dirichlet.setConstantBC(0.)
    bc.w_dirichlet.setConstantBC(0.)

    bc.u_diffusive.uOfXT = u_diffusive
    bc.v_diffusive.setConstantBC(0.)
    bc.w_diffusive.setConstantBC(0.)

    bc.p_dirichlet.uOfXT   = p_dirichlet

    bc.vof_dirichlet.uOfXT = vof_dirichlet
