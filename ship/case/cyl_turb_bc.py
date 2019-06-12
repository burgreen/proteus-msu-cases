# cyl_turb_bc.py

import sys
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus import default_p

def noSlip( bc, condition ):
    """
    Sets no slip conditions at the boundary
    """
    bc.reset()
    bc.BC_type = 'noSlip'

    bc.u_dirichlet.setConstantBC(0.)
    bc.u_advective.uOfXT = None
    #bc.u_diffusive.setConstantBC(0.) # if enabled, no BL forms

    bc.v_dirichlet.setConstantBC(0.)
    bc.v_advective.uOfXT = None
    #bc.v_diffusive.setConstantBC(0.) # if enabled, no BL forms

    bc.w_dirichlet.setConstantBC(0.)
    bc.w_advective.uOfXT = None
    #bc.w_diffusive.setConstantBC(0.) # if enabled, no BL forms

    bc.vof_dirichlet.uOfXT = None
    bc.vof_advective.setConstantBC(0.)
    ##.vof_diffusive.uOfXT = does not exist

    bc.k_dirichlet.setConstantBC(0.)
    #bc.k_advective.
    #bc.k_diffusive.setConstantBC(0.)

    bc.dissipation_dirichlet.setConstantBC(0.0375)  # works for k-omega model
    #bc.dissipation_advective.
    #bc.dissipation_diffusive.setConstantBC(0.)   

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

