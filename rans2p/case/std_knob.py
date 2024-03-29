# std_knob.py

knob = {}

knob['ns'] = {}
knob['ns']['sc_fac']  = 0.5
knob['ns']['sc_lag']  = True
knob['ns']['sge_lag'] = True
knob['ns']['nl_atol'] = 1.e-6
knob['ns']['nl_its']  = 50
knob['ns']['nl_test'] = 'rits'

knob['vof'] = {}
knob['vof']['sc_fac']  = 0.5
knob['vof']['sc_lag']  = True
knob['vof']['sc_uref'] = 1.0
knob['vof']['sc_beta'] = 1.5
knob['vof']['nl_atol'] = 1.e-6
knob['vof']['nl_its']  = 50
knob['vof']['nl_test'] = 'r'

knob['ls'] = {}
knob['ls']['sc_fac']  = 0.5
knob['ls']['sc_lag']  = True
knob['ls']['sc_uref'] = 1.0
knob['ls']['sc_beta'] = 1.5
knob['ls']['nl_atol'] = 1.e-6
knob['ls']['nl_its']  = 50
knob['ls']['nl_test'] = 'r'

knob['rd'] = {}
knob['rd']['sc_fac']  = 0.5
knob['rd']['sc_lag']  = False
knob['rd']['newton']  = True
knob['rd']['nl_atol'] = 1.e-6
knob['rd']['nl_its']  = 25
knob['rd']['nl_test'] = 'rits'

knob['mcorr'] = {}
knob['mcorr']['nl_atol'] = 1.e-6
knob['mcorr']['nl_its']  = 25
knob['mcorr']['nl_test'] = 'r'

knob['turb_k'] = {}
knob['turb_k']['sc_fac']  = 0.5
knob['turb_k']['sc_lag']  = True
knob['turb_k']['sc_uref'] = 1.0
knob['turb_k']['sc_beta'] = 1.0
knob['turb_k']['nl_atol'] = 1.e-6
knob['turb_k']['nl_its']  = 25
knob['turb_k']['nl_test'] = 'rits'

knob['turb_e'] = {}
knob['turb_e']['sc_fac']  = 0.5
knob['turb_e']['sc_lag']  = True
knob['turb_e']['sc_uref'] = 1.0
knob['turb_e']['sc_beta'] = 1.0
knob['turb_e']['nl_atol'] = 1.e-6
knob['turb_e']['nl_its']  = 25
knob['turb_e']['nl_test'] = 'rits'

knob['epsFact'] = {}
knob['epsFact']['std'] = 1.5
# epsFact.rd controls how far away from the interface we "freeze" the result of the 
# level set transport step.  It is known that redistancing can introduce 
# noise into the interface. Values of 0.75 and 1.5 are known to work.
knob['epsFact']['rd']  = 0.33
knob['epsFact']['consrv_diffusion'] = 10.0

knob['system'] = {}
# global var in default_n.py; controls overall timestep size; default = 0.9
knob['system']['runCFL'] = 0.9 
# if positive, imposes Sequential_FixedStep time steps
knob['system']['dt_system_fixed'] = -0.001

import copy

def set_defaults():
    return copy.deepcopy( knob )
