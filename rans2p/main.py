# main.py

from proteus.default_so import *
import main_param as param
import user_param as user_param

eqns = []
eqns.append('rans2p')

if param.modeVF == 1:
  eqns.append('vof')

if param.modeVF == 2:
  eqns.append('vof')
  eqns.append('ls')
  eqns.append('redist')
  eqns.append('ls_consrv')
    
if param.useRANS > 0:
  eqns.append('turb_k')
  eqns.append('turb_e')

pnList = []
for eqn in eqns: pnList.append( (eqn+'_p',eqn+'_n') )

name = "main" 

systemStepControllerType = Sequential_MinAdaptiveModelStep

if user_param.knob['system']['dt_system_fixed'] > 0.0:
  systemStepControllerType = Sequential_FixedStep
  systemStepExact = False
  dt_system_fixed = user_param.knob['system']['dt_system_fixed']

needEBQ_GLOBAL = False
needEBQ        = False

tnList = [0.0] 
if param.dt_init > 0: 
  tnList += [param.dt_init]
  tnList += [param.dt_init+(i+1)*param.dt_fixed for i in range(0,param.dt_fixed_steps)] 
else:
  tnList += [(i+1)*param.dt_fixed for i in range(0,param.dt_fixed_steps)] 
