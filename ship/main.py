# main.py

from proteus.default_so import *
import main_param as param

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

needEBQ_GLOBAL = False
needEBQ        = False

tnList = [0.0] 
if param.dt_init > 0: tnList += [param.dt_init]
tnList += [(i+1)*param.dt_fixed for i in range(0,param.dt_fixed_steps)] 
