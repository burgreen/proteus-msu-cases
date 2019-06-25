#!/usr/bin/env python

# plot-cl-cd.py

# 2019.02.26 initial coding; greg.burgreen@msstate.edu

import sys
import string
import math

#from pylab import *

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#----------------------------------------------------------
def stats( data ):
#----------------------------------------------------------
  avg = 0.0
  for x in data: avg += x
  avg /= len(data)

  std = 0.0
  for x in data: std += (x-avg) * (x-avg)
  std = math.sqrt( std / ( len(data) - 1 ) )

  return (avg, std )

#----------------------------------------------------------
def _parse_res( fp ):
#----------------------------------------------------------
  done = 0
  line = fp.readline()
  while not done:
    if 'Newton' in line:
      done = 1
      s0 = str.split( line, '=' )
      s1 = str.split( s0[1] )
      #s2 = s1[0].rstrip('\r\n,')
      res = float( s1[0] )
      return math.log10(res)
    line = fp.readline()

#----------------------------------------------------------
def s_read_of_output( filename ):
#----------------------------------------------------------
  print('reading ', filename)
  fp = open( filename,'r' )

  line = 1
  residuals = {}
  residuals['i'] = []
  residuals['mom'] = []
  residuals['p'] = []
  residuals['vof'] = []
  residuals['ls'] = []
  residuals['rd'] = []
  residuals['lsc'] = []
  residuals['k'] = []
  residuals['e'] = []
  cnt = 0
  while line:
    line = fp.readline()
    if 'Model' in line:
      i = str.split( line )
      print( i )
      if i[2] == 'Step':
        cnt += 1
        residuals['i'].append( cnt )
      if i[3] == 'rans2p_p':
        residuals['mom'].append( _parse_res( fp ) )
        residuals['p'].append( _parse_res( fp ) )
      if i[3] == 'vof_p':
        residuals['vof'].append( _parse_res( fp ) )
      if i[3] == 'ls_p':
        residuals['ls'].append( _parse_res( fp ) )
      if i[3] == 'redist_p':
        residuals['rd'].append( _parse_res( fp ) )
      if i[3] == 'ls_consrv_p':
        residuals['lsc'].append( _parse_res( fp ) )
      if i[3] == 'turb_k_p':
        residuals['k'].append( _parse_res( fp ) )
      if i[3] == 'turb_e_p':
        residuals['e'].append( _parse_res( fp ) )
  return residuals


  '''
    line = fp.readline().rstrip('\r\n')
    s = str.split( line )

    if len(s) > 0:
      i.append( float(s[0]) )
      cl.append( float(s[cl_idx]) )
      cd.append( float(s[cd_idx]) )
      cnt += 1
      if cnt > monitor_after:
        cl_mon.append(float(s[cl_idx]))
        cd_mon.append(float(s[cd_idx]))

   cl_avg = cl_std = cd_avg = cd_std = 0

   if len(cl_mon) == 0: 
     print('no monitored data: num_its, monitor_after =', cnt, monitor_after )
   else:
     (cl_avg, cl_std) = stats( cl_mon )
     (cd_avg, cd_std) = stats( cd_mon )
 
   return (cl_avg, cl_std, cd_avg, cd_std)
  '''

#----------------------------------------------------------
def _parse_single( fp ):
#----------------------------------------------------------
  done = 0
  line = fp.readline()
  while not done:
    if 'Newton' in line:
      done = 1
      s0 = str.split( line, '=' )
      s1 = str.split( s0[1] )
      #s2 = s1[0].rstrip('\r\n,')
      res = float( s1[0] )
      return math.log10(res)
    line = fp.readline()

#----------------------------------------------------------
def s_read_single( filename, eqn ):
#----------------------------------------------------------
  print('reading ', filename)
  fp = open( filename,'r' )

  line = 1
  residuals = {}
  residuals['i'] = []
  residuals['mom'] = []
  residuals['p'] = []
  residuals['vof'] = []
  residuals['ls'] = []
  residuals['rd'] = []
  residuals['lsc'] = []
  residuals['k'] = []
  residuals['e'] = []
  cnt = 0
  while line:
    line = fp.readline()
    if 'Model' in line:
      i = str.split( line )
      print( i )
      done = 0
      if i[3] == 'rans2p_p' and eqn == 'rans2p':
        while not done:
          residuals['mom'].append( _parse_single( fp ) )
          residuals['p'].append( _parse_single( fp ) )
          line = fp.readline()
          if 'Model' in line: done = 1
          cnt = len(residuals['mom'])
      if i[3] == 'vof_p' and eqn == 'vof':
        while not done:
          residuals['vof'].append( _parse_single( fp ) )
          line = fp.readline()
          if 'Model' in line: done = 1
          cnt = len(residuals['vof'])
      if i[3] == 'ls_p' and eqn == 'ls':
        while not done:
          residuals['ls'].append( _parse_single( fp ) )
          line = fp.readline()
          if 'Model' in line: done = 1
          cnt = len(residuals['ls'])
      if i[3] == 'redist_p' and eqn == 'rd':
        while not done:
          residuals['rd'].append( _parse_single( fp ) )
          line = fp.readline()
          if 'Model' in line: done = 1
          cnt = len(residuals['rd'])
      if i[3] == 'ls_consrv_p' and eqn == 'lsc':
        while not done:
          residuals['lsc'].append( _parse_single( fp ) )
          line = fp.readline()
          if 'Model' in line: done = 1
          cnt = len(residuals['lsc'])
  return residuals, cnt

#-------------------------------------
def main( file ):
#-------------------------------------

  #i = []; cl = [];  cd = []; 

  r = s_read_of_output( file )

  eqn = 'rans2p'
  eqn = 'lsc'
  eqn = 'rd'
  eqn = 'ls'
  eqn = 'vof'
  #r, cnt = s_read_single( file, eqn )
  #for i in range(cnt): r['i'].append(i)
  
  #print('cl avg, stddev', cl_a, cl_sd )
  #print('cd avg, stddev', cd_a, cd_sd )

  fig, ax = plt.subplots()
  if len(r['mom']): ax.plot( np.arange(len(r['mom'])), r['mom'], 'k--', label='mom')
  if len(r['p']):   ax.plot( np.arange(len(r['p'])), r['p'],   'r:',  label='p')
  if len(r['vof']): ax.plot( np.arange(len(r['vof'])), r['vof'], 'k',   label='vof')
  if len(r['ls']):  ax.plot( np.arange(len(r['ls'])), r['ls'],  'r',   label='ls')
  if len(r['rd']):  ax.plot( np.arange(len(r['rd'])), r['rd'],  'g',   label='rd')
  if len(r['lsc']): ax.plot( np.arange(len(r['lsc'])), r['lsc'], 'y',   label='lsc')
  if len(r['k']):   ax.plot( np.arange(len(r['k'])), r['k'],   'g--', label='k')
  if len(r['e']):   ax.plot( np.arange(len(r['e'])), r['e'],   'g:',  label='e')
  #legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')
  legend = ax.legend( loc='upper right', shadow=False )
  
  '''
  pu = plot( residuals['i'], residuals['u'], '-', marker='None', markersize=1, lw=1, c='r' )
  pv = plot( residuals['i'], residuals['v'], '-', marker='None', markersize=1, lw=1, c='k' )
  pw = plot( residuals['i'], residuals['w'], '-', marker='None', markersize=1, lw=1, c='b' )
  pp = plot( residuals['i'], residuals['p'], '-', marker='None', markersize=1, lw=1, c='y' )
  #p2 = plot( i, cd, '-', marker='.', markersize=1, lw=1, c='k' )
  #legend( (pu,pv,pw,pp), ('u','v','w','p'), 'upper center', shadow=False )
  legend( (pu,pv,pw,pp), ('u','v','w','p') )
  '''
  
  plt.xlabel('time steps')
  plt.ylabel('residuals')
  #plt.title('Title')
  plt.grid(True)
  
  #plt.savefig('plot.png')
  #print('- saved plot.png')

  plt.show()

  #matplotlib.axes.Axes.plot
  #matplotlib.pyplot.plot
  #matplotlib.axes.Axes.legend
  #matplotlib.pyplot.legend

#----------------------------------------------------
if __name__ == '__main__':
#----------------------------------------------------
  if len(sys.argv) == 1: print('usage: python plot.py xx.dat [2500]=monitor_after [1]=incomp')
  if len(sys.argv) == 1: exit()

  file = 'undefined'
  if len(sys.argv) > 1: file = str(sys.argv[1])

  #monitor_after = 2500
  #if len(sys.argv) > 2: monitor_after = int(sys.argv[2])

  #incomp = 1
  #if len(sys.argv) > 3: incomp = int(sys.argv[3])

  main( file )
