#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then 
  echo usage: $0 \<num_processors\> \<case_name\> \<output_dir\>   
  echo where \<case_name\> is a case_name.py file existing in the case dir
  echo where \<output_dir\> is the preferred name for proteus output
  echo
  echo For example, $0 4 cyl_1phase z-cyl-1phase 
  exit
fi

if ! which parun > /dev/null 2>&1; then
  echo
  echo parun is not found. You need to:
  echo
  echo '$ source <proteus_dir>/0-setup-proteus.sh'
  echo
  echo 'where <proteus_dir> is a location of a valid Proteus installation'
  exit
fi

echo from case import $2 as user_param > user_param.py

if [ -z "$4" ]; then 
  mpirun -n $1 parun main.py -l 2 -v -O 4-petsc.options -D $3
else
  if [ "$4" == "restart" ]; then 
    mpirun -n $1 parun main.py -l 2 -v -O 4-petsc.options -D $3 -H
  else
    echo Option is not supported: $4
  fi
fi

