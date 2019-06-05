#!/bin/bash

if [ -z "$1" ]; then 
  echo usage: $0 \<output_dir\> 
  echo where \<output_dir\> is the preferred name for proteus output
  echo
  echo For example, $0 zoutput 
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

mpirun -n 3 parun main.py -l 2 -v -O 2-petsc.options.asm -D $1

