#!/bin/bash

source $HOME/singularity/generic_newest.env
for folder in operator conformal propagator
do
  cd $HOME/repository.dev/$folder &&
  mkdir -p build; cd build; rm -rf *;
  cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CXX_FLAGS=-Ofast -DSEPlib_DIR=/opt/SEP/cmake/SEP ..
  make -j 8 install
done
