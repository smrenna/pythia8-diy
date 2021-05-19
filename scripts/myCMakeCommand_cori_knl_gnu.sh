#!/bin/bash

rm -rf CMakeFiles CMakeCache.txt 
DEP=/global/cscratch1/sd/hschulz/HEP-gnu-knl
DIY=/global/cscratch1/sd/hschulz/diy

module swap PrgEnv-intel PrgEnv-gnu

export CRAYPE_LINK_TYPE=dynamic
cmake . \
    -DDIY_INCLUDE_DIRS=${DIY}/include \
    -Ddiy_thread=OFF \
    -DPYTHIA8_DIR=${DEP}/pythia8226 \
    -DRIVET_DIR=${DEP}/local \
    -DHEPMC_DIR=${DEP}/local \
    -DYODA_DIR=${DEP}/local \
    -DCMAKE_CXX_COMPILER=`which CC`

make clean && make 
cp pythia8-diy pythia8-diy-gnu-knl
