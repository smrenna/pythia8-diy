#!/bin/bash

DEP=/home/hschulz/scratch-shared/dep

source /products/setup
setup -B mpich v3_2_1a -q +e17:+prof

cmake . \
    -DDIY_INCLUDE_DIRS=${DEP}/diy/include \
    -Ddiy_thread=OFF \
    -DPYTHIA8_DIR=${DEP}/pythia8226 \
    -DRIVET_DIR=${DEP}/local \
    -DHEPMC_DIR=${DEP}/local \
    -DYODA_DIR=${DEP}/local \
    -DCMAKE_CXX_COMPILER=`which g++`
