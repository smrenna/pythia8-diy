#!/bin/bash

export CRAYPE_LINK_TYPE=dynamic

DIY=/global/u1/x/xju/code/diy
DEP=/global/homes/x/xju/hep_packages/Extra
HIGHFIVE=/global/homes/x/xju/code/HighFive

cmake .. \
	-DCMAKE_C_COMPILER=`which cc` \
	-DCMAKE_CXX_COMPILER=`which CC` \
	-DDIY_INCLUDE_DIRS=${DIY}/include \
	-Ddiy_thread=OFF \
	-DRIVET_DIR=/global/homes/x/xju/code/rivet-dev/local \
	-DHEPMC_DIR=${DEP}/HepMC-2.06.09 \
	-DYODA_DIR=${DEP}/YODA-1.7.0 \
	-DFASTJETS_DIR=${DEP}/fastjet-3.3.0 \
	-DPYTHIA8_DIR=${DEP}/pythia8235 \
	-DHIGHFIVE_DIR=${HIGHFIVE}
