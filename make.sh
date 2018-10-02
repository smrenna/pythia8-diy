#!/bin/bash

export CRAYPE_LINK_TYPE=dynamic

cmake .. \
	-DDIY_INCLUDE_DIRS=/global/u1/x/xju/code/diy/include \
	-Ddiy_thread=OFF \
	-DRIVET_DIR=/global/homes/x/xju/code/rivet-dev/local \
	-DHEPMC_DIR=/global/homes/x/xju/hep_packages/Extra/HepMC-2.06.09 \
	-DYODA_DIR=/global/homes/x/xju/hep_packages/Extra/YODA-1.7.0 \
	-DFASTJETS_DIR=/global/homes/x/xju/hep_packages/Extra/fastjet-3.3.0 \
	-DPYTHIA8_DIR=/global/homes/x/xju/hep_packages/Extra/pythia8235 \
	-DCMAKE_CXX_COMPILER=`which CC` \
	-DHDF5_ROOT=/global/homes/x/xju/hep_packages/Extra/HDF5-1.10.1
	#-DHDF5_DIR=/global/homes/x/xju/hep_packages/HDF5/hdf5-1.10.1/build \
