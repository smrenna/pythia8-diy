#./bin/env bash
# built with gnu8, cmake, mpich, hdf5

BIGDIR=/wclustre/sherpa/mrenna/local
INSTALL_DIR=/wclustre/sherpa/mrenna
rm -rf ./build
mkdir build && cd build 
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
-DCMAKE_CXX_COMPILER=`which g++` \
-DDIY_INCLUDE_DIRS=/wclustre/sherpa/mrenna/pythia8-diy/code/diy/include \
-DPYTHIA8_DIR=/wclustre/sherpa/mrenna/pythia8-diy/code/pythia83 \
-DFASTJET_DIR=$BIGDIR \
-DRIVET_DIR=$BIGDIR \
-DHEPMC_DIR=$BIGDIR \
-DYODA_DIR=$BIGDIR \
-DHDF5_INC=$HDF5_INC \
-DHIGHFIVE_DIR=/wclustre/sherpa/mrenna/pythia8-diy/code/HighFive && make && make install
