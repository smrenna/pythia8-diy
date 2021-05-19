mkdir -p code

mkdir -p code/diy 
cd code/diy 
wget https://github.com/diatomic/diy/archive/3.5.0.tar.gz 
tar x --strip-components 1 -f 3.5.0.tar.gz 
rm 3.5.0.tar.gz
cd -

mkdir -p code/HighFive 
cd code/HighFive 
wget https://github.com/BlueBrain/HighFive/archive/v2.0.tar.gz 
tar x --strip-components 1 -f v2.0.tar.gz 
rm v2.0.tar.gz
cd -

BIGDIR=/wclustre/sherpa/mrenna/local
VERSION=pythia83-vincia-merging-bugfix.tar.gz 
VERSION=pythia83-pythia8305.tar.gz
mkdir -p code/pythia83
cd code/pythia83
wget https://github.com/smrenna/pythia83/raw/main/$VERSION
tar xzf $VERSION --strip-components 1
./configure --with-fastjet3=$BIGDIR --with-hepmc2=$BIGDIR \
   --with-highfive=/wclustre/sherpa/mrenna/pythia8-diy/code/HighFive --with-hdf5=$HDF5_INC
cp ../../patch/LH*.h include/PythiaPlugins/
make -j 11
rm -rf $VERSION
cd -

rm -rf build
mkdir build && cd build 
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD \
-DCMAKE_CXX_COMPILER=`which g++` \
-DDIY_INCLUDE_DIRS=/wclustre/sherpa/mrenna/pythia8-diy/code/diy/include \
-DPYTHIA8_DIR=/wclustre/sherpa/mrenna/pythia8-diy/code/pythia83 \
-DFASTJET_DIR=$BIGDIR \
-DRIVET_DIR=$BIGDIR \
-DHEPMC_DIR=$BIGDIR \
-DYODA_DIR=$BIGDIR \
-DHIGHFIVE_DIR=/wclustre/sherpa/mrenna/pythia8-diy/code/HighFive && make && make install

