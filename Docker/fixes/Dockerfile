FROM docker.io/hepstore/rivet:2.6.1
MAINTAINER Holger Schulz <hschulz@fnal.gov>

WORKDIR /code

RUN dnf install -y hdf5-mpich-devel git cmake rsync
RUN dnf clean all -y
RUN git clone https://github.com/diatomic/diy
RUN git clone https://github.com/BlueBrain/HighFive.git

RUN wget https://bitbucket.org/iamholger/pythia8-diy/downloads/pythia8-diy-0.5.1.tar.gz && tar xzf pythia8-diy-0.5.1.tar.gz
RUN wget https://bitbucket.org/iamholger/lheh5/downloads/lheh5-0.2.0.tar.gz && tar xzf lheh5-0.2.0.tar.gz
RUN wget http://home.thep.lu.se/~torbjorn/pythia8/pythia8235.tgz && tar xzf pythia8235.tgz


ENV LD_LIBRARY_PATH="/usr/lib64/mpich/lib:${LD_LIBRARY_PATH}"
ENV PATH="/usr/lib64/mpich/bin:${PATH}"

RUN cd lheh5-0.2.0 && mkdir build && cd build && cmake .. -DHIGHFIVE_DIR=/code/HighFive -DCMAKE_INSTALL_PREFIX=/usr/local && make && make install

COPY Pythia.cc pythia8235/src/Pythia.cc
COPY MergingHooks.cc  pythia8235/src/MergingHooks.cc
COPY TimeShower.cc    pythia8235/src/TimeShower.cc
COPY SpaceShower.cc   pythia8235/src/SpaceShower.cc 
COPY Info.h pythia8235/include/Pythia8/Info.h
COPY Merging.cc pythia8235/src/Merging.cc
RUN cd pythia8235 && ./configure --prefix=/usr/local --enable-shared --with-fastjet3=/usr/local --with-hepmc2=/usr/local && make -j 8 && make install
COPY NLOMerging.xml /usr/local/share/Pythia8/xmldoc/NLOMerging.xml


RUN cd pythia8-diy-0.5.1 && mkdir build && cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local \
    -DDIY_INCLUDE_DIRS=/code/diy/include \
    -DPYTHIA8_DIR=/usr/local \
    -DFASTJET_DIR=/usr/local \
    -DRIVET_DIR=/usr/local \
    -DHEPMC_DIR=/usr/local \
    -DYODA_DIR=/usr/local \
    -DLHEH5_DIR=/usr/local \
    -DHIGHFIVE_DIR=/code/HighFive && make && make install

RUN cd / && rm -rf /code

WORKDIR /examples
COPY ttbar.cmd ttbar.cmd
COPY merging.cmd merging.cmd

COPY PP_JetRates.info /usr/local/share/Rivet/PP_JetRates.info
COPY PP_JetRates.cc PP_JetRates.cc
COPY PP_JetRates.plot /usr/local/share/Rivet/PP_JetRates.plot
RUN rivet-buildplugin /usr/local/lib/RivetPP_JetRates.so PP_JetRates.cc
COPY merging_LO.cmd merging_LO.cmd
COPY concatenate.py concatenate.py
