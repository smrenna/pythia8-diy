FROM docker.io/hepstore/rivet:2.6.1
MAINTAINER Holger Schulz <hschulz@fnal.gov>

WORKDIR /code

RUN dnf install -y hdf5-mpich-devel git cmake rsync
RUN dnf clean all -y
RUN git clone https://github.com/diatomic/diy
RUN git clone https://github.com/BlueBrain/HighFive.git

RUN wget https://bitbucket.org/iamholger/pythia8-diy/downloads/pythia8-diy-0.1.2.tar.gz && tar xzf pythia8-diy-0.1.2.tar.gz
RUN wget https://bitbucket.org/iamholger/lheh5/downloads/lheh5-0.1.1.tar.gz && tar xzf lheh5-0.1.1.tar.gz
RUN wget http://home.thep.lu.se/~torbjorn/pythia8/pythia8235.tgz && tar xzf pythia8235.tgz


ENV LD_LIBRARY_PATH="/usr/lib64/mpich/lib:${LD_LIBRARY_PATH}"
ENV PATH="/usr/lib64/mpich/bin:${PATH}"

RUN cd lheh5-0.1.1 && mkdir build && cd build && cmake .. -DHIGHFIVE_DIR=/code/HighFive -DCMAKE_INSTALL_PREFIX=/usr/local && make && make install

COPY Pythia.cc pythia8235/src/Pythia.cc 
RUN cd pythia8235 && ./configure --prefix=/usr/local --enable-shared --with-hepmc2=/usr/local && make -j 8 && make install

RUN cd pythia8-diy-0.1.2 && mkdir build && cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local \
    -DDIY_INCLUDE_DIRS=/code/diy/include \
    -DPYTHIA8_DIR=/usr/local \
    -DRIVET_DIR=/usr/local \
    -DHEPMC_DIR=/usr/local \
    -DYODA_DIR=/usr/local \
    -DLHEH5_DIR=/usr/local \
    -DHIGHFIVE_DIR=/code/HighFive && make && make install

RUN cd / && rm -rf /code

WORKDIR /examples
COPY ttbar.cmd ttbar.cmd

