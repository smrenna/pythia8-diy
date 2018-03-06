#!/bin/bash

cmake . \
    -DDIY_INCLUDE_DIRS=/home/hschulz/scidac/diy/include \
    -Ddiy_thread=OFF -DPYTHIA8_DIR=~/src/pythia8230/ \
    -DRIVET_DIR=/home/hschulz/scidac/parallelhep/dep/local \
    -DHEPMC_DIR=/home/hschulz/scidac/parallelhep/dep/local \
    -DYODA_DIR=/home/hschulz/scidac/parallelhep/dep/local \
    -DPYTHIA8_DIR=~/src/pythia8230
