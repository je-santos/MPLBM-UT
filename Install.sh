#!/bin/bash

cd src
unzip palabos-v2.2.1.zip

cd 2-phase_LBM/build
cmake ..
make

cd ../../1-phase_LBM/build
cmake ..
make
