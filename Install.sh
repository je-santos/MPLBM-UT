#!/bin/bash

cd src
unzip palabos-v2.2.1.zip

cd 2-phase_LBM/build
cmake ..
make -j 2

cd ../../1-phase_LBM/build
cmake ..
make -j 2

cd ../../
python -m pip install --upgrade pip
python -m pip install --upgrade setuptools
pip install python/
