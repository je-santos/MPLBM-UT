!#bin/bash

cd src
unzip palabos-v2.zip
cd 2-phase_LBM
make
cd ..
cd 1-phase_LBM
make
