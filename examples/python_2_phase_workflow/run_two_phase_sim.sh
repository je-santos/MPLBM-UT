#!/bin/bash

# Download geometry from Digital Rocks Portal
wget www.digitalrocksportal.org/projects/65/images/71108/download/ -O input/rg_theta60_phi10.raw

# Run single phase python workflow here
python ../../python_utils/workflow_run_2_phase_sim.py input.yml

# Run two-phase simulation
bash input/run_shanchen_sim.sh > sim_log.txt
