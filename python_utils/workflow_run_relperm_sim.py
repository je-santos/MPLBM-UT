import subprocess
import os
import sys
from create_geom_for_rel_perm import *
from create_palabos_input_file import *
from parse_input_file import *

# Steps
# 1) parse 2-phase input file
# 2) create 1-phase palabos input file for relperms
# 3) create geoms for rel perm
# 4) run 1-phase sim to get rel perms

input_file = sys.argv[1]

# 1) Process input file
inputs = parse_input_file(input_file)

if inputs['simulation type'] == '1-phase':
    raise KeyError('Simulation type set to 1-phase...please change to 2-phase.')
sim_directory = inputs['input output']['simulation directory']

# 2) Create geoms for rel perm
print('Creating rel perm geometries...')
user_geom_name = inputs['domain']['geom name']
inputs = create_geom_for_rel_perm(inputs)

# 3) Create simulation input file
print('Creating input file...')
inputs['simulation type'] = 'rel perm'
create_palabos_input_file(inputs)
inputs['domain']['geom name'] = user_geom_name

# 4) Write rel perm bash file
print('Running rel perm simulation...')
num_procs = inputs['simulation']['num procs']
input_dir = inputs['input output']['input folder']
simulation_command = f"mpirun -np {num_procs} ../../src/1-phase_LBM/permeability {input_dir}relperm_input.xml"
file = open(f'{sim_directory}/{input_dir}run_relperm_sim.sh', 'w')
file.write(f'{simulation_command}')
file.close()

