import subprocess
import os
from create_geom_for_palabos import *
from create_palabos_input_file import *
from parse_input_file import *

# Steps
# 1) parse input file
# 2) create geom for palabos
# 3) create palabos input file
# 4) run 1-phase sim
# Please see python_1_phase_workflow in examples directory on how to use this

input_file = sys.argv[1]

# 1) Process input file
inputs = parse_input_file(input_file)

if inputs['simulation type'] == '2-phase':
    raise KeyError('Simulation type set to 2-phase...please change to 1-phase.')
sim_directory = inputs['input output']['simulation directory']

# 2) Create Palabos geometry
print('Creating efficient geometry for Palabos...')
create_geom_for_palabos(inputs)

# 3) Create simulation input file
print('Creating input file...')
create_palabos_input_file(inputs)

# 4) Run 1-phase simulation
print('Running 1-phase simulation...')
num_procs = inputs['simulation']['num procs']
input_dir = inputs['input output']['input folder']
simulation_command = f"mpirun -np {num_procs} ../../src/1-phase_LBM/permeability {input_dir}1_phase_sim_input.xml"
file = open(f'{sim_directory}/{input_dir}run_single_phase_sim.sh', 'w')
file.write(f'{simulation_command}')
file.close()
# simulation_command_subproc = simulation_command.split(' ')
# subprocess.run(simulation_command_subproc)
