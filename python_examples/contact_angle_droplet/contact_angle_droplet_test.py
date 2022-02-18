import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import sys
sys.path.append('../../python_utils/')  # It would be nice to make this a proper package...
from create_geom_for_palabos import *
from create_palabos_input_file import *
from parse_input_file import *


def create_flat_surface(inputs):

    domain = np.zeros((80, 80, 80), dtype=np.uint8)  # empty space
    domain[0:5,:,:] = 1  # surface
    # domain = np.transpose(domain, [1,0,2])
    domain.tofile(f"{inputs['input output']['input folder']}/flat_surface.raw")
    # print(domain.shape)

    print('Creating efficient geometry for Palabos...')
    create_geom_for_palabos(inputs)

    # import vedo as vd
    # vp = vd.Plotter(axes=3, bg='w', bg2='w', size=(1200, 900), offscreen=False)
    # domain = vd.Volume(domain).isosurface(0.5)
    # vp += domain.c('gray')
    # vp.show()

    return


def run_2_phase_sim(inputs):
    # Steps
    # 1) create palabos input file
    # 2) run 2-phase sim

    if inputs['simulation type'] == '1-phase':
        raise KeyError('Simulation type set to 1-phase...please change to 2-phase.')
    sim_directory = inputs['input output']['simulation directory']

    # 1) Create simulation input file
    print('Creating input file...')
    create_palabos_input_file(inputs)

    # 2) Run 2-phase simulation
    print('Running 2-phase simulation...')
    num_procs = inputs['simulation']['num procs']
    input_dir = inputs['input output']['input folder']
    simulation_command = f"mpirun -np {num_procs} ../../src/2-phase_LBM/ShanChen {input_dir}2_phase_sim_input.xml"
    file = open(f'{sim_directory}/{input_dir}run_shanchen_sim.sh', 'w')
    file.write(f'{simulation_command}')
    file.close()

    simulation_command_subproc = f'bash {sim_directory}/{input_dir}run_shanchen_sim.sh'
    subprocess.run(simulation_command_subproc.split(' '))

    return


input_file = 'input.yml'
inputs = parse_input_file(input_file)  # Parse inputs
inputs['input output']['simulation directory'] = os.getcwd()  # Store current working directory
create_flat_surface(inputs)

G_ads_f1_s1 = np.array([-0.4, 0, 0.4])  # G_ads_f1_s1 values corresponding to 156 deg, 90 deg, and 24 deg contact angles
tmp_folders = ['tmp_m0.4/', 'tmp_0/', 'tmp_p0.4/']
for i in range(len(G_ads_f1_s1)):
    inputs['input output']['output folder'] = tmp_folders[i]
    inputs['simulation']['fluid data']['G_ads_f1_s1'] = G_ads_f1_s1[i]
    print(f"\n\nRunning droplet test with G_ads_f1_s1 = {G_ads_f1_s1[i]}...")
    run_2_phase_sim(inputs)  # Run 2 phase sim


