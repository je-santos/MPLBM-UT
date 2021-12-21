import os
import subprocess
import sys
sys.path.append('../../python_utils/')  # It would be nice to make this a proper package...

import numpy as np
import matplotlib.pyplot as plt

from create_geom_for_palabos import *
from create_palabos_input_file import *
from parse_input_file import *
from create_geom_for_rel_perm import *
from parse_palabos_output import *
from create_plots import *


def download_geometry(filename, url):

    download_command = f'wget {url} -O {filename}'
    try:
        subprocess.run(download_command.split(' '))
    except FileNotFoundError:
        raise InterruptedError(f'wget was not found. Try running this in the terminal: \n {download_command}')
    return


def run_2_phase_sim(inputs):
    # Steps
    # 1) create geom for palabos
    # 2) create palabos input file
    # 3) run 2-phase sim

    if inputs['simulation type'] == '1-phase':
        raise KeyError('Simulation type set to 1-phase...please change to 2-phase.')
    sim_directory = inputs['input output']['simulation directory']

    # 1) Create Palabos geometry
    print('Creating efficient geometry for Palabos...')
    create_geom_for_palabos(inputs)

    # 2) Create simulation input file
    print('Creating input file...')
    create_palabos_input_file(inputs)

    # 3) Run 2-phase simulation
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


def run_rel_perm_sim(inputs):

    # Rel Perm Steps
    # 1) create 1-phase palabos input file for relperms
    # 2) create geoms for rel perm
    # 3) run 1-phase sim to get rel perms

    if inputs['simulation type'] == '1-phase':
        raise KeyError('Simulation type set to 1-phase...please change to 2-phase.')
    sim_directory = inputs['input output']['simulation directory']

    # 1) Create geoms for rel perm
    print('Creating rel perm geometries...')
    user_geom_name = inputs['domain']['geom name']
    inputs = create_geom_for_rel_perm(inputs)

    # 2) Create simulation input file
    print('Creating input file...')
    inputs['simulation type'] = 'rel perm'
    create_palabos_input_file(inputs)
    inputs['domain']['geom name'] = user_geom_name

    # 3) Write rel perm bash file
    print('Running rel perm simulation...')
    num_procs = inputs['simulation']['num procs']
    input_dir = inputs['input output']['input folder']
    simulation_command = f"mpirun -np {num_procs} ../../src/1-phase_LBM/permeability {input_dir}relperm_input.xml"
    file = open(f'{sim_directory}/{input_dir}run_relperm_sim.sh', 'w')
    file.write(f'{simulation_command}')
    file.close()

    simulation_command_subproc = f'bash {sim_directory}/{input_dir}run_relperm_sim.sh'
    make_4relperm_folder = 'mkdir tmp/4relperm'
    subprocess.run(make_4relperm_folder.split(' '))
    subprocess.run(simulation_command_subproc.split(' '))

    return


def process_and_plot_results(inputs):

    # Process data
    create_pressure_data_file(inputs)
    create_relperm_data_file(inputs)

    # Load and plot data
    sim_dir = inputs['input output']['simulation directory'] + '/'
    output_dir = inputs['input output']['output folder']
    Sw = np.loadtxt(f'{sim_dir + output_dir}data_Sw.txt')
    Pc = np.loadtxt(f'{sim_dir + output_dir}data_Pc.txt')
    krw = np.loadtxt(f'{sim_dir + output_dir}data_krw.txt')
    krnw = np.loadtxt(f'{sim_dir + output_dir}data_krnw.txt')

    plt.figure(figsize=[10,8])
    plot_capillary_pressure_data(Sw, Pc)
    plt.savefig(sim_dir + 'pc_curve.png', dpi=300)

    plt.figure(figsize=[10,8])
    plot_rel_perm_data(Sw, krw, krnw)
    plt.savefig(sim_dir + 'relperm_curve.png', dpi=300)

    plt.figure(figsize=[6, 10])
    plot_pc_and_rel_perm(Sw, Pc, krw, krnw)
    plt.savefig(sim_dir + 'pc_and_relperm_curve.png', dpi=300)

    plt.show()

    return


drp_url = 'www.digitalrocksportal.org/projects/65/images/71108/download/'
file_name = 'input/rg_theta60_phi10.raw'
download_geometry(file_name, drp_url)

input_file = 'input.yml'
inputs = parse_input_file(input_file)  # Parse inputs
run_2_phase_sim(inputs)  # Run 2 phase sim
run_rel_perm_sim(inputs)  # Run rel perm
process_and_plot_results(inputs)  # Plot results

