import os
import subprocess
import sys
sys.path.append('../../python_utils/')  # It would be nice to make this a proper package...
import numpy as np
import matplotlib.pyplot as plt
import re
import glob
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
        raise InterruptedError(f'wget was not found. Please make sure it is installed.')
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
    output_dir = inputs['input output']['output folder']
    simulation_command = f"mpirun -np {num_procs} ../../src/1-phase_LBM/permeability {input_dir}relperm_input.xml"
    file = open(f'{sim_directory}/{input_dir}run_relperm_sim.sh', 'w')
    file.write(f'{simulation_command}')
    file.close()

    simulation_command_subproc = f'bash {sim_directory}/{input_dir}run_relperm_sim.sh'
    make_4relperm_folder = f'mkdir {sim_directory}/{output_dir}/4relperm'
    subprocess.run(make_4relperm_folder.split(' '))
    subprocess.run(simulation_command_subproc.split(' '))

    return inputs


def process_and_plot_results(inputs):

    # Process data
    # create_pressure_data_file(inputs)
    create_relperm_data_file(inputs)

    # Load and plot data
    sim_dir = inputs['input output']['simulation directory'] + '/'
    output_dir = inputs['input output']['output folder']
    Sw = np.loadtxt(f'{sim_dir + output_dir}data_Sw.txt')
    # Pc = np.loadtxt(f'{sim_dir + output_dir}data_Pc.txt')
    krw = np.loadtxt(f'{sim_dir + output_dir}data_krw.txt')
    krnw = np.loadtxt(f'{sim_dir + output_dir}data_krnw.txt')

    # plt.figure(figsize=[10,8])
    # plot_capillary_pressure_data(Sw, Pc)
    # plt.savefig(sim_dir + 'pc_curve.png', dpi=300)

    plt.figure(figsize=[10,8])
    plot_rel_perm_data(Sw, krw, krnw)
    # plt.yscale('log')
    plt.savefig(sim_dir + 'relperm_curve.png', dpi=300)

    # plt.figure(figsize=[6, 10])
    # plot_pc_and_rel_perm(Sw, Pc, krw, krnw)
    # plt.savefig(sim_dir + 'pc_and_relperm_curve.png', dpi=300)

    # plt.show()

    return


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


# Some notes
# If you want to run these individually/multiple at a time (ie multiple jobs on a supercomputer),
# then you'll want to use this for setup but don't run anything (comment out any subprocess.run()).
# From there, you can write a function that will run the bash files you want in each saturation
# folder as a job (follow your system's instructions on submitting simulation runs/jobs).

# Download Geometry #
#####################
drp_url = 'https://www.digitalrocksportal.org/projects/47/images/73328/download/'
file_name = 'input/finney_pack.raw'
download_geometry(file_name, drp_url)

# Rescale geometry #
####################
print("Rescaling geometry...")
geom = np.fromfile(file_name, dtype='uint8').reshape([500, 500, 500])
scaled_geom = scale_geometry(geom, 0.5, 'uint8')
scaled_geom.tofile('input/finney_pack_half_scale.raw')

# Parse Inputs #
################
input_file = 'input.yml'
inputs = parse_input_file(input_file)  # Parse inputs
inputs['input output']['simulation directory'] = os.getcwd()  # Store current working directory


# Saturation folders #
######################
# Note: this is just going off of domain size since we do not have a morphological drainage algorithm in place.
# Another option that we may implement is using a morphological drainage from PoreSpy to get initial conditions.
# Snw_str = np.array(['10', '30', '50', '70', '90', '95'])
Snw_str = np.array(['95'])
Snw_num = Snw_str.astype('float') / 100
sim_dir = inputs['input output']["simulation directory"]

# Run 2 phase sims #
####################
for i in range(len(Snw_str)):

    # Create output folder
    output_dir = f'{sim_dir}/tmp_{Snw_str[i]}'
    output_dir_exists = os.path.isdir(output_dir)
    if output_dir_exists == False:
        os.makedirs(output_dir)

    # Change inputs
    nx = inputs["domain"]["domain size"]["nx"]
    slices = inputs["domain"]["inlet and outlet layers"]

    inputs["input output"]["output folder"] = f"tmp_{Snw_str[i]}/"
    inputs["simulation"]["fluid 1 init"]["x1"] = 0
    inputs["simulation"]["fluid 1 init"]["x2"] = int(nx * (Snw_num[i]))
    inputs["simulation"]["fluid 2 init"]["x1"] = int(nx * (Snw_num[i])) + 1
    inputs["simulation"]["fluid 2 init"]["x2"] = nx + slices*2

    # run_2_phase_sim(inputs)  # Run 2 phase sim

Snw_str = np.array(['10', '30', '50', '70', '90', '95'])

# Organize outputs #
####################
# Gather final steady state 2 phase configs for rel perm
# Create output folder
rel_perm_dir = f'{sim_dir}/relperm'
rel_perm_dir_exists = os.path.isdir(rel_perm_dir)
if rel_perm_dir_exists == False:
    os.makedirs(rel_perm_dir)

inputs["input output"]["output folder"] = "relperm/"  # Change output folder to match

# Copy last configs to a new output folder for a smooth rel perm run
# for each saturation folder
#   - glob all rho_f1*.dat into a list
#   - natural sort the list
#   - subprocess.run() copy config file to rel_perm_dir
sim_dir = inputs['input output']['simulation directory']
for i in range(len(Snw_str)):
    output_dir = f"{sim_dir}/tmp_{Snw_str[i]}/"
    # Glob all .dat fluid files from simulation output folder
    rho_files_regex = fr'{output_dir}rho_f1*.dat'
    rho_files = glob.glob(rho_files_regex)
    rho_files = natural_sort(rho_files)
    last_rho_config = rho_files[-1]  # Add index so we're not overwriting
    cp_command = f'cp {last_rho_config} {rel_perm_dir}/rho_f1_{i+1}.dat'
    subprocess.run(cp_command.split(' '))

# Run relperm and process results #
###################################
inputs = run_rel_perm_sim(inputs)  # Run rel perm
process_and_plot_results(inputs)  # Plot results

