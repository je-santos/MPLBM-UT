import os
import subprocess
import sys
import skimage.transform as skit
import numpy as np
import matplotlib.pyplot as plt
import mplbm_utils as mplbm


def download_geometry(filename, url):

    download_command = f'wget {url} -O {filename}'
    try:
        subprocess.run(download_command.split(' '))
    except FileNotFoundError:
        raise InterruptedError(f'wget was not found. Please make sure it is installed on your system.')
    return


def create_micromodel(geom, nx, ny, nz, slice_index, x_offset, y_offset):

    # Get slice of interest
    slice = geom[slice_index, y_offset:ny+y_offset, x_offset:nx+x_offset]

    num_slices = nz
    micromodel = np.repeat(slice[np.newaxis, :, :], num_slices, axis=0)

    print(f"Micromodel size = {micromodel.shape}")

    # Check micromodel
    plt.figure()
    plt.imshow(micromodel[0,:,:])
    plt.gca().invert_yaxis()
    plt.show()

    return micromodel


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
    mplbm.create_geom_for_palabos(inputs)

    # 2) Create simulation input file
    print('Creating input file...')
    mplbm.create_palabos_input_file(inputs)

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
    inputs = mplbm.create_geom_for_rel_perm(inputs)

    # 2) Create simulation input file
    print('Creating input file...')
    inputs['simulation type'] = 'rel perm'
    mplbm.create_palabos_input_file(inputs)
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

    return


def process_and_plot_results(inputs):

    # Process data
    mplbm.create_pressure_data_file(inputs)
    mplbm.create_relperm_data_file(inputs)

    # Load and plot data
    sim_dir = inputs['input output']['simulation directory'] + '/'
    output_dir = inputs['input output']['output folder']
    Sw = np.loadtxt(f'{sim_dir + output_dir}data_Sw.txt')
    Pc = np.loadtxt(f'{sim_dir + output_dir}data_Pc.txt')
    krw = np.loadtxt(f'{sim_dir + output_dir}data_krw.txt')
    krnw = np.loadtxt(f'{sim_dir + output_dir}data_krnw.txt')

    plt.figure(figsize=[10,8])
    mplbm.plot_capillary_pressure_data(Sw, Pc)
    plt.savefig(sim_dir + 'pc_curve.png', dpi=300)

    plt.figure(figsize=[10,8])
    mplbm.plot_rel_perm_data(Sw, krw, krnw)
    plt.savefig(sim_dir + 'relperm_curve.png', dpi=300)

    plt.figure(figsize=[6, 10])
    mplbm.plot_pc_and_rel_perm(Sw, Pc, krw, krnw)
    plt.savefig(sim_dir + 'pc_and_relperm_curve.png', dpi=300)

    plt.show()

    return


# Some inputs and download geometry
input_folder = 'input/'
micromodel_name = 'rg_theta30_phi30_micromodel.raw'
data_type = 'uint8'
drp_url = 'https://www.digitalrocksportal.org/projects/65/images/71075/download/'
file_name = f'{input_folder}rg_theta30_phi30.raw'
download_geometry(file_name, drp_url)

# Rescale geometry
print("Rescaling geometry...")
geom = np.fromfile(file_name, dtype=data_type).reshape([501, 501, 501])
scaled_geom = mplbm.scale_geometry(geom, 0.599, data_type)
print(f"New geometry size = {scaled_geom.shape}")

# Create micromodel
print("Creating micromodel...")
micromodel = create_micromodel(scaled_geom, nx=200, ny=150, nz=5, slice_index=234, x_offset=0, y_offset=130)
geom_file = input_folder + micromodel_name
micromodel.flatten().tofile(geom_file)  # Note, use flatten() before writing! Otherwise, data not saved in correct order

# Parse inputs
input_file = 'input.yml'
inputs = mplbm.parse_input_file(input_file)  # Parse inputs
inputs['input output']['simulation directory'] = os.getcwd()  # Store current working directory

run_2_phase_sim(inputs)  # Run 2 phase sim
run_rel_perm_sim(inputs)  # Run rel perm
process_and_plot_results(inputs)  # Plot results

