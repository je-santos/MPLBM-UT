import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import glob
import mplbm_utils as mplbm


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

    return inputs


def process_and_plot_results(inputs):

    # Process data
    # create_pressure_data_file(inputs)
    mplbm.create_relperm_data_file(inputs)

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
    mplbm.plot_rel_perm_data(Sw, krw, krnw)
    # plt.yscale('log')
    plt.savefig(sim_dir + 'relperm_curve.png', dpi=300)

    # plt.figure(figsize=[6, 10])
    # plot_pc_and_rel_perm(Sw, Pc, krw, krnw)
    # plt.savefig(sim_dir + 'pc_and_relperm_curve.png', dpi=300)

    # plt.show()

    return


def process_and_save_porespy_geometries(inputs, image_satn):

    # Process saturation values
    Snw = np.unique(image_satn)
    remove_ind = np.where(Snw <= 0.01)[0]
    Snw = np.delete(Snw, remove_ind)  # skip grains, uninvaded marker, and Snw values less than 1%
    satn_threshold = 0.01  # Remove values closer than 1% together to remove redundant points
    Snw = np.delete(Snw, np.argwhere(np.ediff1d(Snw) <= satn_threshold) + 1)
    print(f"Snw values for sims: {Snw}")

    # Convert to MPLBM conventions and save geometries
    sim_counter = np.arange(1, len(Snw) + 1, 1)
    input_folder = inputs["input output"]["input folder"]
    geom_name = inputs['domain']['geom name']
    for i in range(len(Snw)):
        mplbm_geom = mplbm.convert_porespy_drainage_to_mplbm(image_satn, Snw[i])
        mplbm_geom = mplbm_geom.astype(inputs["geometry"]["data type"])
        mplbm_geom.flatten().tofile(f"{input_folder}{geom_name}_{sim_counter[i]}.raw")

    # Change geometry inputs since now reading the drainage geometry
    # Swap nx and nz because C++ saving (np.tofile()) flipped from Python reading
    inputs['geometry']['geometry size']['Nx'] = inputs['domain']['domain size']['nz']
    inputs['geometry']['geometry size']['Ny'] = inputs['domain']['domain size']['ny']
    inputs['geometry']['geometry size']['Nz'] = inputs['domain']['domain size']['nx']
    # 'swap xz' needs to be opposite of original input since geom altered for sim already
    inputs['domain']['swap xz'] = ~inputs['domain']['swap xz']

    return inputs, Snw, sim_counter


def organize_2_phase_outputs_for_relperm(inputs, sim_counter):
    # Organize outputs 2 phase outputs #
    ####################################
    # Gather final steady state 2 phase configs for rel perm
    # Create output folder
    sim_dir = inputs['input output']['simulation directory']
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
    for i in range(len(sim_counter)):
        output_dir = f"{sim_dir}/tmp_{sim_counter[i]}/"
        # Glob all .dat fluid files from simulation output folder
        rho_files_regex = fr'{output_dir}rho_f1*.dat'
        rho_files = glob.glob(rho_files_regex)
        rho_files = mplbm.natural_sort(rho_files)
        last_rho_config = rho_files[-1]
        cp_command = f'cp {last_rho_config} {rel_perm_dir}/rho_f1_{i + 1}.dat'
        subprocess.run(cp_command.split(' '))

    return inputs

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
scaled_geom = mplbm.scale_geometry(geom, 0.75, 'uint8')
scaled_geom.tofile('input/finney_pack_half_scale.raw')

# Parse Inputs #
################
input_file = 'input.yml'
inputs = mplbm.parse_input_file(input_file)  # Parse inputs
inputs['input output']['simulation directory'] = os.getcwd()  # Store current working directory
sim_dir = inputs['input output']['simulation directory']

# Initial conditions from PoreSpy #
###################################
voxel_size = 1e-6  # m/voxel side
wetting_angle = 180  # in degrees
print('Running PoreSpy drainage sim...')
image_pc, image_satn, pc, snwp = mplbm.run_porespy_drainage(inputs, wetting_angle, voxel_size)
inputs, Snw, sim_counter = process_and_save_porespy_geometries(inputs, image_satn)

# Run 2 phase sims #
####################
geom_name = inputs['domain']['geom name']
for i in range(len(sim_counter)):

    # Create output folder
    output_dir = f'{sim_dir}/tmp_{sim_counter[i]}'
    output_dir_exists = os.path.isdir(output_dir)
    if output_dir_exists == False:
        os.makedirs(output_dir)

    # Update to correct inputs
    inputs["input output"]["output folder"] = f"tmp_{sim_counter[i]}/"
    inputs["geometry"]["file name"] = f"{geom_name}_{sim_counter[i]}.raw"
    inputs["domain"]["geom name"] = f"{geom_name}_{sim_counter[i]}"

    print("\n\n--------------------------------------------")
    print(f"Running sim {i+1} of {len(sim_counter)}")
    print(f"PoreSpy Snw = {Snw[i]}\n")
    run_2_phase_sim(inputs)  # Run 2 phase sim

# Run rel perm and process results #
####################################
inputs = organize_2_phase_outputs_for_relperm(inputs, sim_counter)
inputs = run_rel_perm_sim(inputs)  # Run rel perm
process_and_plot_results(inputs)  # Plot results
