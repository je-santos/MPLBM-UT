import subprocess
import os
import sys
import mplbm_utils as mplbm


def download_geometry(filename, url):

    download_command = f'wget {url} -O {filename}'
    try:
        subprocess.run(download_command.split(' '))
    except FileNotFoundError:
        raise InterruptedError(f'wget was not found. Please make sure it is installed on your system.')
    return


def run_1_phase_sim(inputs):

    if inputs['simulation type'] == '2-phase':
        raise KeyError('Simulation type set to 2-phase...please change to 1-phase.')
    sim_directory = inputs['input output']['simulation directory']

    # 2) Create Palabos geometry
    print('Creating efficient geometry for Palabos...')
    mplbm.create_geom_for_palabos(inputs)

    # 3) Create simulation input file
    print('Creating input file...')
    mplbm.create_palabos_input_file(inputs)

    # 4) Run 1-phase simulationsss   
    print('Running 1-phase simulation...')
    num_procs = inputs['simulation']['num procs']
    input_dir = inputs['input output']['input folder']
    simulation_command = f"mpirun -np {num_procs} ../../src/1-phase_LBM/permeability {input_dir}1_phase_sim_input.xml"
    file = open(f'{sim_directory}/{input_dir}run_single_phase_sim.sh', 'w')
    file.write(f'{simulation_command}')
    file.close()
    simulation_command_subproc = f'bash {sim_directory}/{input_dir}run_single_phase_sim.sh'
    subprocess.run(simulation_command_subproc.split(' '))

    return


drp_url = 'www.digitalrocksportal.org/projects/65/images/71108/download/'
file_name = 'input/rg_theta60_phi10.raw'
download_geometry(file_name, drp_url)

input_file = 'input.yml'
inputs = mplbm.parse_input_file(input_file)  # Parse inputs
inputs['simulation']['fluid init'] = None 
inputs['input output']['simulation directory'] = os.getcwd()  # Store current working directory
run_1_phase_sim(inputs)  # Run 1 phase sim

