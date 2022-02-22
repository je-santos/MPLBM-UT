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


def replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line):

    search_and_replace_command = 'sed -i "/^' + line_to_find_and_replace + r"/c\\" + replacement_line + '" ' + file_to_edit
    os.system(search_and_replace_command)

    return


def create_capillary_tubes(inputs):
    nx = 175
    ny = 175
    nz = 50
    x = np.arange(0,nx)
    y = np.arange(0,ny)
    domain = np.ones((nx, ny), dtype=np.uint8)  # A slice to add tubes to

    # Create tubes
    tube_radii = np.array([3, 4, 5, 6, 11, 16, 21])
    spacing = 5
    for i in range(0, len(tube_radii)):
        r = tube_radii[i]
        cx = spacing + np.max(tube_radii)
        cy = spacing*(i+1) + np.sum(tube_radii[0:i]*2) + r

        tube_mask = (x[np.newaxis, :] - cx)**2 + (y[:, np.newaxis] - cy)**2 < r**2
        domain[tube_mask] = 0  # if in the tube, set to 0 for empty space

    # Get smaller slice for a more efficient geometry
    domain = domain[:, 0:(spacing*2 + np.max(tube_radii)*2)]

    # Extend to 3D
    num_slices = nz
    tubes = np.repeat(domain[np.newaxis,:, :], num_slices, axis=0)
    tubes = tubes.transpose([1,2,0])  # Transpose to get aligned for Palabos properly

    # plt.figure()
    # plt.imshow(tubes[:,:,0])
    # plt.colorbar()
    # plt.show()
    #
    # import vedo as vd
    # vp = vd.Plotter(axes=3, bg='w', bg2='w', size=(1200, 900), offscreen=False)
    # # domain = vd.Volume(domain).isosurface(0.5)
    # tubes3d = vd.Volume(tubes).isosurface(threshold=0.5)
    # vp += tubes3d.c('orange')
    # vp.show()

    tubes.tofile(f"{inputs['input output']['input folder']}/capillary_tubes.raw")
    print(tubes.shape)
    nx = tubes.shape[2]
    ny = tubes.shape[1]
    nz = tubes.shape[0]

    inputs['geometry']['geometry size']['Nx'] = nz
    inputs['geometry']['geometry size']['Ny'] = ny
    inputs['geometry']['geometry size']['Nz'] = nx
    inputs['domain']['domain size']['nx'] = nx
    inputs['domain']['domain size']['ny'] = ny
    inputs['domain']['domain size']['nz'] = nz

    replace_line_in_file('input.yml', '    Nx', f'    Nx: {nx}')
    replace_line_in_file('input.yml', '    Ny', f'    Ny: {ny}')
    replace_line_in_file('input.yml', '    Nz', f'    Nz: {nz}')
    replace_line_in_file('input.yml', '    nx', f'    nx: {nx}')
    replace_line_in_file('input.yml', '    ny', f'    ny: {ny}')
    replace_line_in_file('input.yml', '    nz', f'    nz: {nz}')

    return inputs


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


def process_and_plot_results(inputs, run_name):

    # Process data
    create_pressure_data_file(inputs)

    # Load and plot data
    sim_dir = inputs['input output']['simulation directory'] + '/'
    output_dir = inputs['input output']['output folder']
    Sw = np.loadtxt(f'{sim_dir + output_dir}data_Sw.txt')
    Pc = np.loadtxt(f'{sim_dir + output_dir}data_Pc.txt')
    # krw = np.loadtxt(f'{sim_dir + output_dir}data_krw.txt')
    # krnw = np.loadtxt(f'{sim_dir + output_dir}data_krnw.txt')

    plot_capillary_pressure_data(Sw, Pc, Pc_label=f'{run_name}')

    return


def young_laplace_equation(radii, sigma, theta_w, height):

    # Calculate Young-Laplace solution
    Pc = 2*sigma*np.cos(theta_w)/radii
    total_volume = np.sum(radii**2 * np.pi) * height
    Sw = radii**2*np.pi / total_volume

    # Add entry pressure and saturation of 1
    entry_pressure = np.min(Pc)
    Pc = np.append(Pc, entry_pressure)
    Sw = np.append(Sw, 1)

    return Pc, Sw


input_file = 'input.yml'
inputs = parse_input_file(input_file)  # Parse inputs
inputs['input output']['simulation directory'] = os.getcwd()  # Store current working directory


output_folders = np.array(['1e-2', '1e-3'])  #, '1e-4', '1e-5'])

sigma = 0.15
theta_w = np.radians(180 - 156)
tube_height = 1
tube_radii = np.array([3, 4, 5, 6, 11, 16, 21])
Pc_yl, Sw_yl = young_laplace_equation(tube_radii, sigma, theta_w, tube_height)
plt.figure()

for i in range(len(output_folders)):
    inputs['input output']['output folder'] = output_folders[i] + '/'

    inputs = create_capillary_tubes(inputs)
    inputs['simulation']['convergence'] = float(output_folders[i])
    # run_2_phase_sim(inputs)  # Run 2 phase sim
    # create_geom_for_rel_perm(inputs)  # In order to get saturation files
    process_and_plot_results(inputs, run_name=output_folders[i])  # Plot results

plot_capillary_pressure_data(Sw_yl, Pc_yl, Pc_label='Young-Laplace', line_color='b')

plt.legend()
plt.show()
