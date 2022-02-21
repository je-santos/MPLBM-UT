import matplotlib.pyplot as plt
import numpy as np
import vedo as vd
import pyvista as pv
import os
import glob
import sys
sys.path.append('../../python_utils/')
from parse_input_file import *


def get_rho_files(inputs):

    tmp_folder = inputs['input output']['output folder']

    # Get all the density files
    f1_files_regex = fr'{tmp_folder}rho_f1*.vti'
    f1_files = glob.glob(f1_files_regex)

    # Sort for correct order
    rho_files_list = sorted(f1_files)

    return rho_files_list


def get_slice_of_medium(inputs, slice):

    input_folder = inputs['input output']['input folder']
    nx = inputs['domain']['domain size']['nx']
    ny = inputs['domain']['domain size']['ny']
    nz = inputs['domain']['domain size']['nz']
    n_slices = inputs['domain']['inlet and outlet layers']

    grains = np.fromfile(f"{input_folder}{inputs['geometry']['file name']}", dtype='uint8').reshape(nz, ny, nx)
    grains = np.transpose(grains, [0, 1, 2])
    grains = grains[slice, :, :]

    # medium = pv.read(f"{input_folder}{inputs['geometry']['file name']}")
    # grains = medium.get_array('tag').reshape([nz, ny, nx+n_slices*2])
    # grains = grains[slice, :, n_slices:nx+n_slices]


    return grains


def get_slice_of_fluid(inputs, rho_file, slice):

    nx = inputs['domain']['domain size']['nx']
    ny = inputs['domain']['domain size']['ny']
    nz = inputs['domain']['domain size']['nz']
    n_slices = inputs['domain']['inlet and outlet layers']

    f1_mesh = pv.read(rho_file)
    f1_density = f1_mesh.get_array('Density').reshape([nz, ny, nx+n_slices*2])
    f1_density = f1_density[slice, :, n_slices:nx+n_slices]

    return f1_density


def plot_sim_contours(inputs, rho, medium):

    nx = inputs['domain']['domain size']['nx']
    ny = inputs['domain']['domain size']['ny']

    # Plotting setup
    x = np.arange(0, nx, 1)
    y = np.arange(0, ny, 1)
    X, Y = np.meshgrid(x, y)

    # Plotting
    plt.contourf(X, Y, medium, levels=[0.5, 1], alpha=1, colors='gray')
    plt.contourf(X, Y, rho, levels=[1, 3], alpha=1, colors='red')
    plt.contourf(X, Y, rho, levels=[0, 1], alpha=1, colors='lightblue')

    plt.axis('equal')
    plt.xlim([0,nx])
    plt.ylim([0,ny])

    return


def create_animation(inputs, rho_files_list):

    print('Creating animation...')

    sim_dir = inputs['input output']['simulation directory']
    output_dir = inputs['input output']['output folder']
    anim_dir = f'{sim_dir}/{output_dir}animation'
    anim_dir_exists = os.path.isdir(anim_dir)
    if anim_dir_exists == False:
        os.makedirs(anim_dir)

    grains = get_slice_of_medium(inputs, slice=3)

    for i in range(len(rho_files_list)):
        print(f'Image {i+1} of {len(rho_files_list)}...')
        plt.figure()
        rho = get_slice_of_fluid(inputs, rho_file=rho_files_list[i], slice=3)
        plot_sim_contours(inputs, rho, grains)
        plt.axis('off')
        plt.savefig(f'{anim_dir}/image_{i}.png', dpi=300)
        plt.close()



    return


def update_micromodel_inputs(inputs, rescale):

    nx = inputs['domain']['domain size']['nx']
    ny = inputs['domain']['domain size']['ny']
    nz = inputs['domain']['domain size']['nz']

    # Change inputs to match micromodel geometry
    micromodel_name = inputs['domain']['geom name'] + '_micromodel.raw'
    inputs['geometry']['file name'] = micromodel_name
    inputs['domain']['geom name'] = inputs['domain']['geom name'] + '_micromodel'
    inputs['geometry']['geometry size']['Nx'] = nz
    inputs['geometry']['geometry size']['Ny'] = int(ny * rescale)
    inputs['geometry']['geometry size']['Nz'] = int(nx * rescale)
    inputs['domain']['domain size']['nx'] = int(nx * rescale)
    inputs['domain']['domain size']['ny'] = int(ny * rescale)
    inputs['domain']['domain size']['nz'] = nz

    return inputs

# Get inputs
input_file = 'input.yml'
inputs = parse_input_file(input_file)  # Parse inputs
inputs['input output']['simulation directory'] = os.getcwd()  # Store current working directory

# Update inputs for micromodel
scale = 0.5
inputs = update_micromodel_inputs(inputs, scale)

# Get density files
rho_files_list = get_rho_files(inputs)

create_animation(inputs, rho_files_list)

# # Get slices for contour viz
# rho = get_slice_of_fluid(inputs, rho_file=rho_files_list[-1], slice=3)
# grains = get_slice_of_medium(inputs, slice=3)
#
# # Setup plotter
# plt.figure()
# plot_sim_contours(inputs, rho, grains)
# plt.show()

