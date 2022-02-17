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


def visualize_medium(inputs):

    input_folder = inputs['input output']['input folder']
    Nx = inputs['geometry']['geometry size']['Nx']
    Ny = inputs['geometry']['geometry size']['Ny']
    Nz = inputs['geometry']['geometry size']['Nz']
    n_slices = inputs['domain']['inlet and outlet layers']

    grains = np.fromfile(f"{input_folder}/flat_surface.raw", dtype='uint8').reshape([Nz, Ny, Nx])
    grains = np.transpose(grains, [0,1,2])
    # grains = medium.get_array('tag').
    print(grains.shape)
    # grains = grains[:, :, n_slices:nx - n_slices]
    # grains = grains[:, :, :]
    grains = vd.Volume(grains).isosurface(0.5)

    return grains


def visualize_fluid(inputs, rho_file):

    nx = inputs['domain']['domain size']['nx']
    ny = inputs['domain']['domain size']['ny']
    nz = inputs['domain']['domain size']['nz']
    n_slices = inputs['domain']['inlet and outlet layers']

    f1_mesh = pv.read(rho_file)
    f1_density = f1_mesh.get_array('Density').reshape([nz, ny, nx+n_slices*2])
    f1_density = f1_density[:, :, n_slices:nx - n_slices]
    f1 = vd.Volume(f1_density).isosurface(1.5, largest=True)  # NW fluid
    f2 = vd.Volume(f1_density).isosurface(0.5, largest=True)  # W fluid

    return f1, f2


# Get inputs
input_file = 'input.yml'
inputs = parse_input_file(input_file)  # Parse inputs
inputs['input output']['simulation directory'] = os.getcwd()  # Store current working directory

# Get density files
rho_files_list = get_rho_files(inputs)

# Setup plotter
vp = vd.Plotter(axes=3, bg='w', bg2='w', size=(1200,900), offscreen=False)
# visualize medium
grains = visualize_medium(inputs)
vp += grains.c('gray')
# visualize fluids
index = -1
f1, f2 = visualize_fluid(inputs, rho_file=rho_files_list[index])
vp += f1.c('red')
# vp += f2.c('blue')

vp.show()

