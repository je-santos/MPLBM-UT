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
    f1_density = f1_density[:, :, n_slices:nx - 1]
    f1 = vd.Volume(f1_density).isosurface(1.5, largest=False)  # NW fluid
    f2 = vd.Volume(f1_density).isosurface(0.5, largest=False)  # W fluid

    return f1, f2


def calculate_droplet_curvature(inputs, rho_file):

    nx = inputs['domain']['domain size']['nx']
    ny = inputs['domain']['domain size']['ny']
    nz = inputs['domain']['domain size']['nz']
    n_slices = inputs['domain']['inlet and outlet layers']

    f1_mesh = pv.read(rho_file)
    f1_mesh_subset = f1_mesh.extract_subset([n_slices, nx-n_slices-1, 0, ny-1, 0, nz-1])
    f1_contour = f1_mesh_subset.contour(isosurfaces=[1.5], scalars='Density')
    curvature_array = f1_contour.curvature(curv_type='Mean', progress_bar=False)
    mean_curvature = np.mean(curvature_array)

    # p = pv.Plotter()
    # p.add_mesh(f1_contour, color='red')
    # p.show()

    return mean_curvature


# Get inputs
input_file = 'input.yml'
inputs = parse_input_file(input_file)  # Parse inputs
inputs['input output']['simulation directory'] = os.getcwd()  # Store current working directory

# Choose simulation to visualize
which_sim = 2
G_ads_f1_s1_labels = ['m0.4', '0', 'p0.4']
G_ads_f1_s1_list = np.array([-0.4, 0, 0.4])
tmp_folders = ['tmp_m0.4/', 'tmp_0/', 'tmp_p0.4/']
inputs['input output']['output folder'] = tmp_folders[which_sim]

# Get density files
rho_files_list = get_rho_files(inputs)

# Calculate mean curvature
index = -1  # Choose last simulation output
mean_curvature = calculate_droplet_curvature(inputs, rho_file=rho_files_list[index])
print(mean_curvature)

# Setup plotter
vp = vd.Plotter(axes=0, bg='w', bg2='w', size=(1200,900), offscreen=False)
# visualize medium
grains = visualize_medium(inputs)
vp += grains.lighting('glossy').phong().c('seashell')
# visualize fluids
f1, f2 = visualize_fluid(inputs, rho_file=rho_files_list[index])
vp += f1.subdivide(3).smooth().lighting('glossy').phong().c('lightblue')
vp += vd.shapes.Text2D(f'G_ads_f1_s1 = {G_ads_f1_s1_list[which_sim]}\nMean Curvature = {np.round(mean_curvature, decimals=3)}', pos='top-left')
cam = dict(pos=(90.66, -90.22, 73.85),
           focalPoint=(36.60, 37.96, 5.068),
           viewup=(-0.05450, 0.4543, 0.8892),
           distance=155.2,
           clippingRange=(46.16, 282.1))

vp.show(camera=cam).screenshot(f'droplet_viz_{G_ads_f1_s1_labels[which_sim]}.png', scale=1)

