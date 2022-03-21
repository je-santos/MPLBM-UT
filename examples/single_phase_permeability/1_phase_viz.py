import numpy as np
import vedo as vd
import pyvista as pv
import os
import glob
import sys
import mplbm_utils as mplbm


def get_velocity_files(inputs):

    tmp_folder = inputs['input output']['output folder']

    # Get all the density files
    vel_files_regex = fr'{tmp_folder}vtk_vel*.vti'
    vel_files = glob.glob(vel_files_regex)

    # Sort for correct order
    vel_files_list = sorted(vel_files)

    return vel_files_list


def visualize_medium(inputs):

    input_folder = inputs['input output']['input folder']
    output_folder = inputs['input output']['output folder']
    nx = inputs['domain']['domain size']['nx']
    ny = inputs['domain']['domain size']['ny']
    nz = inputs['domain']['domain size']['nz']
    n_slices = inputs['domain']['inlet and outlet layers']

    medium = pv.read(f"{output_folder}PorousMedium000001.vti")
    grains = medium.get_array('tag').reshape([nz, ny, nx+n_slices*2])
    grains = grains[:, :, n_slices:nx+n_slices]
    grains = vd.Volume(grains).isosurface(0.5)

    return grains


def visualize_velocity(inputs, vel_file):

    nx = inputs['domain']['domain size']['nx']
    ny = inputs['domain']['domain size']['ny']
    nz = inputs['domain']['domain size']['nz']
    n_slices = inputs['domain']['inlet and outlet layers']

    vel_mesh = pv.read(vel_file)
    print(vel_mesh.array_names)
    vel_mesh = vel_mesh.get_array('velocityNorm').reshape([nz, ny, nx+n_slices*2])
    print(np.amax(vel_mesh), np.amin(vel_mesh))
    vel_mesh = vel_mesh[:, :, n_slices:nx]
    vel_thresholds = np.linspace(np.amin(vel_mesh), np.amax(vel_mesh), 20)
    vel = vd.Volume(vel_mesh).isosurface(threshold=vel_thresholds)

    return vel


# Get inputs
input_file = 'input.yml'
inputs = mplbm.parse_input_file(input_file)  # Parse inputs
inputs['input output']['simulation directory'] = os.getcwd()  # Store current working directory
inputs['domain']['inlet and outlet layers'] = 1
# Get density files
vel_files_list = get_velocity_files(inputs)

index = -1  # Choose last simulation output

# Setup plotter
vp = vd.Plotter(axes=9, bg='w', bg2='w', size=(1200,900), offscreen=False)

# visualize medium
grains = visualize_medium(inputs)
vp += grains.lighting('glossy').phong().c('seashell').opacity(0.2)

# visualize velocity
vel = visualize_velocity(inputs, vel_file=vel_files_list[index])
vp += vel.cmap('turbo').lighting('glossy').opacity(0.6).addScalarBar('Velocity [LBM Units]')  # .c('lightblue')

cam = dict(pos=(-85.32, 283.2, 150.6),
           focalPoint=(53.64, 39.77, 36.63),
           viewup=(0.2686, -0.2784, 0.9221),
           distance=302.6,
           clippingRange=(137.9, 487.7))

vp.show(camera=cam)
# vp.show(camera=cam).screenshot(f'velocity_viz.png', scale=1)

