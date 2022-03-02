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

    output_folder = inputs['input output']['output folder']
    nx = inputs['domain']['domain size']['nx']
    ny = inputs['domain']['domain size']['ny']
    nz = inputs['domain']['domain size']['nz']
    n_slices = inputs['domain']['inlet and outlet layers']


    medium = pv.read(f"{output_folder}porousMedium.vti")
    medium = medium.get_array('tag').reshape([nz, ny, nx+n_slices*2])
    medium = medium[:, :, n_slices:nx]
    # medium = medium[:, :, :]
    grains = vd.Volume(medium).isosurface(0.5)

    return grains


def visualize_fluid(inputs, rho_file):

    nx = inputs['domain']['domain size']['nx']
    ny = inputs['domain']['domain size']['ny']
    nz = inputs['domain']['domain size']['nz']
    n_slices = inputs['domain']['inlet and outlet layers']

    f1_mesh = pv.read(rho_file)
    f1_density = f1_mesh.get_array('Density').reshape([nz, ny, nx+n_slices*2])
    f1_density = f1_density[:, :, n_slices:nx]
    f1 = vd.Volume(f1_density).isosurface([1.5,3])  # NW fluid
    f2 = vd.Volume(f1_density).isosurface([0.01,1.5])  # W fluid

    return f1, f2


def create_animation(inputs, rho_files_list, resolution_scale, restart):

    print('Creating animation...')

    sim_dir = inputs['input output']['simulation directory']
    output_dir = inputs['input output']['output folder']
    anim_dir = f'{sim_dir}/{output_dir}animation'
    anim_dir_exists = os.path.isdir(anim_dir)
    if anim_dir_exists == False:
        os.makedirs(anim_dir)

    # Loop through all rho files
    for i in range(len(rho_files_list)):

        if restart == False:
            current_image = f'{anim_dir}/image_{i}.png'
            current_image_exists = os.path.isfile(current_image)
            if current_image_exists == True:
                continue

        print(f'Image {i+1} of {len(rho_files_list)}...')

        # Create vedo plotter
        vp = vd.Plotter(axes=9, bg='w', bg2='w', size=(1200, 900), offscreen=False)

        # visualize medium
        grains = visualize_medium(inputs)
        vp += grains.smooth().lighting('glossy').c('seashell').opacity(0.1)

        # visualize fluids
        f1, f2 = visualize_fluid(inputs, rho_file=rho_files_list[i])
        vp += f1.c([191, 84, 0]).opacity(
            1).smooth()  # .lighting('plastic')  # subdivide(3).smooth().lighting('glossy').phong()
        vp += f2.c('lightblue').opacity(
            0.05).smooth()  # .lighting('plastic')  # subdivide(3).smooth().lighting('glossy').phong()

        cam = dict(pos=(245.1, -145.1, 139.7),
                   focalPoint=(58.32, 30.99, 36.04),
                   viewup=(-0.2395, 0.2914, 0.9262),
                   distance=276.9,
                   clippingRange=(146.4, 466.1))

        vp.show(camera=cam, interactive=False).screenshot(f'{anim_dir}/image_{i}.png', scale=resolution_scale)
        vp.close()

    return


# Get inputs
input_file = 'input.yml'
inputs = parse_input_file(input_file)  # Parse inputs
inputs['input output']['simulation directory'] = os.getcwd()  # Store current working directory

# update output directory
which_sim = 0
Sw_str = np.array(['10', '30', '50', '70', '90'])
inputs['input output']['output folder'] = f"tmp_{Sw_str[which_sim]}/"

# Get density files
rho_files_list = get_rho_files(inputs)

create_animation(inputs, rho_files_list, resolution_scale=1, restart=False)

exit()

index = -1  # Choose last simulation output

# Setup plotter
vp = vd.Plotter(axes=9, bg='w', bg2='w', size=(1200,900), offscreen=False)

# visualize medium
grains = visualize_medium(inputs)
vp += grains.smooth().lighting('glossy').c('seashell').opacity(0.1)

# visualize fluids
f1, f2 = visualize_fluid(inputs, rho_file=rho_files_list[index])
vp += f1.c([191, 84, 0]).opacity(1).smooth()  # .lighting('plastic')  # subdivide(3).smooth().lighting('glossy').phong()
vp += f2.c('lightblue').opacity(0.05).smooth()  # .lighting('plastic')  # subdivide(3).smooth().lighting('glossy').phong()

cam = dict(pos=(245.1, -145.1, 139.7),
           focalPoint=(58.32, 30.99, 36.04),
           viewup=(-0.2395, 0.2914, 0.9262),
           distance=276.9,
           clippingRange=(146.4, 466.1))

vp.show(camera=cam)
