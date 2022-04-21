import numpy as np
import vedo as vd
import pyvista as pv
import glob
import mplbm_utils as mplbm
from create_gif_and_mp4 import *


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


def create_animation(inputs, rho_files_list, cam, resolution_scale, restart):

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
        # vp += grains.smooth().lighting('glossy').c('seashell').opacity(0.1)

        # visualize fluids
        f1, f2 = visualize_fluid(inputs, rho_file=rho_files_list[i])
        vp += f1.c([191, 84, 0]).opacity(0.9).smooth()  # .lighting('plastic')  # subdivide(3).smooth().lighting('glossy').phong()
        vp += f2.c('lightblue').opacity(0.2).smooth()  # .lighting('plastic')  # subdivide(3).smooth().lighting('glossy').phong()

        vp.show(camera=cam, interactive=False).screenshot(f'{anim_dir}/image_{i}.png', scale=resolution_scale)
        vp.close()

    anim_dir = inputs['input output']['output folder'] + 'animation/'
    save_name = inputs['domain']['geom name']
    create_gif(anim_dir, save_name)
    create_mp4(anim_dir, save_name, speed_factor=0.5)  # speed_factor = 1 means no slow down or speed up

    return


# Get inputs
input_file = 'input.yml'
inputs = mplbm.parse_input_file(input_file)  # Parse inputs
inputs['input output']['simulation directory'] = os.getcwd()  # Store current working directory

# update output directory
which_sim = 1
Snw_str = np.array([1,2,3,4,5,6,7,8])
inputs['input output']['output folder'] = f"tmp_{Snw_str[which_sim]}/"

# Get density files
rho_files_list = get_rho_files(inputs)

# set camera
cam = dict(pos=(328.4, -216.7, 190.2),
           focalPoint=(74.81, 58.93, 34.88),
           viewup=(-0.2424, 0.2973, 0.9235),
           distance=405.4,
           clippingRange=(195.6, 646.7))

# For animation
# create_animation(inputs, rho_files_list, cam, resolution_scale=1, restart=True)

# For single frame
vp = vd.Plotter(axes=9, bg='w', bg2='w', size=(1200, 900), offscreen=False)
index = -1  # Choose last simulation output
grains = visualize_medium(inputs)  # visualize medium
f1, f2 = visualize_fluid(inputs, rho_file=rho_files_list[index])  # Viz for fluid
grains = grains.c('seashell').opacity(0.1).smooth()
f1 = f1.c([191, 84, 0]).opacity(0.95).smooth()
f2 = f2.c('lightblue').opacity(0.2).smooth()
# vp += grains
vp += f1
vp += f2
vp.show(camera=cam, interactive=True)  # .screenshot(f'{anim_dir}/image_{i}.png', scale=resolution_scale)


# Visualize PoreSpy Drainage (visualize the initial conditions)
porespy_satn_geom = np.load(f'input/finney_pack_satn_image.npy')

Snw = np.load(f'input/finney_pack_snwp_data.npy')
remove_ind = np.where(Snw <= 0.01)[0]
Snw = np.delete(Snw, remove_ind)  # skip grains, uninvaded marker, and Snw values less than 1%
satn_threshold = 0.01  # Remove values closer than 1% together to remove redundant points
Snw = np.delete(Snw, np.argwhere(np.ediff1d(Snw) <= satn_threshold) + 1)

porespy_to_mplbm = mplbm.convert_porespy_drainage_to_mplbm(porespy_satn_geom, Snw[which_sim])

vp = vd.Plotter(axes=1, bg='w', bg2='w', size=(1200, 900), offscreen=False)
mplbm_geom = vd.Volume(porespy_to_mplbm)
nw = mplbm_geom.isosurface(threshold=[0, 3])
grain = mplbm_geom.isosurface(threshold=[0, 1])
vp += nw.c('red').opacity(1)
vp += grain.c('shell').opacity(0.1)
vp.show(camera=cam)

