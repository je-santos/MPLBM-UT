import numpy as np
from skimage import measure
import skimage.transform as skit
from scipy.ndimage.morphology import distance_transform_edt as edist
import re
import os
from edt import edt
import porespy as ps


def create_geom_edist(rock, args, nw_fluid_mask):

    if args.swapXZ:
        rock = rock.transpose([2, 1, 0])
        
    if args.scale_2:
        NotImplementedError('Feature not yet implemented')

    erock = edist(rock)

    # make sure all the BCs have bounce back nodes
    erock[0, :, :] = 1
    erock[:, 0, :] = 1
    erock[:, :, 0] = 1

    erock[-1, :, :] = 1
    erock[:, -1, :] = 1
    erock[:, :, -1] = 1

    # re open the pores
    erock[rock==0] = 0

    # Get the final matrix [0,1,2]
    erock[(erock>0)*(erock<2)] = 1
    erock[erock>1] = 2


    if args.add_mesh:
        NotImplementedError('Feature not yet implemented')
        
    if args.num_slices:
        erock = np.pad(erock, [(args.num_slices,args.num_slices), (0,0), (0,0)])

    if args.print_size:
        size = erock.shape
        geom_name = f'{args.name}_{size[0]}_{size[1]}_{size[2]}'
    else:
        geom_name = args.name

    # Save
    erock = erock.astype(np.int16)
    erock[erock == 0] = 2608  # pore space / w fluid
    erock[erock == 1] = 2609  # boundary
    erock[erock == 2] = 2610  # grains
    erock[nw_fluid_mask == 3] = 2611  # add nw fluid back in if needed
    erock = erock.astype(np.int16)

    return erock, geom_name


def create_nw_fluid_mask(rock, args):

    # Save indices for NW phase; can't do Euclidean distance properly with them.
    # Also need to take into account: (1) transpose and (2) number of slices added

    rock_tmp = np.copy(rock)
    if args.swapXZ:
        rock_tmp = rock_tmp.transpose([2, 1, 0])
    if args.num_slices:
        if args.set_inlet_outlet_fluids == True:
            if args.inlet_fluid == 'fluid 1':
                inlet_fluid = 3  # Set to fluid 1, NW phase
            elif args.inlet_fluid == 'fluid 2':
                inlet_fluid = 0  # Set to fluid 2, W phase
            else:
                raise ValueError('Please make sure inlet fluid set to "fluid 1" or "fluid 2"')
            if args.outlet_fluid == 'fluid 1':
                outlet_fluid = 3  # Set to fluid 1, NW phase
            elif args.outlet_fluid == 'fluid 2':
                outlet_fluid = 0  # Set to fluid 2, W phase
            else:
                raise ValueError('Please make sure outlet fluid set to "fluid 1" or "fluid 2"')
        else:
            inlet_fluid = 0
            outlet_fluid = 0
        rock_tmp = np.pad(rock_tmp, [(args.num_slices, args.num_slices),(0, 0),(0, 0)],
                          'constant', constant_values=(inlet_fluid, outlet_fluid))

        # There's currently an instability when NW fluid is right up against the geom...
        # The best solution will likely be to add a mesh up against the inlet and outlet?
        # n = args.num_slices  # inlet index
        # rock_tmp[n-2:n-1,:,:] = 0  # Layer of W fluid at start of geom for stability
        # rock_tmp[-n:-n+1, :, :] = 0  # Layer of W fluid at end of geom for stability

        # import matplotlib.pyplot as plt
        # plt.figure(figsize=[3,3])
        # plt.imshow(rock_tmp[:,:,40])
        # plt.colorbar()
        # plt.show()

    fluid_mask = np.where(rock_tmp == 3, rock_tmp, 0)  # Save NW whole block to preserve orientation
    rock = np.where(rock == 3, 0, rock)  # Finally, remove Nw phase from original image for rest of processing

    return rock, fluid_mask


def erase_regions(rock):
    # find connected-comps
    blobs_labels = measure.label(rock, background=1, connectivity=1)

    #vols = [np.sum(blobs_labels==label) for label in range(np.max(blobs_labels))]
    # it seems that label 1 is the largest comp, but gotta check

    # delete non-connected regions
    rock[blobs_labels>1] = 0
    
    return rock


def run_porespy_drainage(inputs, wetting_angle, voxel_size):
    # This function is just running PoreSpy drainage simulation on the image.
    # Much of this is from the drainage simulation example in the PoreSpy docs,
    # but there are a few modifications to make it compatible with MPLBM.

    sim_dir = inputs['input output']['simulation directory']
    input_dir = inputs['input output']['input folder']
    geom_file_name = inputs['geometry']['file name']
    data_type = inputs['geometry']['data type']
    geom_file = sim_dir + '/' + input_dir + geom_file_name
    Nx = inputs['geometry']['geometry size']['Nx']
    Ny = inputs['geometry']['geometry size']['Ny']
    Nz = inputs['geometry']['geometry size']['Nz']
    nx = inputs['domain']['domain size']['nx']
    ny = inputs['domain']['domain size']['ny']
    nz = inputs['domain']['domain size']['nz']
    swap_xz = inputs['domain']['swap xz']
    geom_name = inputs['domain']['geom name']

    image = np.fromfile(geom_file, dtype=data_type).reshape([Nx, Ny, Nz])
    image = image[0:nz, 0:ny, 0:nx]

    # Take into account user specified orientation
    if swap_xz == True:
        image = image.transpose([2, 1, 0])

    image = ~np.array(image, dtype=bool)  # Convert to bool and invert pores and grains for PoreSpy format
    inlets = np.zeros_like(image)  # Add inlets
    inlets[:,:,0] = True  # Make sure this is in XZ plane for correct orientation with Palabos
    sigma = 0.15  # This is the value of sigma found from MPLBM experiments (See Young-Laplace example)
    dt = edt(image)  # Get distance transform
    pc = -2 * sigma * np.cos(np.deg2rad(wetting_angle)) / (dt * voxel_size)  # Use Washburn equation for Pc values
    drn = ps.simulations.drainage(pc=pc, im=image, inlets=inlets, voxel_size=voxel_size, g=0)

    np.save(f'{sim_dir}/{input_dir}{geom_name}_pc_image.npy', drn.im_pc)
    np.save(f'{sim_dir}/{input_dir}{geom_name}_satn_image.npy', drn.im_satn)
    np.save(f'{sim_dir}/{input_dir}{geom_name}_pc_data.npy', drn.pc)
    np.save(f'{sim_dir}/{input_dir}{geom_name}_snwp_data.npy', drn.snwp)

    return drn.im_pc, drn.im_satn, drn.pc, drn.snwp


def convert_porespy_drainage_to_mplbm(image_satn, Snw):
    # Segment and convert from porespy to mplbm notation
    # 1) porespy all -1 --> mplbm 0 (pores not invaded)
    # 2) porespy all 0 --> mplbm 1 (grain)
    # 3) porespy saturation of interest and below --> mplbm 3 (nw phase, invaded)
    # 4) porespy above satn of interest --> mplbm 0 (w phase, not invaded)

    mplbm_geom = np.zeros_like(image_satn)
    mplbm_geom[image_satn == -1] = 0
    mplbm_geom[image_satn == 0] = 1
    mplbm_geom[(image_satn <= Snw) & (image_satn > 0)] = 3
    mplbm_geom[image_satn > Snw] = 0

    return mplbm_geom


def scale_geometry(geom, rescale_factor, data_type):

    geom_shape = np.array(geom.shape)
    scaled_geom_shape = np.array(geom.shape)*rescale_factor
    print(f'Scaling geometry from {geom_shape} to {scaled_geom_shape.astype(int)}')

    # Rescale geometry
    geom = skit.rescale(geom, rescale_factor, anti_aliasing=False,
                        order=0)  # order=0 means nearest neighbor interpolation (keeps image binary)

    # Ensure image has 0 as pore space and 1 as grains
    geom = edist(geom)
    geom[geom==0] = 0
    geom[geom>0] = 1

    # Change to specified data type
    geom = geom.astype(data_type)

    return geom


def natural_sort(l):

    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]

    return sorted(l, key=alphanum_key)


def find_line_in_file(file_name, line_to_match, data_to_add_line_index):

    file = open(file_name)
    data = np.array([])

    for line in file:
        if line_to_match in line:
            line_split = line.split()
            data = np.append(data, float(line_split[data_to_add_line_index]))

    file.close()

    return data


def replace_line_in_file(file_to_edit, line_to_find_and_replace, replacement_line):

    search_and_replace_command = 'sed -i "/^' + line_to_find_and_replace + r"/c\\" + replacement_line + '" ' + file_to_edit
    os.system(search_and_replace_command)

    return

