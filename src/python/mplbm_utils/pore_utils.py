import numpy as np
from skimage import measure
import skimage.transform as skit
from scipy.ndimage.morphology import distance_transform_edt as edist
import re
import os


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
    print('Original unique values: ', np.unique(rock))
    nw_indices = np.where(rock == 3)[0]

    rock_tmp = np.copy(rock)
    if args.swapXZ:
        rock_tmp = rock_tmp.transpose([2, 1, 0])
    if args.num_slices:
        rock_tmp = np.pad(rock_tmp, [(args.num_slices, args.num_slices), (0, 0), (0, 0)])

    fluid_mask = np.where(rock_tmp == 3, rock_tmp, 0)  # Save NW whole block to preserve orientation
    rock[nw_indices] = 0  # Finally, remove Nw phase from original image for rest of processing

    return rock, fluid_mask


def erase_regions(rock):
    # find connected-comps
    blobs_labels = measure.label(rock, background=1, connectivity=1)

    #vols = [np.sum(blobs_labels==label) for label in range(np.max(blobs_labels))]
    # it seems that label 1 is the largest comp, but gotta check

    # delete non-connected regions
    rock[blobs_labels>1] = 0
    
    return rock


def scale_geometry(geom, rescale_factor, data_type):

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

