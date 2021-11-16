import sys
from argparse import Namespace

import numpy as np

sys.path.append('../../pre-processing') # I'll do this the proper way soon
from pore_utils import *
from create_single_phase_input_file import *


# read-in file
rock = np.fromfile('input/grid_25.bin', dtype='int8').reshape([480,480,480])/3

# select a subset for simultion
subsize = 250
rock    = rock[0:subsize,0:subsize,0:subsize]

# geom inputs
geom            = Namespace()
geom.name       = 'carbonate4Palabos'
geom.print_size = True
geom.add_mesh   = False # add a neutral-wet mesh at the end of the domain
geom.num_slices = 1     # add n empty slices at the beginning and end of domain
                        # for pressure bcs
geom.swapXZ     = True  # Swap x and z data if needed to ensure Palabos simulation in Z-direction              
geom.scale_2    = False # Double the grain (pore) size if needed to prevent single pixel throats
                        # for tight/ low porosity geometries                   

rock     = erase_regions(rock)
rock4sim = create_geom_edist(rock,geom) # provides an efficient geometry for simulation

# Create input file
input_file_name = "input_perm.xml"
num_slices = geom.num_slices*2
geom_name = f'{geom.name}_{subsize + num_slices}_{subsize}_{subsize}'
domain_size = [subsize + num_slices, subsize, subsize]  # Domain size x, y, z
periodic = ["false", "true", "true"]  # Periodicity in x, y, z
input_folder = "input/"
output_folder = "tmp/"
io_folders = [input_folder, output_folder]
num_sims = 1
sim_pressure = 0.0005
sim_max_iter = 1000000
sim_convergence = 0.005
sim_settings = [num_sims, sim_pressure, sim_max_iter, sim_convergence]
save_vtks = "true"
create_single_phase_input_file(input_file_name, geom_name, domain_size, periodic, io_folders, sim_settings, save_vtks)


