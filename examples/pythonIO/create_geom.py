import sys
from argparse import Namespace

import numpy as np

sys.path.append('../../pre-processing') # I'll do this the proper way soon
from pore_utils import *




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
geom.num_slices = 1     # add n empty slices at the beggining and end of domain 
                        # for pressure bcs
geom.swapXZ     = True  # Swap x and z data if needed to ensure Palabos simulation in Z-direction              
geom.scale_2    = False # Double the grain (pore) size if needed to prevent single pixel throats
                        # for tight/ low porosity geometries                   

rock     = erase_regions(rock)
rock4sim = create_geom_edist(rock,geom) # provides an efficient geometry for simulation



