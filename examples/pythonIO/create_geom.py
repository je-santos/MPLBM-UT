from argparse import Namespace

import numpy as np
import matplotlib.pyplot as plt
from skimage import measure



def create_geom_edist(rock, args):
    from scipy.ndimage.morphology import distance_transform_edt as edist
    
    if args.swapXZ:
        rock = rock.transpose([2, 1, 0])
        
    if args.scale_2:
        NotImplementedError('Feature not yet implemented')

    erock   = edist(rock)
    
    # make sure all the BCs have BB nodes
    erock[0,:,:]  = 1
    erock[:,0,:]  = 1
    erock[:,:,0]  = 1
    
    erock[-1,:,:] = 1
    erock[:,-1,:] = 1
    erock[:,:,-1] = 1
    
    # re open the pores
    erock[rock==0]=0;
    
    # Get the final matrix [0,1,2]
    erock[(erock>0)*(erock<2)] = 1
    erock[erock>1] = 2
    
    if args.add_mesh:
        NotImplementedError('Feature not yet implemented')
        
    if args.num_slices:
        erock = np.pad(erock,[(args.num_slices,args.num_slices),(0,0),(0,0)])
        
    
    if args.print_size:
        size = erock.shape
        geom_name = f'{args.name}_{size[0]}_{size[1]}_{size[2]}.dat'
    else:
        geom_name = args.name
        
    np.savetxt(f'input/{geom_name}.dat', erock.flatten())
    return erock



# read-in file
rock = np.fromfile('input/grid_25.bin', dtype='int8').reshape([480,480,480])/3

# select a subset for simultion
subsize = 250
rock    = rock[0:subsize,0:subsize,0:subsize]

# find connected-comps
blobs_labels = measure.label(rock, background=1, connectivity=1)

#vols = [np.sum(blobs_labels==label) for label in range(np.max(blobs_labels))]
# it seems that label 1 is the largest comp, but gotta check

# delete non-connected regions
rock[blobs_labels>1] = 0


geom            = Namespace()
geom.name       = 'carbonate4Palabos'
geom.print_size = True
geom.add_mesh   = False # add a neutral-wet mesh at the end of the domain
geom.num_slices = 1     # add n empty slices at the beggining and end of domain 
                        # for pressure bcs
geom.swapXZ     = True    # Swap x and z data if needed to ensure Palabos simulation in Z-direction              
geom.scale_2    = False   # Double the grain (pore) size if needed to prevent single pixel throats
                          # for tight/ low porosity geometries                   

rock4sim = create_geom_edist(rock,geom) # provides an efficient geometry for simulation





