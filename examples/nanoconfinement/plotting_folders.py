from vedo import *
import numpy as np
import vedo as vd
from hdf5storage import loadmat
import pyvista as pv

import glob
import os

import matplotlib.pyplot as plt    
    
def legoplot(im, im_name, folder):
    
    im_size = im.shape[0]
    
    vol = vd.Volume(im)
    vol.addScalarBar3D()
    
    #lego = vol.legosurface(vmin=1, vmax=2)
    lego = vol.legosurface(vmin=0.9, vmax=2)
    
    
    if im_size == 256:
        cam = dict(pos=(-304, 570, 645),
                   focalPoint=(148, 76.8, 108),
                   viewup=(0, 1.00, 0),
                   distance=859,
                   clippingRange=(360, 1.37e+3))
    elif im_size == 480:
        cam = dict(pos=(-499, 1.02e+3, 1.13e+3),
                   focalPoint=(303, 140, 178),
                   viewup=(0, 1.00, 0),
                   distance=1.52e+3,
                   clippingRange=(556, 2.45e+3))
    
    vd.show(lego, camera=cam, interactive=False).screenshot(
                                 f'{folder}{im_name}_lego.png',scale=4).close()
        
    
def sliceplot(im, im_name, folder):    
    grid = pv.UniformGrid()
    grid.dimensions = np.array(im.shape) + 1
    grid.origin = (0, 0, 0)   # The bottom left corner of the data set
    grid.spacing = (1, 1, 1)  # These are the cell sizes along each axis
    grid.cell_arrays["values"] = im.flatten(order="F")  # Flatten the array!
    slices = grid.slice_orthogonal(x=1,y=1,z=1)
    cmap = plt.cm.get_cmap("plasma", 2)
    sargs = dict(height=0.00025, vertical=True, position_x=0.05, position_y=0.05)

    slices.plot(cmap=cmap, background='white', 
                scalar_bar_args=sargs, 
                screenshot= f'{folder}{im_name}_slices.png',notebook=False)
    




folder = 'binary_volumes_real/'
mat_files = [os.path.basename(x) for x in glob.glob(folder + '/*.mat')]
#bin_image = loadmat(folder + mat_files[0])['bin']
ggg
for file in mat_files:
    bin_image = loadmat(folder + file)['bin']
    
    im_name = file.replace('.mat','')
    legoplot(bin_image,  im_name, folder)
    sliceplot(bin_image, im_name, folder)