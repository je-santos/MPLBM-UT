import pyvista.examples as examples
import pyvista as pv
import numpy as np
import os
import imageio
import re
import time
import sys


t = time.time()  # start time

i = sys.argv[1]  # iteration count from bash script
i = int(i)
# print(i)
# print(type(i))

cwd = os.getcwd()
f1_storage_dir = cwd + '/fluid_data/rho_f1_data/'
f2_storage_dir = cwd + '/fluid_data/rho_f2_data/'
f1_list = os.listdir(f1_storage_dir)
f2_list = os.listdir(f2_storage_dir)
# Sort lists for correct order
f1_list = sorted(f1_list)
f2_list = sorted(f2_list)

# Slices added 
n_slices_start = 4
n_slices_end = 5

# Progress update
counter = i + 1
total_images = len(f1_list)

print("Processing images...")
# for i in range(len(f1_list)):
image_time = time.time()  # start time

# Create plotter with screen off after first iteration.
p = pv.Plotter()  # off_screen=True

f1_mesh = pv.read(f1_storage_dir + f1_list[i])
f2_mesh = pv.read(f2_storage_dir + f2_list[i])

# Get dimensions
x_dim = f1_mesh.dimensions[0]
y_dim = f1_mesh.dimensions[1]
z_dim = f1_mesh.dimensions[2]

# Clip added slices
f1_mesh = f1_mesh.extract_subset([n_slices_start,x_dim-n_slices_end-1, 0,y_dim-1, 0,z_dim-1])
f2_mesh = f2_mesh.extract_subset([n_slices_start,x_dim-n_slices_end-1, 0,y_dim-1, 0,z_dim-1])

f1_contours = f1_mesh.contour(isosurfaces=np.linspace(0.5,2,5), scalars='Density')
f2_contours = f2_mesh.contour(isosurfaces=np.linspace(0.5,2,5), scalars='Density')
#f1_contours = f1_contours.extract_geometry().smooth(n_iter=200)
#f2_contours = f2_contours.extract_geometry().smooth(n_iter=200)

p.set_background(color='w')
p.add_mesh(f1_mesh.outline(), color="k")
p.add_mesh(f1_contours, opacity=0.1, color='r')
p.add_mesh(f2_contours, opacity=0.1, color=(132/255, 210/255, 232/255))
p.show_grid(color='k')


# camera = [(609.531061491956, 605.531061491956, 605.531061491956),
#           (128.5, 124.5, 124.5),
#           (0.0, 0.0, 1.0)]

camera = [(638.3648346692544, 638.1479077135746, 529.3073566909535),
          (127.14199155932997, 139.84240950470522, 99.75470392571974),
          (-0.34745354867374556, -0.3818154900807505, 0.8564420371796433)]

cpos = p.show(cpos=camera, screenshot='gif_images/image_0' + str(i) + '.png')
# print(cpos)

elapsed_image_time = time.time() - image_time
progress = 'Image frame ' + str(counter) + ' of ' + str(total_images) + ' processed in ' + str(np.around(elapsed_image_time,3)) + ' s.'
counter = counter + 1
print(progress)  # To see camera position


# def natural_sort(l):
#     convert = lambda text: int(text) if text.isdigit() else text.lower()
#     alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
#     return sorted(l, key=alphanum_key)
#
#
# print("Creating GIF...")
# images = []
# filenames = os.listdir(cwd + '/gif_images/')
# filenames = natural_sort(filenames)
#
# for filename in filenames:
#     images.append(imageio.imread('gif_images/' + filename))
# imageio.mimsave('pore_throat_LBM.gif', images)
#
# print("Done!")

elapsed_time = time.time() - t
time_message = "Everything took " + str(np.around(elapsed_time, 3)) + ' s'
# print(time_message)
