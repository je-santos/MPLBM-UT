import pyvista as pv
import vedo as vd
import glob
import numpy as np

# Slices added
n_slices_start = 4
n_slices_end = 4

# Geometry size
nz = 75
ny = 75
nx = 60 + 2*n_slices_start

tmp_folder = '/home/alexgigliotti/Documents/ut/research/mplbm-ut/MPLBM-UT/examples/python_2_phase_workflow/tmp/'

f1_files_regex = fr'{tmp_folder}rho_f1*.vti'
f1_files = glob.glob(f1_files_regex)

f2_files_regex = fr'{tmp_folder}rho_f2*.vti'
f2_files = glob.glob(f2_files_regex)

# Sort lists for correct order
f1_list = sorted(f1_files)
f2_list = sorted(f2_files)

# print(f1_list)
# print(f2_list)

index = 3

medium = pv.read(tmp_folder + 'porousMedium.vti')
grains = medium.get_array('tag').reshape([nz, ny, nx])
grains = grains[:, :, n_slices_start:nx-n_slices_end]
grains = grains[:,:,:]
print(grains.shape)
grains = vd.Volume(grains).isosurface(1)

plt = vd.Plotter(axes=9, bg='w', bg2='w', size=(1200,900), offscreen=False)  # Create a vedo plotter; axes=9 is just an outline
plt += grains.c("gray").opacity(0.4)
cam = dict(pos=(179, 313, 136),
           focalPoint=(45.0, 49.5, 44.0),
           viewup=(-0.163, -0.252, 0.954),
           distance=309,
           clippingRange=(158, 501))

plt.show(camera=cam, interactive=True)  # .screenshot(f'animation_and_plotting/gif_images/test{index}.png', scale=4)
plt.close()
# exit()

# for index in range(len(f1_list)):
f1_mesh = pv.read(f1_list[index])
f1_density = f1_mesh.get_array('Density').reshape([nz, ny, nx])
f1_density = f1_density[:, :, n_slices_start:nx-n_slices_end]
oil = vd.Volume(f1_density).isosurface(1, largest=False)  # .addCurvatureScalars(method=1)  # .smoothLaplacian(niter=1000)

f2_mesh = pv.read(f2_list[index])
f2_density = f2_mesh.get_array('Density').reshape([nz, ny, nx])
f2_density = f2_density[:, :, n_slices_start:nx-n_slices_end]
water = vd.Volume(f2_density).isosurface(0.405, largest=False)  # .smoothLaplacian(niter=1000)

axes_opts = dict(xtitle='', ytitle='', ztitle='', yzGrid=False,
                 xyGrid=False, xyFrameLine=2, yzFrameLine=2, zxFrameLine=2)  # Can create advanced axes options
plt = vd.Plotter(axes=9, bg='w', bg2='w', size=(1200,900), offscreen=False)  # Create a vedo plotter; axes=9 is just an outline

# plt += grains.c("gray").opacity(0.4)
plt += water.c([132, 210, 232]).opacity(0.3)
plt += oil.c([191, 84, 0]).opacity(1)

cam = dict(pos=(179, 313, 136),
           focalPoint=(45.0, 49.5, 44.0),
           viewup=(-0.163, -0.252, 0.954),
           distance=309,
           clippingRange=(158, 501))

plt.show(camera=cam, interactive=True)  # .screenshot(f'animation_and_plotting/gif_images/test{index}.png', scale=4)
plt.close()
