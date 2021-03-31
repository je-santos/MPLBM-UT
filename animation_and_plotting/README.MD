# **Animation and Plotting Code**
The Python code in this folder was designed for two main tasks:
1) Creating 3D animations/gifs of the simulation and 3D visualizations of the fluid configuration files
2) Creating plots for the relative permeability and (for unsteady-state simulations) the capillary pressure

The idea is that you can copy this whole folder into your simulation directory when you want to visualize something. In order for the code to work, you need the following dependencies:
- Python 3.7.x in a venv or Anaconda environment with the following modules:
  - [PyVista](https://github.com/pyvista/pyvista) (and related dependencies)
  - Matplotlib
  - Numpy
  - Imageio

## Creating Relative Permeability and Capillary Pressure Plots
1) Once your simulation is complete, open `plot_results.py`.
2) Copy and paste your non-wetting saturation (f1 saturation from `info_4kr.txt`), wetting phase relative permeability, non-wetting phase relative permeability, and (if needed) capillary pressure data into the following variables respectively:
    - `Snw`
    - `krw`
    - `krnw`
    - `rho2_start` is the starting density of the simulation and `rho2_end` is the last density of the simulaiton. The actual capillary pressure, `Pc`, will be calculated in dimensionless units.
3) Edit the plotting details and customize as you need! If you do not need capillary pressure plots, feel free to comment those lines out.
4) When you run the script (unedited), three plots will be saved:
    - Relative permeability plot as `rel_perm_graph.png`
    - Capillary pressure plot as `pc_graph.png`
    - A combined relative permeability and capillary pressure plot as `rel_perm_and_pc_graph.png` 

## Creating 3D Visualizations and Animations
These steps describe how to create an animation, but a single visualization is also described once you get to step 3.
1) Once your simulation is complete, create the following directories inside the `animation_and_plotting` folder:
    - `cd your_simulation_directory/animation_and_plotting/`
    - `mkdir /fluid_data/rho_f1_data/`
    - `mkdir /fluid_data/rho_f2_data/`
2) From `/tmp` in your simulation directory, move the 'rho_f1' and 'rho_f2' vti files into their respecitve directories
    - `cd your_simulation_directory/`
    - `mv tmp/rho_f1*.vti animation_and_plotting/fluid_data/rho_f1_data/`
    - `mv tmp/rho_f2*.vti animation_and_plotting/fluid_data/rho_f2_data/`
3) Open `pyvista_lbm_animation.py` and at the top of the file, change `n_slices_start` and `n_slices_end` to reflect the number of additional slices added at the beginning and end of the simultion domain respectively.
    - Familiarity with the PyVista module will be helpful when making edits and customizing this code. You may also need to adjust the camera position; comment it out from being set in the code if you would like PyVista to automatically try and figure it out.
    - Leaving this code as is will work for creating animations, but it can be easily modified to visualize a single fluid configuration by doing the following:
      - Comment out the line reading `i = sys.argv[1]` and add below it `i = image_number_of_interest`. For example, if you want the first pair of density files visualized, set `i = 0`. As of now, you do have to manually count (index 0) to the file of interest if you would like a single image.
      - In the line reading `p = pv.Plotter(off_screen=True)`, delete the argument `off_screen=True`
      - From here, you can simply run this as a python script and the visualization should appear
4) Open `create_animation.sh` and edit the for loop to reflect the number of images needed to loop through all the vti files.
    - **IMPORTANT:** This script is set up with multithreading for speed! Please adjust the `num_proc` variable so that your system does not run out of RAM! Each thread will open 2 vti files (one for f1 and one for f2) which can be quite large. To be safe, try with `num_proc=1` or `num_proc=2` at first and adjust to get faster results. To get a feel for the orders of magnitude of memory needed, a 400^3 simulation with `num_proc=4` used about ~15 GB RAM (each thread using about 3.75 GB RAM).
    - Note that the images are processed in pairs of rho_f1 and rho_f2, and the bash script is zero indexed. So, if you have 10 rho_f1 files and 10 rho_f2 files and you want an animation of all the files, enter `{0..9}` in the for loop
5) Run the animation script: `bash create_animation.sh`
