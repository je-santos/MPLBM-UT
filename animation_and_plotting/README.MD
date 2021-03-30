# **Animation and Plotting Code**
The Python code in this folder was designed for two main tasks:
1) Creating 3D animations/gifs of the simulation and 3D visualizations of the fluid configuration files
2) Creating plots for the relative permeability and (for unsteady-state simulations) the capillary pressure

The idea is that you can copy this whole folder into your simulation directory when you want to visualize something. In order for the code to work, you need the following dependencies:
- Python 3.7.x in a venv or Anaconda environment with the following modules:
  - PyVista (and related dependencies)
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



