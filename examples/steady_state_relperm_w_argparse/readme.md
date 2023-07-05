This Python script runs a steady-state relative permeability simulation given a specific saturation ratio condition. It works in a similar manner to the 'steady_state_relperm' example but leverages Python's `argparse` to define simulation conditions, offering more flexibility and code implementability for users.

The program works in the following steps:

1. Processes the 3D geometry of the model
2. Removes isolated pores
3. Initializes a simulation matrix with wetting and non-wetting fluids given a saturation ratio
4. Saves the simulation matrix
5. Plots the simulation matrix
6. Creates an input file for the simulation
7. Runs the 2-phase simulation

The script leverages multiprocessing and MPI for the simulation, making it suitable for large-scale simulations and high-performance computing environments.

## Simulation Parameters

The script accepts a wide range of parameters to configure the simulation. Here's a brief overview:

- `--simulation_dir`: The full path to the simulation directory
- `--input_dir`: The directory containing the input files
- `--output_dir`: The directory for storing output files
- `--file_name`: The name of the input geometry file
- `--data_type`: The data type of the input file
- `--raw_geometry_size`: The size of the geometry (Nx, Ny, Nz)
- `--sim_geometry_size`: The size of the simulation domain (nx, ny, nz)
- `--use_2d_stack`: Option to use 2D stacked geometry for simulation
- `--slice_no_2d`: The slice number in the y-direction for 2D image stacking
- `--periodic_boundary`: Option to use periodic boundary conditions
- `--inlet_and_outlet_layers`: Number of inlet and outlet layers (3-4 recommended for 2-phase simulation)
- `--add_mesh`: Option to add neutral mesh to the simulation
- `--swap_xz`: Option to swap x and z dimensions
- `--double_geom_resolution`: Option to double the geometry resolution
- `--wetting_f_saturation_ratio`: The saturation ratio of the wetting fluid
- `--num_procs`: The number of processors for the simulation
- `--restart_sim`: Option to continue from a previous saved state
- `--rho_f1`, `--rho_f2`: The densities of the two fluids
- `--force_f1`, `--force_f2`: The forces applied to the two fluids
- `--pressure_bc`: Option to use pressure difference boundary conditions
- `--minimum_radius`: The minimum radius for the simulation
- `--num_pressure_steps`: The number of pressure steps in the simulation
- `--fluid_init`: The fluid initialization method
- `--inlet_fluid`, `--outlet_fluid`: The type of fluid at the inlet and outlet respectively
- `--fluid_1_init`, `--fluid_2_init`: The initialization parameters for the two fluids
- `--Gc`: The Gc parameter for the simulation
- `--omega_f1`, `--omega_f2`: The omega parameters for the two fluids
- `--G_ads_f1_s1`, `--G_ads_f1_s2`, `--G_ads_f1_s3`, `--G_ads_f1_s4`: The adsorption parameters for the first fluid
- `--convergence`: The convergence threshold for the simulation
- `--convergence_iter`: The number of iterations between each convergence check
- `--max_iterations`: The maximum number of iterations per pressure step
- `--save_sim`: Option to save restart files
- `--save_iter`: The number of iterations between each save
- `--gif_iter`: The number of iterations between each GIF save
- `--vtk_iter`: The number of iterations between each VTK save
- `--rho_f2_vtk`: Option to save the VTK of the second fluid density
- `--print_geom`: Option to print the geometry at the start of the simulation
- `--print_stl`: Option to print the STL of the geometry at the start of the simulation
- `--pressure`: The pressure for one-phase simulations
- `--max_iterations_relperm`: The maximum number of iterations for relative permeability
- `--convergence_relperm`: The convergence for relative permeability
- `--save_vtks`: Option to save VTK files of the velocity

## Usage

The script can be used as follows:

```bash
python 2_phase_sim.py --simulation_type 2-phase --input_dir "input/" --output_dir "tmp/" --file_name "finney_pack.raw" --raw_geometry_size 500 500 500 --sim_geometry_size 128 10 128 --num_procs 6
