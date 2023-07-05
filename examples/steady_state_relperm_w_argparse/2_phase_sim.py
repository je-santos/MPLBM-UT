import os
import subprocess
import sys
import mplbm_utils as mplbm
import argparse as args
import warnings
from numba import NumbaDeprecationWarning

# Ignore the specific NumbaDeprecationWarning warning
warnings.filterwarnings("ignore", category=NumbaDeprecationWarning)

parser = args.ArgumentParser()

# Simulation type
parser.add_argument('--simulation_type', type=str, default="2-phase", help='Simulation type: "1-phase" or "2-phase".')

# Input output parameters
parser.add_argument('--simulation_dir', type=str, default="", help='Full path to simulation directory.')
parser.add_argument('--input_dir', type=str, default="input/", help='Input folder.')
parser.add_argument('--output_dir', type=str, default="tmp/", help='Output folder.')

# Geometry parameters
#parser.add_argument('--file_name', type=str, default="rg_theta30_phi30.raw", help='Name of the input geometry file.')
parser.add_argument('--file_name', type=str, default="finney_pack.raw", help='Name of the input geometry file.')
parser.add_argument('--data_type', type=str, default="uint8", help='Data type.')
#parser.add_argument('--raw_geometry_size', type=int, nargs=3, default=[501, 501, 501], help='Geometry size (Nx, Ny, Nz).')
parser.add_argument('--raw_geometry_size', type=int, nargs=3, default=[500, 500, 500], help='Geometry size (Nx, Ny, Nz).')

# Domain parameters
#parser.add_argument('--geom_name', type=str, default="finney_pack", help='Name of .dat file, rename from original if you\'d like. Do not include the file extension.')
parser.add_argument('--sim_geometry_size', type=int, nargs=3, default=[128, 10, 128], help='Domain size (nx, ny, nz).')
parser.add_argument('--use_2d_stack', type=bool, default=True, help='Used 2d stacked geometry for simulation')
parser.add_argument('--slice_no_2d', type=int, default=0, help='slice number in y direction in case of stacking 2D image for sim.')
parser.add_argument('--periodic_boundary', type=bool, nargs=3, default=[True, True, True], help='Periodic boundary (x, y, z).')
parser.add_argument('--inlet_and_outlet_layers', type=int, default=0, help='Inlet and outlet layers, 3-4 recommended for 2-phase sim.')
parser.add_argument('--add_mesh', type=bool, default=False, help='Add neutral mesh, by default False.')
parser.add_argument('--swap_xz', type=bool, default=False, help='Swap xz, False by default.')
parser.add_argument('--double_geom_resolution', type=bool, default=False, help='Double geom resolution, False by default.')

# Simulation parameters
parser.add_argument('--wetting_f_saturation_ratio', type=float, default=0.5, help='Wetting fluid saturation ratio')
parser.add_argument('--num_procs', type=int, default=6, help='Number of processes.')
parser.add_argument('--restart_sim', type=bool, default=False, help='Set to true if you would like to continue from a previous saved state.')
parser.add_argument('--rho_f1', type=float, default=2, help='Fluid density 1.')
parser.add_argument('--rho_f2', type=float, default=2, help='Fluid density 2.')
parser.add_argument('--force_f1', type=float, default=1e-4, help='Force f1.')
parser.add_argument('--force_f2', type=float, default=1e-4, help='Force f2.')
parser.add_argument('--pressure_bc', type=bool, default=False, help='Pressure difference boundary conditions.')
parser.add_argument('--minimum_radius', type=int, default=3, help='Minimum radius.')
parser.add_argument('--num_pressure_steps', type=int, default=1, help='Number of pressure steps.')
parser.add_argument('--fluid_init', type=str, default='geom', help='Fluid initialization method.')
parser.add_argument('--inlet_fluid', type=str, default='fluid 2', help='Inlet fluid.')
parser.add_argument('--outlet_fluid', type=str, default='fluid 2', help='Outlet fluid.')
parser.add_argument('--fluid_1_init', type=int, nargs=6, default=[1, 2, 1, 75, 1, 75], help='Fluid 1 initialization parameters.')
parser.add_argument('--fluid_2_init', type=int, nargs=6, default=[3, 150, 1, 75, 1, 75], help='Fluid 2 initialization parameters.')
parser.add_argument('--Gc', type=float, default=0.9, help='Gc parameter.')
parser.add_argument('--omega_f1', type=float, default=1, help='Omega f1 parameter.')
parser.add_argument('--omega_f2', type=float, default=1, help='Omega f2 parameter.')
parser.add_argument('--G_ads_f1_s1', type=float, default=-0.4, help='G ads f1 s1 parameter.')
parser.add_argument('--G_ads_f1_s2', type=float, default=0, help='G ads f1 s2 parameter.')
parser.add_argument('--G_ads_f1_s3', type=float, default=0, help='G ads f1 s3 parameter.')
parser.add_argument('--G_ads_f1_s4', type=float, default=0, help='G ads f1 s4 parameter.')
parser.add_argument('--convergence', type=float, default=1e-4, help='Convergence threshold.')
parser.add_argument('--convergence_iter', type=int, default=1000, help='How often to check for convergence.')
parser.add_argument('--max_iterations', type=int, default=500000, help='Max iterations per Pc step.')
parser.add_argument('--save_sim', type=bool, default=True, help='Save restart files.')
parser.add_argument('--save_iter', type=int, default=20000, help='How often to save restart files.')
parser.add_argument('--gif_iter', type=int, default=2000, help='How often to save gifs.')
parser.add_argument('--vtk_iter', type=int, default=2000, help='How often to save vtk files.')
parser.add_argument('--rho_f2_vtk', type=bool, default=False, help='When True, saves rho f1 and f2 vtks. If False, only saves rho f1 vtk.')
parser.add_argument('--print_geom', type=bool, default=True, help='Create vtk of geometry at beginning.')
parser.add_argument('--print_stl', type=bool, default=False, help='Create stl of geometry at beginning.')

# Relative permeability parameters for 1-phase sims
parser.add_argument('--pressure', type=float, default=0.0005, help='Pressure.')
parser.add_argument('--max_iterations_relperm', type=int, default=5000000, help='Max iterations for rel perm.')
parser.add_argument('--convergence_relperm', type=float, default=1e-6, help='Convergence for rel perm.')
parser.add_argument('--save_vtks', type=bool, default=True, help='Save velocity vtks.')
args = parser.parse_args()

rock = mplbm.process_geometry(args) # Load geometry for simulation
rock = mplbm.remove_isolated_pores(rock) # Remove isolated pores in the geometry
sim_matrix = mplbm.initialize_simulation_matrix(rock, wetting_saturation_ratio=args.wetting_f_saturation_ratio) # Initialize wetting/non-wetting fluid given saturation ratio
sim_matrix.flatten().tofile(f'{args.input_dir}{args.file_name[:-4]}_micromodel.dat') # Save the simulation matrix
mplbm.plot_simulation_matrix(sim_matrix, args) # Plot the simulation matrix
mplbm.create_two_phase_input_file_2(args) # Save input file for simulation
print('Running 2-phase simulation...')
simulation_command = f"mpirun -np {args.num_procs} ../../src/2-phase_LBM/ShanChen {args.input_dir}2_phase_sim_input.xml"
with open(f'{args.input_dir}run_shanchen_sim.sh', 'w') as file:
    file.write(simulation_command)
simulation_command_subproc = f'bash {args.input_dir}run_shanchen_sim.sh'
subprocess.run(simulation_command_subproc.split(' '))

