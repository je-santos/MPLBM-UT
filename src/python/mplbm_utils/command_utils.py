import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.ndimage import distance_transform_edt as edist, convolve, binary_dilation, label
from scipy.ndimage import label, sum as ndi_sum
from .pore_utils import scale_geometry

def load_geometry(args):
    geom = np.fromfile(f'{args.input_dir}{args.file_name}', dtype='uint8').reshape(args.raw_geometry_size)
    scaled_geom = scale_geometry(geom, 1, 'uint8')
    return scaled_geom

def stack_geometry(geom, args):
    if args.use_2d_stack==False:
        start_idx = [geom.shape[i]//2 - args.sim_geometry_size[i]//2 for i in range(3)]
        stop_idx = [start_idx[i] + args.sim_geometry_size[i] for i in range(3)]
        micromodel = geom[start_idx[0]:stop_idx[0], start_idx[1]:stop_idx[1], start_idx[2]:stop_idx[2]]
    else:
        slice_no = args.slice_no_2d
        stack_no = args.sim_geometry_size[1]
        geom_slice = geom[slice_no*args.sim_geometry_size[0]:(slice_no+1)*args.sim_geometry_size[0], 
                          slice_no, 
                          slice_no*args.sim_geometry_size[2]:(slice_no+1)*args.sim_geometry_size[2]]
        micromodel = np.repeat(geom_slice[:, np.newaxis, :], stack_no, axis=1)
    return micromodel

def process_geometry(args):
    geom = load_geometry(args)
    micromodel = stack_geometry(geom, args)
    return micromodel

def remove_isolated_pores(rock_matrix):
    pore_mask = (rock_matrix == 0)
    labeled_array, num_features = label(pore_mask)
    pore_sizes = ndi_sum(pore_mask, labeled_array, range(num_features + 1))
    isolated_mask = pore_sizes[labeled_array] == 1
    rock_matrix[isolated_mask] = 2
    return rock_matrix

def initialize_simulation_matrix(rock, wetting_saturation_ratio=0.5):
    erock = edist(rock)    
    # Ensure all the BCs have bounce back nodes
    erock[0, :, :] = 1
    erock[:, 0, :] = 1
    erock[: :, 0] = 1
    erock[-1, :, :] = 1
    erock[:, -1, :] = 1
    erock[:, :, -1] = 1    
    # Re open the pores
    erock[rock == 0] = 0    
    # Get the final matrix [0,1,2]
    erock[(erock > 0) & (erock < 2)] = 1
    erock[erock > 1] = 2    
    pore_voxel_num = np.count_nonzero(erock==0)
    wetting_voxel_num = int(pore_voxel_num*wetting_saturation_ratio)
    non_wetting_voxel_num = pore_voxel_num - wetting_voxel_num    
    # Fill non-wetting ratio from the farthest points from the grains 
    inverse_edist = edist(np.where(rock==1, 0, 1))
    sorted_indices = np.unravel_index(np.argsort(-inverse_edist, axis=None), inverse_edist.shape)    
    # Assuming 'arr_to_change' is the array you want to change
    for i in range(non_wetting_voxel_num):
        erock[sorted_indices[0][i], sorted_indices[1][i], sorted_indices[2][i]] = 3    
    # Assign new values to wetting phase, boundary, inside solid, non-wetting phase
    erock = erock.astype(np.int16)
    erock[erock == 0] = 2608  # pore space / wetting fluid
    erock[erock == 1] = 2609  # boundary
    erock[erock == 2] = 2610  # grains
    erock[erock == 3] = 2611  # non-wetting fluid
    erock = erock.astype(np.int16)
    return erock

def plot_simulation_matrix(sim_matrix, args):
    """
    Function to plot a simulation matrix.

    Args:
    sim_matrix (np.array): The simulation matrix.
    args (object): An object with attributes used for the simulation.
    """
    # Plot the simulation matrix
    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
    cmap = mcolors.ListedColormap(['red', 'blue', 'green', 'yellow'])
    bounds = [2608, 2609, 2610, 2611, 2612]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    img = ax.imshow(sim_matrix[:,args.sim_geometry_size[1]//2,:], cmap=cmap, norm=norm)
    # Create a colorbar
    cbar = fig.colorbar(img, cmap=cmap, norm=norm, boundaries=bounds, ticks=[2608.5, 2609.5, 2610.5, 2611.5])
    cbar.ax.set_yticklabels(['Pore Space / Wetting Fluid', 'Boundary', 'Grains', 'Non-wetting Fluid'])  
    # If the directory doesn't exist, create it
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    # Save the plot
    plt.savefig(f'{args.output_dir}Initial_fluid_configuration.png')
    plt.close(fig)  # Close the figure
    
def create_two_phase_input_file_2(args):
    # Generate input xml file
    geom_name = f'{args.input_dir}{args.file_name[:-4]}_micromodel.dat'
    nx = args.sim_geometry_size[0]
    ny = args.sim_geometry_size[1]
    nz = args.sim_geometry_size[2]
    num_layers = args.inlet_and_outlet_layers
    domain_size = args.sim_geometry_size
    periodic = args.periodic_boundary
    
    # Parse i/o inputs
    input_folder = args.input_dir
    output_folder = args.output_dir
    input_xml_file = '2_phase_sim_input.xml'
    
    restart_sim = args.restart_sim
    
    if args.pressure_bc == True:
        minimum_radius = args.minimum_radius
        num_pc_steps = args.num_pressure_steps
    else:
        minimum_radius = 1
        num_pc_steps = 0
    
    load_fluid_type = args.fluid_init
    if load_fluid_type == 'geom':
        load_fluid_from_geom = True
    else:
        load_fluid_from_geom = False
    
    # If the directory doesn't exist, create it
    if not os.path.exists(args.input_dir):
        os.makedirs(args.input_dir)
        
    # Create/open input file
    file = open(f'{input_folder}{input_xml_file}', 'w+')
    file.write('<?xml version="1.0" ?>\n\n')  # Write xml header
    
    # Restart sim?
    file.write(f'<load_savedstated> {restart_sim} </load_savedstated>\n\n')
    
    # Write geometry section
    file.write('<geometry>\n')
    # Geometry name
    file.write(f'\t<file_geom> {geom_name} </file_geom>\n')
    # Geometry size
    file.write(f'\t<size> <x> {nx} </x> <y> {ny} </y> <z> {nz} </z> </size>\n')
    # Periodicity
    file.write(f'\t<per>\n')
    file.write(f'\t\t<fluid1> <x> {periodic[0]} </x> <y> {periodic[1]} </y> <z> {periodic[2]} </z> </fluid1>\n')
    file.write(f'\t\t<fluid2> <x> {periodic[0]} </x> <y> {periodic[1]} </y> <z> {periodic[2]} </z> </fluid2>\n')
    file.write(f'\t</per>\n')
    file.write('</geometry>\n\n')
    
    # Write initial position of fluids
    file.write(f'<init>\n')
    file.write(f'\t<fluid_from_geom> {load_fluid_from_geom} </fluid_from_geom>\n')
    file.write(f'\t<fluid1>\n')
    file.write(f'\t\t <x1> {args.fluid_1_init[0]} </x1> <y1> {args.fluid_1_init[1]} </y1> <z1> {args.fluid_1_init[2]} </z1>\n')
    file.write(f'\t\t <x2> {args.fluid_1_init[3]} </x2> <y2> {args.fluid_1_init[4]} </y2> <z2> {args.fluid_1_init[5]} </z2>\n')
    file.write(f'\t</fluid1>\n')
    file.write(f'\t<fluid2>\n')
    file.write(f'\t\t <x1> {args.fluid_2_init[0]} </x1> <y1> {args.fluid_2_init[1]} </y1> <z1> {args.fluid_2_init[2]} </z1>\n')
    file.write(f'\t\t <x2> {args.fluid_2_init[3]} </x2> <y2> {args.fluid_2_init[4]} </y2> <z2> {args.fluid_2_init[5]} </z2>\n')
    file.write(f'\t</fluid2>\n')
    file.write('</init>\n\n')
    
    # Write fluid data
    file.write('<fluids>\n')   
    file.write(f'\t<Gc> {args.Gc} </Gc>\n')
    file.write(f'\t<omega_f1> {args.omega_f1} </omega_f1>\n')
    file.write(f'\t<omega_f2> {args.omega_f2} </omega_f2>\n')
    file.write(f'\t<force_f1> {args.force_f1} </force_f1>\n')
    file.write(f'\t<force_f2> {args.force_f2} </force_f2>\n')
    file.write(f'\t<G_ads_f1_s1> {args.G_ads_f1_s1} </G_ads_f1_s1>\n')
    file.write(f'\t<G_ads_f1_s2> {args.G_ads_f1_s2} </G_ads_f1_s2>\n')
    file.write(f'\t<G_ads_f1_s3> {args.G_ads_f1_s3} </G_ads_f1_s3>\n')
    file.write(f'\t<G_ads_f1_s4> {args.G_ads_f1_s4} </G_ads_f1_s4>\n')    
    file.write(f'\t<rho_f1> {args.rho_f1} </rho_f1>\n')
    file.write(f'\t<rho_f2> {args.rho_f2} </rho_f2>\n')    
    file.write(f'\t<pressure_bc> {args.pressure_bc} </pressure_bc>\n')
    file.write(f'\t<rho_f1_i> {args.rho_f1} </rho_f1_i>\n')
    file.write(f'\t<rho_f2_i> {args.rho_f2} </rho_f2_i>\n')
    file.write(f'\t<num_pc_steps> {num_pc_steps} </num_pc_steps>\n')
    file.write(f'\t<min_radius> {minimum_radius} </min_radius>\n')
    file.write(f'\t<rho_d> 0.06 </rho_d>\n')    
    file.write('</fluids>\n\n')    
    # Write output section
    file.write('<output>\n')    
    file.write(f'\t<out_folder> {output_folder} </out_folder>\n')
    file.write(f'\t<save_it> {args.save_iter} </save_it>\n')
    file.write(f'\t<save_sim> {args.save_sim} </save_sim>\n')
    file.write(f'\t<convergence> {args.convergence} </convergence>\n')
    file.write(f'\t<it_max> {args.max_iterations} </it_max>\n')
    file.write(f'\t<it_conv> {args.convergence_iter} </it_conv>\n')
    file.write(f'\t<it_gif> {args.gif_iter} </it_gif>\n')
    file.write(f'\t<rho_vtk> {args.rho_f2_vtk} </rho_vtk>\n')
    file.write(f'\t<it_vtk> {args.vtk_iter} </it_vtk>\n')
    file.write(f'\t<print_geom> {args.print_geom} </print_geom>\n')
    file.write(f'\t<print_stl> {args.print_stl} </print_stl>\n')    
    file.write('</output>')    
    file.close()    
    return    