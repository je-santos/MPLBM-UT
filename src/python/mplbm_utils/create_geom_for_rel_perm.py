import numpy as np
import glob
from argparse import Namespace
from .pore_utils import natural_sort, erase_regions, create_geom_edist, create_nw_fluid_mask


def create_geom_for_rel_perm(inputs):

    sim_dir = inputs['input output']['simulation directory']
    input_dir = inputs['input output']['input folder']
    output_dir = inputs['input output']['output folder']

    # Glob all .dat fluid files from simulation output folder
    rho_files_regex = fr'{sim_dir + "/" + output_dir}rho_f1*.dat'
    rho_files = glob.glob(rho_files_regex)
    rho_files = natural_sort(rho_files)

    # Save original geometry with 1 slice added for absolute permeability
    print("Creating original geometry for absolute permeability...")
    # read-in file
    data_type = inputs['geometry']['data type']
    Nx = inputs['geometry']['geometry size']['Nx']
    Ny = inputs['geometry']['geometry size']['Ny']
    Nz = inputs['geometry']['geometry size']['Nz']
    nx = inputs['domain']['domain size']['nx']
    ny = inputs['domain']['domain size']['ny']
    nz = inputs['domain']['domain size']['nz']
    geom_file_name = inputs['geometry']['file name']
    geom_name = inputs['domain']['geom name']
    geom_file = fr'{sim_dir}' + "/" + f'{input_dir + geom_file_name}'

    rock = np.fromfile(geom_file, dtype=data_type).reshape([Nx, Ny, Nz])
    # select a subset for simulation
    rock = rock[0:nz, 0:ny, 0:nx]
    geom = Namespace()
    geom.name = geom_name
    geom.print_size = True
    geom.add_mesh = False
    geom.num_slices = 2
    geom.swapXZ = inputs['domain']['swap xz']  # THIS NEEDS TO BE SAME AS ORIGINAL USER INPUT FOUND IN INPUT.YML
    geom.scale_2 = False
    geom.set_inlet_outlet_fluids = False  # For rel perms we don't want to add different inlet/outlet fluids

    rock, nw_fluid_mask = create_nw_fluid_mask(rock, geom)
    rock = rock/3
    rock = erase_regions(rock)
    # For absolute permeability, we want the whole geometry with only 1 fluid phase: nw_fluid_mask=[]
    rock4sim, geom_name = create_geom_edist(rock, geom, nw_fluid_mask=[])  # provides an efficient geometry for simulation
    inputs['domain']['geom name'] = geom_name
    rock4sim.flatten().tofile(sim_dir + '/' + input_dir + f'{geom_name}.dat')  # Save geometry

    # Calculate Pore Volume #
    #########################
    data = np.loadtxt(rho_files[0], dtype=np.float32)

    num_slices = inputs['domain']['inlet and outlet layers']
    # Update nx with inlet and outlet layers
    nx = inputs['domain']['domain size']['nx'] + num_slices * 2
    data = data.reshape([nx, ny, nz])  # Reshape
    data = data[num_slices:nx - num_slices, :, :]  # Remove slices

    domain = np.where(data > 0, data, -1)
    pores = np.where(domain < 0, domain, 0)
    pore_vol = len(np.where(pores == 0)[0])

    Sw_array = np.zeros([len(rho_files),])

    # For each file
    # 1) Read the file
    # 2) Reshape and remove slices
    # 3) Separate fluid 1 and save
    # 4) Separate fluid 2 and save

    for i in range(len(rho_files)):
        # i = 15
        print(f'Preparing geometry {i+1} of {len(rho_files)} for rel perm...')

        # 1) read the file #
        ####################
        data = np.loadtxt(rho_files[i], dtype=np.float32)
        # print(data.shape)

        # 2) Reshape, remove slices, calculate saturation #
        ###################################################
        num_slices = inputs['domain']['inlet and outlet layers']
        # Update nx with inlet and outlet layers
        nx = inputs['domain']['domain size']['nx'] + num_slices*2
        data = data.reshape([nx, ny, nz])  # Reshape
        data = data[num_slices:nx - num_slices, :, :]  # Remove slices

        # 3) Separate fluid 1 #
        #######################
        # Segmentation
        rho_1 = inputs['simulation']['rho_f1']
        tol = rho_1 - 1  # In order to account for interface and new boundaries that will be created
        rho1_mask = rho_1 - tol
        rho1_data = np.where(data >= rho1_mask, data, 2)  # Set to grain index
        rho1_data = np.where(data <= rho1_mask, rho1_data, 0)  # Set to pore index

        geom_name = 'f1_for_kr'
        geom = Namespace()
        geom.name = geom_name
        geom.print_size = False
        geom.add_mesh = False
        geom.num_slices = 2
        geom.swapXZ = False  # THIS NEEDS TO BE FALSE so user specified orientation remains the same
        geom.scale_2 = False
        geom.set_inlet_outlet_fluids = False  # For rel perms we don't want to add different inlet/outlet fluids

        rock = rho1_data
        # nw_fluid_mask=[] because we don't want to add nw fluid to the current config
        rock4sim, geom_name = create_geom_edist(rock, geom, nw_fluid_mask=[])

        rock4sim.flatten().tofile(sim_dir + '/' + input_dir + f'{geom_name}_{i+1}.dat')  # Save geometry

        # 4) Separate fluid 2 #
        #######################
        # Segmentation
        rho2_data = np.where((data > 0) & (data <= rho1_mask), data, 2)  # Set to grain index
        rho2_data = np.where(rho2_data >= rho1_mask, rho2_data, 0)  # Set to pore index

        # Calculate saturation of fluid 2 (Sw)
        f2_vol = len(np.where(rho2_data == 0)[0])
        Sw = f2_vol/pore_vol
        print(f'Sw = {Sw}')
        Sw_array[i] = Sw

        geom_name = 'f2_for_kr'
        geom = Namespace()
        geom.name = geom_name
        geom.print_size = False
        geom.add_mesh = False  # add a neutral-wet mesh at the end of the domain
        geom.num_slices = 2
        geom.swapXZ = False  # THIS NEEDS TO BE FALSE so user specified orientation remains the same
        geom.scale_2 = False
        geom.set_inlet_outlet_fluids = False  # For rel perms we don't want to add different inlet/outlet fluids

        rock = rho2_data
        # nw_fluid_mask=[] because we don't want to add nw fluid to the current config
        rock4sim, geom_name = create_geom_edist(rock, geom, nw_fluid_mask=[])
        rock4sim.flatten().tofile(sim_dir + '/' + input_dir + f'{geom_name}_{i+1}.dat')  # Save geometry

    # Save Sw array (Don't save very first Sw because that is an equilibration step, not used)
    # if inputs['simulation']['num pressure steps'] > 0:
    #     Sw_array[0] = 1  # If using pressure bcs, set very first pressure to 1. It's a test without pressure difference

    Sw_array = np.insert(Sw_array, 0, 1)

    np.savetxt(sim_dir + '/' + output_dir + 'data_Sw.txt', Sw_array)

    # Update inputs
    inputs['rel perm']['num_geoms'] = len(rho_files)*2 + 1

    return inputs

