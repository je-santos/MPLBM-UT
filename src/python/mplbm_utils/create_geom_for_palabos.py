from argparse import Namespace
import numpy as np
from .pore_utils import erase_regions, create_geom_edist, create_nw_fluid_mask


def create_geom_for_palabos(inputs):

    sim_dir = inputs['input output']['simulation directory']
    input_dir = inputs['input output']['input folder']
    geom_file_name = inputs['geometry']['file name']
    data_type = inputs['geometry']['data type']
    geom_file = sim_dir + '/' + input_dir + geom_file_name
    Nx = inputs['geometry']['geometry size']['Nx']
    Ny = inputs['geometry']['geometry size']['Ny']
    Nz = inputs['geometry']['geometry size']['Nz']
    nx = inputs['domain']['domain size']['nx']
    ny = inputs['domain']['domain size']['ny']
    nz = inputs['domain']['domain size']['nz']
    geom_name = inputs['domain']['geom name']

    # read-in file
    rock = np.fromfile(geom_file, dtype=data_type).reshape([Nx, Ny, Nz])
    # select a subset for simulation
    rock = rock[0:nz, 0:ny, 0:nx]

    # geom inputs
    geom            = Namespace()
    geom.name       = geom_name
    geom.print_size = True
    geom.add_mesh   = False  # add a neutral-wet mesh at the end of the domain
    geom.num_slices = inputs['domain']['inlet and outlet layers']  # add n empty slices at the beginning and end of domain
                            # for pressure bcs
    geom.swapXZ     = inputs['domain']['swap xz']  # Swap x and z data if needed to ensure Palabos simulation in Z-direction
    geom.scale_2    = inputs['domain']['double geom resolution']  # Double the grain (pore) size if needed to prevent single pixel throats
                            # for tight/ low porosity geometries

    if inputs['simulation']['fluid init'] == 'geom':
        geom.set_inlet_outlet_fluids = True
        geom.inlet_fluid = inputs['simulation']['inlet fluid']
        geom.outlet_fluid = inputs['simulation']['outlet fluid']
    else:
        geom.set_inlet_outlet_fluids = False

    rock, nw_fluid_mask = create_nw_fluid_mask(rock, geom)
    rock = rock/3  # For proper erase regions and edist
    rock = erase_regions(rock)
    rock4sim, geom_name = create_geom_edist(rock, geom, nw_fluid_mask)  # provides an efficient geometry for simulation
    inputs['domain']['geom name'] = geom_name  # Update the geom name for later
    rock4sim.flatten().tofile(sim_dir + '/' + input_dir + f'{geom_name}.dat')  # Save geometry
    # np.savetxt(sim_dir + '/' + input_dir + f'{geom_name}.dat', rock4sim.flatten('C'))

    return

