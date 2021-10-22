def create_single_phase_input_file(input_file_name, geom_name, domain_size, periodic, io_folders, sim_settings, save_vtks):
    """
    :param input_file_name: Name of input file to open or create
    :type input_file_name: str
    :param geom_name: The name of the .dat file that contains the entire geometry (do not include ".dat" extension)
    :type geom_name: str
    :param domain_size: Size (in voxels) of the simulation domain. ([x, y, z])
    :type domain_size: list
    :param periodic: Turn periodicity true or false in x, y, z. (list of strings, "True" or "False")
    :type periodic: list
    :param io_folders: List containing the input folder and output folder as strings
    :type io_folders: list
    :param sim_settings: List containing - number of geometries/sims (int), pressure (float), max iterations (int), convergence (float)
    :type sim_settings: list
    :param save_vtks: "True" or "False", save vtks for medium and velocity.
    :type save_vtks: str
    :return:
    """

    # Parse geometry inputs
    x_size = domain_size[0]
    y_size = domain_size[1]
    z_size = domain_size[2]
    periodic_x = periodic[0]
    periodic_y = periodic[1]
    periodic_z = periodic[2]

    # Parse i/o inputs
    input_folder = io_folders[0]
    output_folder = io_folders[1]

    # Parse simulation inputs
    num_geoms_or_sims = sim_settings[0]
    pressure = sim_settings[1]
    max_iter = sim_settings[2]
    convergence = sim_settings[3]

    # Create/open input file
    file = open(f'{input_file_name}', 'w')
    file.write('<?xml version="1.0" ?>\n\n')  # Write xml header

    # Write geometry section
    file.write('<geometry>\n')
    # Geometry name
    file.write(f'\t<file_geom> {geom_name} </file_geom>\n')
    # Geometry size
    file.write(f'\t<size> <x> {x_size} </x> <y> {y_size} </y> <z> {z_size} </z> </size>\n')
    # Periodicity
    file.write(f'\t<per> <x> {periodic_x} </x> <y> {periodic_y} </y> <z> {periodic_z} </z> </per>\n')
    file.write('</geometry>\n\n')

    # Write i/o section
    file.write('<folder>\n')
    # Input folder name
    file.write(f'\t<in_f> {input_folder} </in_f>\n')
    # Output folder name
    file.write(f'\t<out_f> {output_folder} </out_f>\n')
    file.write('</folder>\n\n')

    # Write simulation section
    file.write('<simulations>\n')
    # Number of sims/geometries
    file.write(f'\t<num> {num_geoms_or_sims} </num>\n')
    # Pressure
    file.write(f'\t<press> {pressure} </press>\n')
    # Max simulation iterations
    file.write(f'\t<iter> {max_iter} </iter>\n')
    # Convergence
    file.write(f'\t<conv> {convergence} </conv>\n')
    # Save vtks
    file.write(f'\t<vtk_out> {save_vtks} </vtk_out>\n')
    file.write('</simulations>')

    file.close()

    return
