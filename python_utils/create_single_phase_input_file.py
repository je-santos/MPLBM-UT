def create_single_phase_input_file(inputs):

    if inputs['simulation type'] == '1-phase':
        input_file_name = '1_phase_sim_input.xml'
    elif inputs['simulation type'] == '2-phase':
        input_file_name = '2_phase_sim_input.xml'

    geom_name = inputs['domain']['geom name']
    nx = inputs['domain']['domain size']['nx']
    ny = inputs['domain']['domain size']['ny']
    nz = inputs['domain']['domain size']['nz']
    num_layers = inputs['domain']['inlet and outlet layers']
    domain_size = [nx+(2*num_layers), ny, nz]
    per_x = inputs['domain']['periodic boundary']['x']
    per_y = inputs['domain']['periodic boundary']['y']
    per_z = inputs['domain']['periodic boundary']['z']
    periodic = [per_x, per_y, per_z]
    io_folders = [inputs['input output']['input folder'], inputs['input output']['output folder']]
    num_geoms_or_sims = inputs['simulation']['num geoms']
    pressure = inputs['simulation']['pressure']
    max_iter = inputs['simulation']['max iterations']
    convergence = inputs['simulation']['convergence']
    save_vtks = inputs['simulation']['save vtks']

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

    # Create/open input file
    file = open(f'{input_folder}{input_file_name}', 'w')
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
