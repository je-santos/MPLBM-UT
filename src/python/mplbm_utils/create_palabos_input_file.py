def create_palabos_input_file(inputs):

    if inputs['simulation type'] == '1-phase':
        input_file_name = '1_phase_sim_input.xml'
        create_one_phase_input_file(inputs, input_file_name)

    elif inputs['simulation type'] == '2-phase':
        input_file_name = '2_phase_sim_input.xml'
        create_two_phase_input_file(inputs, input_file_name)

    elif inputs['simulation type'] == 'rel perm':
        input_file_name = 'relperm_input.xml'
        create_relperm_input_file(inputs, input_file_name)

    return


def create_one_phase_input_file(inputs, input_file_name):

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
    # perm_model = inputs['simulation']['perm model']
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
    # # Permeability model
    # file.write(f'\t<perm_model> {perm_model} </perm_model>\n')
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


def create_two_phase_input_file(inputs, input_file_name):
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

    restart_sim = inputs['simulation']['restart sim']

    rho_f1 = inputs['simulation']['rho_f1']
    rho_f2 = inputs['simulation']['rho_f2']

    force_f1 = inputs['simulation']['force_f1']
    force_f2 = inputs['simulation']['force_f2']

    pressure_bc = inputs['simulation']['pressure bc']
    if pressure_bc == True:
        minimum_radius = inputs['simulation']['minimum radius']
        num_pc_steps = inputs['simulation']['num pressure steps']
    else:
        minimum_radius = 1
        num_pc_steps = 0

    load_fluid_type = inputs["simulation"]["fluid init"]
    if load_fluid_type == 'geom':
        load_fluid_from_geom = True
    else:
        load_fluid_from_geom = False

    fluid1_x1 = inputs['simulation']['fluid 1 init']['x1']
    fluid1_x2 = inputs['simulation']['fluid 1 init']['x2']
    fluid1_y1 = inputs['simulation']['fluid 1 init']['y1']
    fluid1_y2 = inputs['simulation']['fluid 1 init']['y2']
    fluid1_z1 = inputs['simulation']['fluid 1 init']['z1']
    fluid1_z2 = inputs['simulation']['fluid 1 init']['z2']
    fluid2_x1 = inputs['simulation']['fluid 2 init']['x1']
    fluid2_x2 = inputs['simulation']['fluid 2 init']['x2']
    fluid2_y1 = inputs['simulation']['fluid 2 init']['y1']
    fluid2_y2 = inputs['simulation']['fluid 2 init']['y2']
    fluid2_z1 = inputs['simulation']['fluid 2 init']['z1']
    fluid2_z2 = inputs['simulation']['fluid 2 init']['z2']

    Gc = inputs['simulation']['fluid data']['Gc']
    omega_f1 = inputs['simulation']['fluid data']['omega_f1']
    omega_f2 = inputs['simulation']['fluid data']['omega_f2']
    G_ads_f1_s1 = inputs['simulation']['fluid data']['G_ads_f1_s1']
    G_ads_f1_s2 = inputs['simulation']['fluid data']['G_ads_f1_s2']
    G_ads_f1_s3 = inputs['simulation']['fluid data']['G_ads_f1_s3']
    G_ads_f1_s4 = inputs['simulation']['fluid data']['G_ads_f1_s4']

    convergence = inputs['simulation']['convergence']
    convergence_iter = inputs['simulation']['convergence iter']
    max_iter = inputs['simulation']['max iterations']
    save_sim = inputs['simulation']['save sim']
    save_iter = inputs['simulation']['save iter']
    gif_iter = inputs['simulation']['gif iter']
    vtk_iter = inputs['simulation']['vtk iter']
    rho_f2_vtk = inputs['simulation']['rho_f2_vtk']
    print_geom = inputs['simulation']['print geom']
    print_stl = inputs['simulation']['print stl']

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

    # Restart sim?
    file.write(f'<load_savedstated> {restart_sim} </load_savedstated>\n\n')

    # Write geometry section
    file.write('<geometry>\n')
    # Geometry name
    file.write(f'\t<file_geom> {input_folder}{geom_name}.dat </file_geom>\n')
    # Geometry size
    file.write(f'\t<size> <x> {x_size} </x> <y> {y_size} </y> <z> {z_size} </z> </size>\n')
    # Periodicity
    file.write(f'\t<per>\n')
    file.write(f'\t\t<fluid1> <x> {periodic_x} </x> <y> {periodic_y} </y> <z> {periodic_z} </z> </fluid1>\n')
    file.write(f'\t\t<fluid2> <x> {periodic_x} </x> <y> {periodic_y} </y> <z> {periodic_z} </z> </fluid2>\n')
    file.write(f'\t</per>\n')
    file.write('</geometry>\n\n')

    # Write initial position of fluids
    file.write(f'<init>\n')
    file.write(f'\t<fluid_from_geom> {load_fluid_from_geom} </fluid_from_geom>\n')
    file.write(f'\t<fluid1>\n')
    file.write(f'\t\t <x1> {fluid1_x1} </x1> <y1> {fluid1_y1} </y1> <z1> {fluid1_z1} </z1>\n')
    file.write(f'\t\t <x2> {fluid1_x2} </x2> <y2> {fluid1_y2} </y2> <z2> {fluid1_z2} </z2>\n')
    file.write(f'\t</fluid1>\n')
    file.write(f'\t<fluid2>\n')
    file.write(f'\t\t <x1> {fluid2_x1} </x1> <y1> {fluid2_y1} </y1> <z1> {fluid2_z1} </z1>\n')
    file.write(f'\t\t <x2> {fluid2_x2} </x2> <y2> {fluid2_y2} </y2> <z2> {fluid2_z2} </z2>\n')
    file.write(f'\t</fluid2>\n')
    file.write('</init>\n\n')

    # Write fluid data
    file.write('<fluids>\n')

    file.write(f'\t<Gc> {Gc} </Gc>\n')
    file.write(f'\t<omega_f1> {omega_f1} </omega_f1>\n')
    file.write(f'\t<omega_f2> {omega_f2} </omega_f2>\n')
    file.write(f'\t<force_f1> {force_f1} </force_f1>\n')
    file.write(f'\t<force_f2> {force_f2} </force_f2>\n')
    file.write(f'\t<G_ads_f1_s1> {G_ads_f1_s1} </G_ads_f1_s1>\n')
    file.write(f'\t<G_ads_f1_s2> {G_ads_f1_s2} </G_ads_f1_s2>\n')
    file.write(f'\t<G_ads_f1_s3> {G_ads_f1_s3} </G_ads_f1_s3>\n')
    file.write(f'\t<G_ads_f1_s4> {G_ads_f1_s4} </G_ads_f1_s4>\n')

    file.write(f'\t<rho_f1> {rho_f1} </rho_f1>\n')
    file.write(f'\t<rho_f2> {rho_f2} </rho_f2>\n')

    file.write(f'\t<pressure_bc> {pressure_bc} </pressure_bc>\n')
    file.write(f'\t<rho_f1_i> {rho_f1} </rho_f1_i>\n')
    file.write(f'\t<rho_f2_i> {rho_f2} </rho_f2_i>\n')
    file.write(f'\t<num_pc_steps> {num_pc_steps} </num_pc_steps>\n')
    file.write(f'\t<min_radius> {minimum_radius} </min_radius>\n')
    file.write(f'\t<rho_d> 0.06 </rho_d>\n')

    file.write('</fluids>\n\n')

    # Write output section
    file.write('<output>\n')

    file.write(f'\t<out_folder> {output_folder} </out_folder>\n')
    file.write(f'\t<save_it> {save_iter} </save_it>\n')
    file.write(f'\t<save_sim> {save_sim} </save_sim>\n')
    file.write(f'\t<convergence> {convergence} </convergence>\n')
    file.write(f'\t<it_max> {max_iter} </it_max>\n')
    file.write(f'\t<it_conv> {convergence_iter} </it_conv>\n')
    file.write(f'\t<it_gif> {gif_iter} </it_gif>\n')
    file.write(f'\t<rho_vtk> {rho_f2_vtk} </rho_vtk>\n')
    file.write(f'\t<it_vtk> {vtk_iter} </it_vtk>\n')
    file.write(f'\t<print_geom> {print_geom} </print_geom>\n')
    file.write(f'\t<print_stl> {print_stl} </print_stl>\n')

    file.write('</output>')

    file.close()

    return


def create_relperm_input_file(inputs, input_file_name):

    geom_name = inputs['domain']['geom name']
    nx = inputs['domain']['domain size']['nx']
    ny = inputs['domain']['domain size']['ny']
    nz = inputs['domain']['domain size']['nz']
    num_layers = 2
    domain_size = [nx+(2*num_layers), ny, nz]
    per_x = inputs['domain']['periodic boundary']['x']
    per_y = inputs['domain']['periodic boundary']['y']
    per_z = inputs['domain']['periodic boundary']['z']
    periodic = [per_x, per_y, per_z]
    io_folders = [inputs['input output']['input folder'], inputs['input output']['output folder']]
    # perm_model = inputs['simulation']['perm model']
    num_geoms_or_sims = inputs['rel perm']['num_geoms']
    pressure = inputs['rel perm']['pressure']
    max_iter = inputs['rel perm']['max iterations']
    convergence = inputs['rel perm']['convergence']
    save_vtks = inputs['rel perm']['save vtks']

    # Parse geometry inputs
    x_size = domain_size[0]
    y_size = domain_size[1]
    z_size = domain_size[2]
    periodic_x = periodic[0]
    periodic_y = periodic[1]
    periodic_z = periodic[2]

    # Parse i/o inputs
    input_folder = io_folders[0]
    output_folder = io_folders[1] + '4relperm/'

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
    # # Permeability model
    # file.write(f'\t<perm_model> {perm_model} </perm_model>\n')
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