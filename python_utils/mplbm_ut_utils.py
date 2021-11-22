import parse_input_file
import pore_utils
import create_geom_for_palabos
import create_palabos_input_file
import yaml


class MplbmUTUtils():

    def __init__(self, input_file):

        self.parse_input_file(input_file)

        return

    def parse_input_file(self, input_file):

        # Reads inputs from a yaml file or from a Python dictionary
        if isinstance(input_file, str):
            with open(input_file) as f:
                inputs = yaml.full_load(f)
        else:
            inputs = input_file

        # Check for inputs

        #########################
        # Check simulation type #
        #########################
        if 'simulation type' in inputs:

            if inputs['simulation type'] == '1-phase':
                print("Simulation Type: " + inputs['simulation type'])

            elif inputs['simulation type'] == '2-phase':
                raise NotImplementedError("Automated 2-phase workflow is still in progress!")

            else:
                raise ValueError("Please check your simulation type: it should be either '1-phase' or '2-phase'.")

        else:
            raise KeyError("Please add 'simulation type' entry to input file.")

        ####################
        # Check i/o inputs #
        ####################
        if 'input output' in inputs:

            if 'simulation directory' in inputs['input output']:
                print("Simulation Directory: " + inputs['input output']['simulation directory'])
            else:
                raise KeyError("Please add 'simulation directory' entry to 'input output'.")

            if 'input folder' in inputs['input output']:
                print("Input Folder: " + inputs['input output']['input folder'])
            else:
                raise KeyError("Please add 'input folder' entry to 'input output'.")

            if 'output folder' in inputs['input output']:
                print("Output Folder: " + inputs['input output']['output folder'])
            else:
                raise KeyError("Please add 'output folder' entry to 'input output'.")

        else:
            raise KeyError("Please add 'input output' entry to input file.")

        #########################
        # Check geometry inputs #
        #########################
        if 'geometry' in inputs:
            if 'file name' in inputs['geometry']:
                print("Geometry File: " + inputs['geometry']['file name'])
            else:
                raise KeyError("Please add 'file name' entry to 'geometry'.")

            if 'data type' in inputs['geometry']:
                print("Geometry Data Type: " + inputs['geometry']['data type'])
            else:
                raise KeyError("Please add 'data type' entry to 'geometry'.")

            if 'geometry size' in inputs['geometry']:
                if 'Nx' and 'Ny' and 'Nz' in inputs['geometry']['geometry size']:
                    Nx = inputs['geometry']['geometry size']['Nx']
                    Ny = inputs['geometry']['geometry size']['Ny']
                    Nz = inputs['geometry']['geometry size']['Nz']
                    print("Geometry Size: Nx=" + str(Nx) + ", Ny=" + str(Ny) + ", Nz=" + str(Nz))
                else:
                    raise KeyError("Plase make sure Nx, Ny, and Nz are in the 'geometry size' entry")
            else:
                raise KeyError("Please add 'geometry size' entry to 'geometry'.")

        else:
            raise KeyError("Please add 'geometry' entry to input file.")

        #######################
        # Check domain inputs #
        #######################
        if 'domain' in inputs:

            if 'geom name' in inputs['domain']:
                print("Geometry Name: " + inputs['domain']['geom name'])
            else:
                inputs['domain']['geom name'] = inputs['geometry']['file name'].split('.')[0]
                print("Geometry Name: " + inputs['domain']['geom name'])

            if 'domain size' in inputs['domain']:
                if 'nx' and 'ny' and 'nz' in inputs['domain']['domain size']:
                    nx = inputs['domain']['domain size']['nx']
                    ny = inputs['domain']['domain size']['ny']
                    nz = inputs['domain']['domain size']['nz']
                    print("Domain Size: nx=" + str(nx) + ", ny=" + str(ny) + ", nz=" + str(nz))
                else:
                    raise KeyError("Please make sure nx, ny, and nz are in the 'domain size' entry")

            if 'periodic boundary' in inputs['domain']:
                if 'x' and 'y' and 'z' in inputs['domain']['periodic boundary']:
                    pass
                else:
                    raise KeyError("Please make sure x, y, and z are in the 'periodic boundary' entry")
            else:
                raise KeyError("Please make sure 'periodic boundary' is in the 'domain' entry")

            if 'inlet and outlet layers' in inputs['domain']:
                pass
            else:
                inputs['domain']['inlet and outlet layers'] = 1

            if 'add mesh' in inputs['domain']:
                pass
            else:
                inputs['domain']['add mesh'] = False

            if 'swap xz' in inputs['domain']:
                pass
            else:
                inputs['domain']['swap xz'] = False

            if 'double geom resolution' in inputs['domain']:
                pass
            else:
                inputs['domain']['double geom resolution'] = False

        else:
            raise KeyError("Please add 'domain' entry to input file.")

        ###########################
        # Check simulation inputs #
        ###########################
        if 'simulation' in inputs:

            if inputs['simulation type'] == '1-phase':
                if 'num geoms' in inputs['simulation']:
                    print("Number of Geometries: " + str(inputs['simulation']['num geoms']))
                else:
                    raise KeyError("Please add 'num geoms' to 'simulation' entry.")

                if 'pressure' in inputs['simulation']:
                    print("Applied Pressure: " + str(inputs['simulation']['pressure']))
                else:
                    raise KeyError("Please add 'pressure' to 'simulation' entry.")

                if 'max iterations' in inputs['simulation']:
                    print("Maximum Iterations: " + str(inputs['simulation']['max iterations']))
                else:
                    raise KeyError("Please add 'max iterations' to 'simulation' entry.")

                if 'convergence' in inputs['simulation']:
                    print("Convergence: " + str(inputs['simulation']['convergence']))
                else:
                    raise KeyError("Please add 'convergence' to 'simulation' entry.")

                if 'save vtks' in inputs['simulation']:
                    print("Saving Velocity and Medium vtks: " + str(inputs['simulation']['save vtks']))
                else:
                    raise KeyError("Please add 'save vtks' to 'simulation' entry.")

            elif inputs['simulation type'] == '2-phase':
                raise NotImplementedError('2-phase coming soon!')

        else:
            raise KeyError("Please add 'simulation' entry to input file.")

        ##############################
        # Check visualization inputs #
        ##############################
        if 'visualization' in inputs:
            pass
        else:
            pass

        ###############################
        # Store all inputs as objects #
        ###############################
        # Simulation type
        self.simulation_type = inputs['simulation type']

        # Input output
        self.sim_dir = inputs['input output']['simulation directory']
        self.input_dir = inputs['input output']['input folder']
        self.output_dir = inputs['input output']['output folder']
        self.io_folders = [self.input_dir, self.output_dir]

        # Geometry
        self.geom_file_name = inputs['geometry']['file name']
        self.data_type = inputs['geometry']['data type']
        self.geom_file = self.sim_dir + '/' + self.input_dir + self.geom_file_name
        self.Nx = inputs['geometry']['geometry size']['Nx']
        self.Ny = inputs['geometry']['geometry size']['Ny']
        self.Nz = inputs['geometry']['geometry size']['Nz']

        # Domain
        self.geom_name = inputs['domain']['geom name']
        self.nx = inputs['domain']['domain size']['nx']
        self.ny = inputs['domain']['domain size']['ny']
        self.nz = inputs['domain']['domain size']['nz']
        self.domain_size = [self.nx + (2 * self.num_layers), self.ny, self.nz]
        self.per_x = inputs['domain']['periodic boundary']['x']
        self.per_y = inputs['domain']['periodic boundary']['y']
        self.per_z = inputs['domain']['periodic boundary']['z']
        self.periodic = [self.per_x, self.per_y, self.per_z]
        self.num_layers = inputs['domain']['inlet and outlet layers']
        self.add_mesh = inputs['domain']['add mesh']
        self.swap_xz = inputs['domain']['swap xz']
        self.double_geom_res = inputs['domain']['double geom resolution']

        # Simulation
        self.num_procs = inputs['simulation']['num procs']
        self.num_geoms_or_sims = inputs['simulation']['num geoms']
        self.pressure = inputs['simulation']['pressure']
        self.max_iter = inputs['simulation']['max iterations']
        self.convergence = inputs['simulation']['convergence']
        self.save_vtks = inputs['simulation']['save vtks']

        # Visualization
        # Coming soon

        return

