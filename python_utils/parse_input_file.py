import yaml


def parse_input_file(input_file):

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
            print("Simulation Type: " + inputs['simulation type'])

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

            if 'restart sim' in inputs['simulation']:
                print("Restart Simulation: " + str(inputs['simulation']['restart sim']))
            else:
                raise KeyError("Please add 'restart sim' to 'simulation' entry.")

            if 'rho_f1' in inputs['simulation']:
                print("Fluid 1 Equilibrium Density: " + str(inputs['simulation']['rho_f1']))
            else:
                raise KeyError("Please add 'rho_f1' to 'simulation' entry.")

            if 'rho_f2' in inputs['simulation']:
                print("Fluid 2 Equilibrium Density: " + str(inputs['simulation']['rho_f2']))
            else:
                raise KeyError("Please add 'rho_f2' to 'simulation' entry.")

            if 'force_f1' in inputs['simulation']:
                print("Fluid 1 Force: " + str(inputs['simulation']['force_f1']))
            else:
                raise KeyError("Please add 'force_f1' to 'simulation' entry.")

            if 'force_f2' in inputs['simulation']:
                print("Fluid 2 Force: " + str(inputs['simulation']['force_f2']))
            else:
                raise KeyError("Please add 'force_f2' to 'simulation' entry.")

            if 'pressure bc' in inputs['simulation']:
                print("Using Pressure Boundary Conditions: " + str(inputs['simulation']['pressure bc']))
            else:
                raise KeyError("Please add 'force_f2' to 'simulation' entry.")

            if inputs['simulation']['pressure bc'] == True:

                if 'minimum radius' in inputs['simulation']:
                    print("Minimum radius of invasion (voxels): " + str(inputs['simulation']['minimum radius']))
                else:
                    raise KeyError("Please add 'minimum radius' to 'simulation' entry if using pressure boundary conditions.")

                if 'num pressure steps' in inputs['simulation']:
                    print("Number of capillary pressure steps: " + str(inputs['simulation']['num pressure steps']))
                else:
                    raise KeyError("Please add 'num pressure steps' to 'simulation' entry if using pressure boundary conditions.")

            if 'fluid init' in inputs['simulation']:
                if inputs['simulation']['fluid init'] == 'drainage':
                    print("Initial fluid configuration setup for drainage simulation.")
                    inputs['simulation']['fluid 1 init']['x1'] = 1
                    inputs['simulation']['fluid 1 init']['x2'] = 2
                    inputs['simulation']['fluid 1 init']['y1'] = 1
                    inputs['simulation']['fluid 1 init']['y2'] = ny
                    inputs['simulation']['fluid 1 init']['z1'] = 1
                    inputs['simulation']['fluid 1 init']['z2'] = nz

                    inputs['simulation']['fluid 2 init']['x1'] = 3
                    inputs['simulation']['fluid 2 init']['x2'] = nx + 2*inputs['domain']['inlet and outlet layers']
                    inputs['simulation']['fluid 2 init']['y1'] = 1
                    inputs['simulation']['fluid 2 init']['y2'] = ny
                    inputs['simulation']['fluid 2 init']['z1'] = 1
                    inputs['simulation']['fluid 2 init']['z2'] = nz

                elif inputs['simulation']['fluid init'] == 'custom':
                    print("Initial fluid configuration setup with custom settings.")
                else:
                    raise KeyError("Please set 'fluid init' in 'simulation' entry to 'drainage' or 'custom'.")
            else:
                raise KeyError("Please add 'fluid init' to 'simulation' entry.")

            if 'Gc' in inputs['simulation']['fluid data']:
                pass
            else:
                raise KeyError("Please add 'Gc' to 'fluid data' under 'simulation' entry.")

            if 'omega_f1' in inputs['simulation']['fluid data']:
                pass
            else:
                raise KeyError("Please add 'omega_f1' to 'fluid data' under 'simulation' entry.")

            if 'omega_f2' in inputs['simulation']['fluid data']:
                pass
            else:
                raise KeyError("Please add 'omega_f2' to 'fluid data' under 'simulation' entry.")

            if 'G_ads_f1_s1' in inputs['simulation']['fluid data']:
                pass
            else:
                raise KeyError("Please add 'G_ads_f1_s1' to 'fluid data' under 'simulation' entry.")

            if 'G_ads_f1_s2' in inputs['simulation']['fluid data']:
                pass
            else:
                raise KeyError("Please add 'G_ads_f1_s2' to 'fluid data' under 'simulation' entry.")

            if 'G_ads_f1_s3' in inputs['simulation']['fluid data']:
                pass
            else:
                raise KeyError("Please add 'G_ads_f1_s3' to 'fluid data' under 'simulation' entry.")

            if 'G_ads_f1_s4' in inputs['simulation']['fluid data']:
                pass
            else:
                raise KeyError("Please add 'G_ads_f1_s4' to 'fluid data' under 'simulation' entry.")

            if 'convergence' in inputs['simulation']:
                pass
            else:
                raise KeyError("Please add 'convergence' to 'simulation' entry.")

            if 'convergence iter' in inputs['simulation']:
                pass
            else:
                raise KeyError("Please add 'convergence iter' to 'simulation' entry.")

            if 'max iterations' in inputs['simulation']:
                pass
            else:
                raise KeyError("Please add 'max iterations' to 'simulation' entry.")

            if 'save sim' in inputs['simulation']:
                print("Save simulation state: " + str(inputs['simulation']['save sim']))
            else:
                raise KeyError("Please add 'max iterations' to 'simulation' entry.")

            if 'save iter' in inputs['simulation']:
                if inputs['simulation']['save sim']==True:
                    print("Saving simulation state every " + str(inputs['simulation']['save iter']) + " iterations.")
            else:
                raise KeyError("Please add 'save iter' to 'simulation' entry if saving simulation state.")

            if 'gif iter' in inputs['simulation']:
                pass
            else:
                raise KeyError("Please add 'gif iter' to 'simulation' entry.")

            if 'vtk iter' in inputs['simulation']:
                print("Saving density vtks every " + str(inputs['simulation']['vtk iter']) + " iterations.")
            else:
                raise KeyError("Please add 'vtk iter' to 'simulation' entry.")

            if 'rho_f2_vtk' in inputs['simulation']:
                print("Saving density fluid 1 and fluid 2 vtks: " + str(inputs['simulation']['rho_f2_vtk']))
            else:
                raise KeyError("Please add 'rho_f2_vtk' to 'simulation' entry.")

            if 'print geom' in inputs['simulation']:
                pass
            else:
                inputs['simulation']['print geom'] = True

            if 'print stl' in inputs['simulation']:
                pass
            else:
                inputs['simulation']['print stl'] = False

    else:
        raise KeyError("Please add 'simulation' entry to input file.")

    ##############################
    # Check visualization inputs #
    ##############################
    if 'visualization' in inputs:
        pass
    else:
        pass

    return inputs

