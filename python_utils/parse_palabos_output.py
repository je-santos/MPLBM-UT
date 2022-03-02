import numpy as np


def find_line_in_file(file_name, line_to_match, data_to_add_line_index):

    file = open(file_name)
    data = np.array([])

    for line in file:
        if line_to_match in line:
            line_split = line.split()
            data = np.append(data, float(line_split[data_to_add_line_index]))

    file.close()

    return data


def create_pressure_data_file(inputs):

    sim_dir = inputs['input output']['simulation directory']
    output_dir = inputs['input output']['output folder']

    # Parse pressure data
    pc_data_file = sim_dir + '/' + output_dir + 'output.dat'
    pc_data = find_line_in_file(pc_data_file, 'Pressure difference =', 3)
    # pc_data[0] = 0

    np.savetxt(sim_dir + '/' + output_dir + 'data_Pc.txt', pc_data)

    return


def create_relperm_data_file(inputs):

    sim_dir = inputs['input output']['simulation directory']
    output_dir = inputs['input output']['output folder']
    num_runs = inputs['rel perm']['num_geoms']
    print((num_runs-1)/2)

    # Parse rel perm data
    kr_data_file = sim_dir + '/' + output_dir + r'4relperm/relPerm&vels.txt'
    k_absolute = find_line_in_file(kr_data_file, 'Absolute Permeability', 3)
    kr_data = find_line_in_file(kr_data_file, 'Relative Permeability', 3)
    divide_kr_data = int((num_runs-1)/2) + 1
    krw_data = kr_data[0:divide_kr_data]
    krnw_data = kr_data[divide_kr_data:]
    krnw_data = np.insert(krnw_data, 0, 0)

    np.savetxt(sim_dir + '/' + output_dir + 'data_krw.txt', krw_data)
    np.savetxt(sim_dir + '/' + output_dir + 'data_krnw.txt', krnw_data)
    np.savetxt(sim_dir + '/' + output_dir + 'data_k.txt', k_absolute)

    return


