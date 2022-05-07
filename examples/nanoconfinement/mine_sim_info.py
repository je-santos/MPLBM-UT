import numpy as np
from hdf5storage import loadmat

import csv
import glob
import os
from shutil import copyfile

import matplotlib.pyplot as plt    
    


def full_dict(folder='domains_real/', results={}):
    #results = {}
    #folder = 'domains_real/'
    o_files = [os.path.basename(x) for x in glob.glob(folder + '/*.o*')]
    
    for file in o_files:
        sim_name = file.split('.',1)[0]
        try:
            results[sim_name]
            results[sim_name]['num_sims'] += 1
        except KeyError:
            results[sim_name] = {}
            results[sim_name]['num_sims'] = 1 
            
        num_sim = results[sim_name]['num_sims']
        results[sim_name][f'convergence_{num_sim}'] = False
    
        text = open(folder + file)
        conv = None
        for line in text:
            if line.find('Velocity convergence')>-1:
                conv = float(line.split()[2])
            if line.find('simulation has reached velocity convergence')>-1:
                results[sim_name][f'convergence_{num_sim}'] = True
                count = 0
            if results[sim_name][f'convergence_{num_sim}'] == True:
                count += 1
                if count == 3: 
                    it = [int(i) for i in line.split() if i.isdigit()][0]
                    results[sim_name][f'it_{num_sim}'] = it
                    #it = line
                if count == 4: #sim time
                    #print(line)
                    simtime = float(line.split()[3])/3600
                    results[sim_name][f'time_{num_sim}'] = simtime 
        results[sim_name][f'conv_{num_sim}'] = conv
    return results
    

def clean_dict(folders = ['domains_real/']):
    
    results = {}
    #[full_dict(data, results) for data in ['domains_real/','domains_syn/']]
    [full_dict(data, results) for data in folders]
    
    sorted_results = {}
    for sim in results.keys():
        num_sims = results[sim]['num_sims']
        convergence = False
        for num in range(1,num_sims+1):
            if results[sim][f'convergence_{num}']:
                convergence = True
        if convergence:
            max_it, min_eps, max_time = 0,100,0
            for num in range(1,num_sims+1):
                if results[sim][f'convergence_{num}']:
                    max_it =  np.max([min_eps,results[sim][f'it_{num}']])
                    min_eps = np.min([min_eps,results[sim][f'conv_{num}']])
                    max_time = np.max([min_eps,results[sim][f'time_{num}']])
            sorted_results[sim]={'max_it':max_it, 'min_eps': min_eps, 
                                 'max_time':max_time}
    return sorted_results
            
            

def print_resend(results):
    for sim in results.keys():
        num_sims = results[sim]['num_sims']
        convergence = False
        for num in range(1,num_sims+1):
            if results[sim][f'convergence_{num}']:
                convergence = True
        if not convergence:
            best_conv = np.Inf
            max_its   = 0
            for num in range(1,num_sims+1):
                try:
                    results[sim][f'it_{num}']
                except KeyError:
                    results[sim][f'it_{num}'] = 0
                if not results[sim][f'conv_{num}']:
                    pass
                else:
                    best_conv = np.min([best_conv, results[sim][f'conv_{num}']])
                    max_its   = np.max([max_its,   results[sim][f'it_{num}']])
            print(f'{sim} at {max_its} its achieved c = {best_conv}')
            

def create_log():
    folder = ['domains_syn/']
    sorted_res = clean_dict(folders = folder)
    
    for sim in sorted_res.keys():
    #for sim in b[-50:]:
        print_dict = {}
        
        print_dict['pressure [MPa]'] = int(sim.split('_')[-1])
        print_dict['convergence'] = sorted_res[sim]['min_eps']
        print_dict['iterations'] = sorted_res[sim]['max_it']
        print_dict['num_of_cores'] = 144
        print_dict['simulation_time'] = sorted_res[sim]['max_time']
        
        try:
            vz = loadmat(f'results_{folder[0].split("_")[1]}{sim}.mat')['uz']
        except FileNotFoundError:
            print(f'Sim {sim} not found')
            continue
        except ValueError:
            print(f'Sim {sim} has an error')
            continue
            
        
        mv = np.mean(vz)
        perm = mv*(1.9119e-10)*(1/0.5e-9)**2  # to Darcy*(1/9.869233e-13)
        print_dict['apparent permeability (lus)'] = perm
        
        with open(f'results_{folder[0].split("_")[1]}{sim}_log.csv', 'w', newline='') as csv_file:  
            writer = csv.writer(csv_file)
            for key, value in print_dict.items():
               writer.writerow([key, value])
        
        
        
        if print_dict['pressure [MPa]'] == 20:
            print_dict.pop('apparent permeability (lus)')
            print_dict.pop('pressure [MPa]')
            print_dict['permeability (lus)'] = perm*(1e-6)**2
            
            new_name = '_'.join( sim.split("_")[:-1] )
            with open(f'results_{folder[0].split("_")[1]}{new_name}_log.csv', 'w', newline='') as csv_file:  
                writer = csv.writer(csv_file)
                for key, value in print_dict.items():
                   writer.writerow([key, value])
        
        
        

def create_folders_from_feats():
    
    # change loc to hard drive
    os.chdir(r'F:\nanolbm_bu')
    folder = ['features_syn/','features_real/']
    feats = [glob.glob(f + '*.mat') for f in folder]
    feats = [item for sublist in feats for item in sublist]
    feat_files = [os.path.basename(file) for file in feats]
    feat_names = [f.split('.')[0] for f in feat_files]
    
    sx = [st.replace('features','results') for st in feats]
    
    sims={
    'sims0'   : sx,
    'sims1'  : [st.replace('.','_1.') for st in sx],
    'sims2'  : [st.replace('.','_2.') for st in sx],
    'sims5'  : [st.replace('.','_5.') for st in sx],
    'sims10' : [st.replace('.','_10.') for st in sx],
    'sims20' : [st.replace('.','_20.') for st in sx],
    }
    
    
    for num in range(len(feats)):
        sim_folder = feat_names[num]
        if feats[num].find('syn') > -1:
            sim_folder = '374_0' + sim_folder
        try:
            os.mkdir(fr'../DRP_Upload/{sim_folder}')
        except FileExistsError:
            pass
        try:
            png = feats[num].split('.')[0]+'.png'
            copyfile(feats[num],fr'../DRP_Upload/{sim_folder}/features.mat')
            copyfile(png,fr'../DRP_Upload/{sim_folder}/features.png')
        except:
            print('help')
         
        for p in [0,1,2,5,10,20]:
            sim_name = f'P_{p}_MPa'
            if p == 0:
                sim_name = 'LBM'
            log_loc = sims[f'sims{p}'][num].replace('.mat','_log.csv')
            if os.path.isfile(log_loc):
                try:
                    copyfile(sims[f'sims{p}'][num],f'../DRP_Upload/{sim_folder}/{sim_name}.mat')
                    copyfile(log_loc,fr'../DRP_Upload/{sim_folder}/{sim_name}.csv')
                except FileNotFoundError:
                    print(f'{feat_names[num]} file not found')
            else:
                print(f'log {log_loc} not found')
            
       


def create_folders():
    # create log files
    
    # check all the feats real and fake
    # create folders with good names
    # copy features with new name and pngs
    # copy sims and log if converged with new name
    pass
 