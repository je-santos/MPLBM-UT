import glob
import numpy as np


all_dirs = glob.glob('lbm_results/*')

X = []
y = []
for sim_dir in all_dirs:
    sim = glob.glob(sim_dir+'/*.dat')
    if sim:
        print('list')
        y.append(np.loadtxt(sim[0]).reshape(3,258,128,3))
        
        
        # the x is located in input
    