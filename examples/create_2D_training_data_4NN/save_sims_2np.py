import re
import glob
import numpy as np


all_dirs = glob.glob('lbm_results/*')

X = []
y = []
for sim_dir in all_dirs:
    sim = glob.glob(sim_dir+'/*.dat')
    if sim:
        
        domain_num = re.findall(r'\d+',sim[0])[0]
        domain = glob.glob('input/')
        y.append(np.loadtxt(sim[0]).reshape(258,128,3,3)[1:-1,:,1,0].T) #Y,X,slice,vel
        domain_loc = glob.glob(f'input/{domain_num}_data*')[0]
        X.append(((np.loadtxt(domain_loc).reshape(258,128,3)[1:-1,:,1].T)>0)*1)
        

X_np = np.array(X).astype(np.int8)
y_np = np.array(y).astype(np.float16)


np.save('numpys/X',X_np)
np.save('numpys/y',y_np)