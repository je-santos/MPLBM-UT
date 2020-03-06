import numpy as np
import matplotlib.pyplot as plt

pc = []
r  = []

dRho = 0.05
dPc  = dRho/3


for pc_inc in range(1,21):
    
    if pc == []:     
        pc.append(dPc)
    else:
        pc.append(pc[-1]+dPc)
    
    r.append( 0.825/(pc[-1]*3) )

        
conv_inc = {}
conv_inc['1e0']  = [3,3,3,9,9,11,16]
conv_inc['1e-1'] = [2,3,3,5,5,8,14]
conv_inc['1e-2'] = [2,2,3,4,4,5,7]
conv_inc['1e-3'] = [2,2,3,4,4,5,6]