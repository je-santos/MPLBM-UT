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


tube_radii = [20.5,15.5,10,5.5,4.5,3.5,2.5]
    
conv_inc = {}
conv_inc['1e0']  = [3,3,3,9,9,11,16]
conv_inc['1e-1'] = [2,3,3,5,5,8,14]
conv_inc['1e-2'] = [2,2,3,4,4,5,7]
conv_inc['1e-3'] = [2,2,3,4,4,5,6]

conv_inc['Y-L'] = [2,2,3,4,5,6,8]

plt.plot([tube_radii[i] for i in range(len(conv_inc['1e0']))],
            [pc[i] for i in conv_inc['Y-L']])
plt.plot([tube_radii[i] for i in range(len(conv_inc['1e0']))],
            [pc[i] for i in conv_inc['1e0']],'.--', markersize=10)
plt.plot([tube_radii[i] for i in range(len(conv_inc['1e0']))],
            [pc[i] for i in conv_inc['1e-1']],'.--', markersize=10)
plt.plot([tube_radii[i] for i in range(len(conv_inc['1e0']))],
            [pc[i] for i in conv_inc['1e-2']],'.--', markersize=10)
plt.plot([tube_radii[i] for i in range(len(conv_inc['1e0']))],
            [pc[i] for i in conv_inc['1e-3']],'.--', markersize=10)

plt.legend(['Young-Laplace','1e0','1e-1','1e-2','1e-3'])

plt.xlabel('Approximate tube radius [voxels]')
plt.ylabel('Intrusion P$_c$ [lattice units]')
plt.title('Effect of convergence criterion on the solution')