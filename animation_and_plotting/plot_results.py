import numpy as np
import matplotlib.pyplot as plt


# Input saturations
Snw = np.array([0.006194, 0.006363, 0.006530, 0.006679, 0.300658, 0.301136,
                0.431331, 0.500419, 0.508805, 0.513808, 0.514832, 0.514835,
                0.514289])
Sw = 1 - Snw
no_perc_index = np.where(Sw>0.99)[0]
#print(Sw)
Sw[no_perc_index] = 1
#print(Sw)

# Input permeabilities
k = 0.410636;
krw =  np.array([1, 0.547822, 0.537374, 0.531291, 0.527665, 0.0212693, 0.0212781,
                 0.00963508, 0.00735022, 0.00681436, 0.00660228, 0.00677459, 
                 0.00668022])
krnw = np.array([0, 5.92137e-16, 6.05393e-16, 6.31209e-16, 6.505e-16, 6.10992e-16,
                 6.17531e-16, 0.291641, 0.459982, 0.492054, 0.496082, 0.500978, 
                 0.502009])
#print(krw)
krw[no_perc_index] = 1
#print(krw)

# Calculate Capillary Pressure, Pc
rho2_start = 2
rho2_end = 1.4
pressure_step = 0.05
print(len(krw))
d_rho = np.arange(0, (rho2_start-rho2_end), pressure_step)
print(len(d_rho))
print(len(Sw))
print(d_rho)
Pc = d_rho/3  # From MPLBM-UT docs

# Plots

# Rel perm
def plot_rel_perm():
    plt.plot(Sw, krw, 'bo-', linewidth=3, markersize=7)
    plt.plot(Sw, krnw, 'ro-', linewidth=3, markersize=7)
    plt.xlim([0,1])
    plt.grid()
    plt.xlabel(r'$S_{silicate}$', fontsize=label_font_size)
    plt.ylabel(r'$k_r$ [LBM Units]', fontsize=label_font_size)
    plt.title(r'Relative Permeability', fontsize=label_font_size)
    plt.legend([r'$k_{r\ silicate}$', r'$k_{r\ metal-sulphide}$'], fontsize=label_font_size)
    #plt.show()
    
    return


# Pc
def plot_pc():
    plt.plot(Sw, Pc, 'ko-', linewidth=3, markersize=7)
    plt.xlim([0,1])
    plt.grid()
    plt.xlabel(r'$S_{silicate}$', fontsize=label_font_size)
    plt.ylabel(r'$P_c$ [LBM Units]', fontsize=label_font_size)
    plt.title(r'Capillary Pressure', fontsize=title_font_size)
    #plt.show()
    
    return


label_font_size = 18
title_font_size = 18

# Plot Rel Perm
f1 = plt.figure()  # dpi=200, figsize=[6,5]
plot_rel_perm()
f1.savefig('rel_perm_graph.png', dpi=300)

# Plot Pc
f2 = plt.figure()  # dpi=200, figsize=[6,5]
plot_pc()
f2.savefig('pc_graph.png', dpi=300)

# Plot Rel perm and Pc
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
f3 = plt.figure(dpi=120, figsize=[6,10])
# plt.suptitle(r'250$^3$ LBM Real Grains, $\phi = 10\%$, $\theta = 60^o$', fontsize=20)
plt.subplot(2,1,1)
plot_pc()

plt.subplot(2,1,2)
plot_rel_perm()
plt.tight_layout()
f3.savefig('rel_perm_and_pc_graph.png', dpi=300)

# plt.show()
