This file contains the sequence of codes to follow from start to end for calculation of Pc curves, percolation pathway and relative permeability of all grain packs/image slices.

1) Create the .dat file geometry for Palabos. 
a) If the geometry is a matrix in MATLAB, use createLBM_MAT2DAT.m
b) If the gemoetry is an image sequence, use createDAT.m
Both codes add 2 blank slices at beginning and end of geometry in YZ plane and a unit mesh slice 3rd from end.

2) Use control_sc2.m to create a string of inputs for the 2-phase LBM simulation.

3) Make/Run 2-phase LBM simulation file ShanChen_v4_new.cpp in bash using input string generated to simulation capillary drainage. 

4) Read the generated vtk files using read_save_fluids.m It will calculate the wetting saturation for all vtk files, will convert the fluid configurations (1 and 2) to .dat files for 1-phase LBM simulation and will find the vtk file where breakthrough occurs and the percolation path/tortuosity. You can choose if you want to generate fluid geometries or calculate precolation path at breakthrough

5) Make/Run 1-phase LBM simulation file permeability_rel.cpp in bash to calculate absolute and relative permeabilities of the flow from fluid geometries.