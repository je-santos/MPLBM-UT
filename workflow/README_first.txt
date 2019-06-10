This folder contains the sequence of codes to follow from start to end for calculation of Pc curves, percolation pathway and relative permeability of all grain packs/image slices.

Workflow created by Abhishek Bihani and Javier Santos
Collaborators: Christopher Landry, Hugh Daigle and Masa Prodanovic

Pre-processing (MATLAB):

1) To create the .dat file geometry for Palabos, 
a) If the geometry is a matrix in MATLAB, use createLBM_MAT2DAT.m
b) If the geometry is an image sequence, use createDAT.m 
(Credits: Degruyter et al. - http://www.palabos.org/documentation/tutorial/permeability.html)

Both codes add 2 blank slices at beginning and end of geometry in YZ plane, (optional) adds a unit mesh slice 3rd from outlet.

2-Phase LBM Simulation (PALABOS):

2) Update geometry and simulation details in xml file.
3) Make/Run 2-phase LBM simulation file ShanChen_new.cpp in bash using input string generated to simulate capillary drainage. 
(For running in TACC update batch file details and run from work folder)

Post-processing (MATLAB):

4) Read the generated vtk files using read_save_fluids.m.
(It will calculate the wetting saturation for all vtk files, will convert the fluid configurations (1 and 2) to .dat files for 1-phase LBM simulation and will find the vtk file where breakthrough occurs and the percolation path/tortuosity. You can choose if you want to generate fluid geometries or calculate precolation path at breakthrough.)
(Credits: for skeletonization code- Kollmannsberger Philip, for Dijkstra's algorithm code- Kirk Joseph)

1-Phase LBM Simulation (PALABOS):

5) Make/Run 1-phase LBM simulation file permeability_rel.cpp in bash to calculate absolute and relative permeabilities of the flow from fluid geometries.
(Credits: Modified code from Degruyter et al.)

Optional: Run porethroat_dist.m to calculate pore and throat distribution of geometry from vtk file created during 2-phase LBM simulation.