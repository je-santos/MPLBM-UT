## Explanation of inputs for 2 phase flow XML file (input_sc.xml)

<img src="https://github.com/je-santos/MultiphasePorousMediaPalabos/blob/master/illustrations/LBM%20geometry%203D.png" align="middle" width="600" height="400" alt="geometry inputs explanation">
 
Figure 1- Geometry setup for 2-phase flow simulations


**The inputs are explained in the same sequence as the xml file**

**load_savedstated** This input gives the user the option of loading a previous incomplete simulation (True), or starting a new simulation (False). 

**geometry** There are multiple inputs required under this heading.

**file_geom** The name of the geometry for simulation. PALABOS requires the input geometry be in .DAT file format which can be created using the [pre-processing steps](https://github.com/je-santos/MultiphasePorousMediaPalabos/tree/master/pre-processing) The geometry file should be placed in the same folder as the [2-phase simulation code](https://github.com/je-santos/MultiphasePorousMediaPalabos/tree/master/src/2-phase_LBM) or if placed elsewhere, the path should be modified to point to the geometry.

**size** This input requires the size (in voxels) in the X, Y, and Z directions of the geometry (Figure 1).
 
*Note:* 

1) PALABOS conducts simulations in X-direction, so please double-check X and Z directions of the geometry.
2) If the blank slices and mesh are added in the pre-proccesing step, 5 voxels will be added to the original geometry in the X-direction.
	   
**init** There are multiple inputs required under this heading.

**fluid1** This input requires the initial positions of fluid 1 (usually invading fluid). As shown in Figure 1, these inputs are x1 to x2, y1 to y2, and z1 to z2 given in orange color. Fluid 1 has one edge as the YZ plane at x=x1 at inlet side of the geometry and the second edge at x=x2 (which is also the left edge of fluid interface at t=0). As the entire fluid in the Y-Z space is filled with fluid 1 between x=x1 and x=x2, the Y and Z direction limits are the geometry limits.

**fluid2** This input requires the initial positions of fluid 2 (usually defending fluid). As shown in Figure 1, these inputs are x1 to x2, y1 to y2, and z1 to z2 given in blue color. Fluid 2 has one edge as the YZ plane at x=x1 (which is the right edge of fluid interface at t=0) and the second edge at x=x2 at outlet side of the geometry. As the entire fluid in the Y-Z space is filled with fluid 1 between x=x1 and x=x2, the Y and Z direction limits are again the geometry limits.
	   
**fluids** There are multiple inputs required under this heading.

**Gc**
**omega_f1**
**omega_f2**
**force_f1**
**force_f2**
**Wetting forces**
(G_ads_f1_s1, G_ads_f1_s2, G_ads_f1_s3, G_ads_f1_s4)

**rho_f1** This input takes the intial density of fluid 1 throughout the fluid
**rho_f2** This input takes the intial density of fluid 2 throughout the fluid
**pressure_bc** This input asks if a pressure gradient will be applied in the geometry. True for flow simulations and False for equilibrium calculations.
**num_pc** 
**rho_f1_i**
**rho_f2_i**
**rho_f2_f**
**rho_d**
**drho_f2**

**output** There are multiple inputs required under this heading.
**out_folder**
**save_sim**
**convergence**
**it_max**
**it_conv**
**it_info**
**it_gif**
**rho_vtk**
**it_vtk**
**print_geom**
**print_stl**



