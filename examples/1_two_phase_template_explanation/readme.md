## Explanation of inputs for 2 phase flow XML file (input_sc.xml)

<img src="https://github.com/je-santos/MultiphasePorousMediaPalabos/blob/master/illustrations/LBM%20geometry%203D.png" align="middle" width="600" height="400" alt="geometry inputs explanation">

Figure 1- Geometry setup for 2-phase flow simulations


**The inputs are explained in the same sequence as the xml file**

**load_savedstated:** This input gives the user the option of loading a previous incomplete simulation (True), or starting a new simulation (False).

**geometry:** There are multiple inputs required under this heading.

**file_geom:** This input asks for the name of the geometry for simulation. PALABOS requires the input geometry be in .dat file format which can be created using the [pre-processing steps](https://github.com/je-santos/MultiphasePorousMediaPalabos/tree/master/pre-processing) The geometry file should be placed in the same folder as the [2-phase simulation code](https://github.com/je-santos/MultiphasePorousMediaPalabos/tree/master/src/2-phase_LBM) or if placed elsewhere, the path should be modified to point to the geometry.

**size:** This input requires the size (in voxels) in the X, Y, and Z directions of the geometry (Figure 1).

*Note:*

1) PALABOS conducts simulations in X-direction, so please double-check X and Z directions of the geometry.
2) If the blank slices and a mesh are added in the pre-proccesing step, the original size in the X-direction would be larger.

**init:** There are multiple inputs required under this heading.

**fluid1:** This input requires the initial positions of fluid 1 (usually invading fluid). As shown in Figure 1, these inputs are x1 to x2, y1 to y2, and z1 to z2 given in orange color. Fluid 1 has one edge as the YZ plane at x=x1 at inlet side of the geometry and the second edge at x=x2 (which is also the left edge of fluid interface at t=0). As the entire fluid in the Y-Z space is filled with fluid 1 between x=x1 and x=x2, the Y and Z direction limits are the geometry limits.

**fluid2:** This input requires the initial positions of fluid 2 (usually defending fluid). As shown in Figure 1, these inputs are x1 to x2, y1 to y2, and z1 to z2 given in blue color. Fluid 2 has one edge as the YZ plane at x=x1 (which is the right edge of fluid interface at t=0) and the second edge at x=x2 at outlet side of the geometry. As the entire fluid in the Y-Z space is filled with fluid 1 between x=x1 and x=x2, the Y and Z direction limits are again the geometry limits.

**fluids:** There are multiple inputs required under this heading.

**Gc:** Interparticle (cohesion) force. This input controls the fluid-fluid interfacial tension. This value assures phase separation. A stable separation is reached by Gc > 1/(rho_f1+rho_f2). A value of 0.9 is recommended.

**omega_f1 and omega_f2:** This parameter is used to calculate the kinematic viscosity of fluid 1 and 2, using:  v=( 1 / omega_fi - 0.5 ) / c^2.

**force_f1 and force_f2** If this term is different than zero, a driving force will be added to each fluid (i.e. gravitational) in the x-direction. The pressure boundary conditions are suggested to be turned off and periodicity should be enabled (to reach a steady-state).

**Wetting forces:**
(G_ads_f1_s1, G_ads_f1_s2, G_ads_f1_s3, G_ads_f1_s4): These terms refer to the interaction force between the fluids and the solid walls. This code has the option to add 4 different wetting conditions ( 4 different solid surfaces ), but more could be added with ease. In the 3D image, the voxels labeled with 1, 3, 5, 6 are assigned G_ads_f1_s1, G_ads_f1_s2, G_ads_f1_s3, G_ads_f1_s4, respectively (2 is reserved for inside solids, 4 for the neutral-wet mesh and 0 for the fluids). The contact angle is calculated as: cos(theta) = 4*G_ads_f1_si/( Gc*( rho_f1-rho_d ) )

<img src="https://latex.codecogs.com/svg.latex?\Large&space;cos(\theta)=\frac{4G_{ads_{f1,si}}-b\pm\sqrt{b^2-4ac}}{G_c(rho_{f1}-rho_d)}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

**rho_f1:** This input takes the initial density of fluid 1 throughout the geometry.

**rho_f2:** This input takes the initial density of fluid 2 throughout the geometry.

*Note:* A stable value for both densities is 2. High density ratios tend to be numerically unstable.

**pressure_bc:** This input asks if a pressure gradient will be applied in the geometry. Please use True for un-steady state flow simulations and False for steady state calculations.

**rho_f1_i:** This input takes the initial density of fluid 1 at the inlet pressure boundary and is kept constant.

**rho_f2_i:** This input takes the initial density of fluid 2 at the outlet pressure boundary at the beginning of the simulation.

**rho_f2_f:** This input takes the final density of fluid 2 at the outlet pressure boundary at the end of the simulation. The difference between the inlet and the outlet pressure boundaries decides the capillary pressure.

**rho_d:** This input takes the dissolved density of one phase in the other (both fluid1 and fluid2). The default value may be kept 0.06.

**drho_f2:** This input takes the decrement in the pressure of fluid 2 at the outlet pressure boundary at every step (capillary pressure change) in the simulation. A range of 0.01 to 0.1 may be input depending on balance between sensitivity / computational time, as smaller decrement will require a longer time but will have greater sensitivity to measure change in fluid movement.

The pressure in the Shan-Chen model is calculated as:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;P(x)=\frac{\rho_1(x)+\rho_2(x)}{3}+G_c\frac{\rho_1(x)\rho_2(x)}{3}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />


**output:** There are multiple inputs required under this heading.

**out_folder:** This input takes the name of the folder where all the output files will be stored.

**save_sim:** This input asks if user wants to save the simulation lattices after every capillary pressure decrement for both fluid 1 and fluid 2. The saved files are large (> 1 GB) but are overwritten after every decrement and may be used to restart the simulation from that step.

**convergence:** This input takes the value of the convergence criterion for the simulation. The convergence value is inversely proportional to the computational time and the accuracy. See the examples for best practices.

**it_max:** This input takes the value of maximum iterations allowed at a particular capillary pressure if the convergence is not reached.

**it_conv:** This input takes the value of number of iterations after which to check if convergence criterion is satisfied. This takes a certain computational overhead, since general statistics have to be computed.

**it_gif:** This input takes the value of number of iterations after which 2D images of the geometry crossection (.gif) showing the density and current fluid configuration are to be saved.

**rho_vtk:** This input asks if 3D geometry files (.vtk) for both fluids 1 and 2 are to be saved: True, or only for for Fluid 1: False.

**it_vtk:** This input takes the  number of iterations after which the 3D geometry files (.vtk) showing the density and current fluid configuration are to be saved. If this number is greater than 100000 no vtk files will be output (which keeps the output lighter and the run faster).

**print_geom:** This input asks if a 3D geometry file (.vtk) is to be saved at the beginning of the simulation.

**print_stl:** This input asks if a 3D geometry file (.stl) is to be saved at the beginning of the simulation. This file may be used for 3D printing the geometry or for viewing.

*Note:* The .vtk files may be viewed in [PARAVIEW](https://www.paraview.org/).
