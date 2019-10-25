## Explanation of inputs for 2 phase flow XML file (input_sc.xml)

<img src="https://github.com/je-santos/MultiphasePorousMediaPalabos/blob/master/illustrations/LBM%20geometry%203D.png" align="middle" width="600" height="400" alt="geometry inputs explanation">
 
Figure 1- Geometry setup for 2-phase flow simulations


**The inputs are explained in the same sequence as the xml file**

**<load_savedstated>:** This input gives the user the option of loading a previous incomplete simulation (True), or starting a new simulation (False). 

**<geometry>:** There are multiple inputs required under this heading.

**<file_geom>:** The name of the geometry for simulation. PALABOS requires the input geometry be in .DAT file format which can be created using the [pre-processing steps](https://github.com/je-santos/MultiphasePorousMediaPalabos/tree/master/pre-processing) The geometry file should be placed in the same folder as the [2-phase simulation code](https://github.com/je-santos/MultiphasePorousMediaPalabos/tree/master/src/2-phase_LBM) or if placed elsewhere, the path should be modified to point to the geometry.

**<size>:** This input requires the size (in voxels) in the X, Y, and Z directions of the geometry (Figure 1).
 
*Note:* 1) PALABOS conducts simulations in X-direction. 
       2) If the blank slices and mesh are added in the pre-proccesing steps, 5 voxels will be added to the original geometry in the X-direction.
	   
	   
	   
	   




