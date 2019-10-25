================================================================================
Multiphase LBM for Porous Media using the Palabos library
================================================================================

This workflow can be used for calculation of capillary pressure curves, relative permeability and the percolation pathway of grain packs/image slices

Workflow created by Abhishek Bihani and Javier Santos

Collaborators: Christopher Landry, Hugh Daigle and Masa Prodanovic

.. image:: https://github.com/je-santos/MultiphasePorousMediaPalabos/blob/master/illustrations/nw_flow.jpg
:Percolating path of a non-wetting fluid (rock and wetting fluid not shown).

.. contents::


################################################################################
Overview
################################################################################

The Palabos library makes the engine flexible enought to be run on a single core of a personal laptop or in hundreds of nodes in a supercomputer.


################################################################################
Installation
################################################################################

.. code-block:: bash

  ./Install.sh

MATLAB or Octave are used to run the domain construction and post-processing tools.

################################################################################
Running a simulation
################################################################################

A) Pre-processing (MATLAB):

- To create the .dat file geometry for Palabos,
a) If the geometry is a matrix in MATLAB, use createLBM_MAT2DAT.m
b) If the geometry is an image sequence, use createDAT.m

Both codes add 2 blank slices at beginning and end of geometry in YZ plane, (optional) adds a unit mesh slice 3rd from outlet

B) 2-Phase LBM Simulation (PALABOS):

- Update geometry and simulation parameters in input xml file
- Modify the Makefile to point out to the directory where Palabos is saved and make
- Run 2-phase LBM simulation file ShanChen.cpp in bash to simulate capillary drainage

C) Post-processing (MATLAB):

- Read the generated vtk files using read_save_fluids.m
(It will calculate the wetting saturation for all vtk files, will convert the fluid configurations (1 and 2) to .dat files for 1-phase LBM simulation and will find the vtk file where breakthrough occurs and the percolation path/tortuosity. You can choose if you want to generate fluid geometries or calculate percolation path at breakthrough)

D) 1-Phase LBM Simulation (PALABOS):

- Update geometry and simulation parameters in input xml file
- Modify the Makefile to point out to the directory where Palabos is saved and make
- Make/Run 1-phase LBM simulation file permeability.cpp in bash to calculate absolute and relative permeabilities of the flow from fluid geometries

Optional: Run porethroat_dist.m to calculate pore and throat size distribution of geometry from vtk file created during 2-phase LBM simulation

################################################################################
Numerical validations
################################################################################
1

################################################################################
Examples
################################################################################

1

################################################################################
FAQ
################################################################################

- Segmentation errors
- Image is loaded properly but it doesn't look right: This is likely because Palabos engine switches the X and Z coordinates
- SCons compilation: create a conda enviroment with python2:

.. code-block:: bash
  conda create --name py2 python=2.7
  ./Install.sh


and install from there

If you use our workflow, please cite as: Santos J., Bihani A., Landry C., Daigle H., and Prodanovic M. Multi-Phase Porous Media for Palabos. Git code (2019). https://github.com/je-santos/MultiphasePorousMediaPalabos

################################################################################
Publications
################################################################################

1 Urtec
2 AWR
3 NeuralNets
4 a couple under review

################################################################################
Credits/References
################################################################################

1. Geometry creation and permeability code modified from Degruyter et al. - http://www.palabos.org/documentation/tutorial/permeability.html
2. Skeletonization code modified from Kollmannsberger Philip- https://www.mathworks.com/matlabcentral/profile/authors/4556277-philip-kollmannsberger
3. Dijkstra's algorithm code modified from Kirk Joseph- https://www.mathworks.com/matlabcentral/fileexchange/12850-dijkstra-s-shortest-path-algorithm


We welcome contributions
