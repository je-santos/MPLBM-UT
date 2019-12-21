================================================================================
Multiphase LBM Toolbox: Permeable media analysis using the Palabos library and in-house codes
================================================================================

This workflow calculates capillary pressure curves, relative permeability, tortuosity, contact angles, and the percolation pathway of grain packs/image slices

Workflow created by Javier E. Santos and Abhishek Bihani 

Collaborators: Christopher Landry, Hugh Daigle and Masa Prodanovic

.. figure:: /illustrations/percolation.png
    :align: right
    :alt: alternate text
    :figclass: align-right
    
    Percolating path of a non-wetting fluid (rock and wetting fluid not shown).

.. contents::


################################################################################
Overview
################################################################################

The Palabos library makes the engine flexible enought to be run on a single core of a personal laptop or in hundreds of nodes of a supercomputer.

################################################################################
Requirements
################################################################################

- Matlab or Octave (Python alternative coming soon)

- Unix system (the Windows bash or the terminal in Mac)

################################################################################
Installation
################################################################################

.. code-block:: bash

  ./Install.sh

################################################################################
Running a simulation
################################################################################

Please refer to the unsteady state example for a complete workflow


----------------------------------------------------------------------------



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

- Young's equation

- Young-Laplace

- Washburn

################################################################################
Example description
################################################################################

- Unsteady state
- Steady state: an video example is shown in 

(`youtube <https://www.youtube.com/watch?v=wc8ZxwejcHk>`__)

################################################################################
FAQ
################################################################################

Segmentation error: Image dimensions are not correct, try switching the dimensions.

-------------------------------------------------------------------------------------

Image is loaded properly but it doesn't look right: This is likely because Palabos engine switches the X and Z coordinates

---------------------------------------------------------------------------------------------


SCons compilation error: create a conda enviroment with python2 (Palabos needs it):

.. code-block:: bash

  conda create --name py2 python=2.7
  ./Install.sh


-----------------------------------------------------------------------------------------------------------

Why am I seeing the same line printed multiple times? / Why is the code so slow?: Probably MPI is not installed in your system, this could be solved by:

.. code-block:: bash

  sudo apt-get install build-essential
  sudo apt-get install gfortran
  cd /tmp
  wget http://www.mpich.org/static/downloads/1.4.1/mpich2-1.4.1.tar.gz
  tar xzvf mpich2-1.4.1.tar.gz
  cd mpich2-1.4.1/
  ./configure
  make
  sudo make install

################################################################################
Citing the workflow
################################################################################

If you use our workflow, please cite as: 

Santos J., Bihani A., Landry C., Multi-Phase Porous Media for Palabos. Git code (2019). https://github.com/je-santos/MultiphasePorousMediaPalabos


################################################################################
Author's Publications
################################################################################
1. Bihani A., Daigle H., Santos J., Landry C., Prodanović M., Milliken K. (2019). H44B-06: Insight into the Sealing Capacity of Mudrocks determined using a Digital Rock Physics Workflow. AGU Fall Meeting, 9-13 December, San Francisco, USA.

2. Santos, J. E., Prodanovic, M., Landry, C. J., & Jo, H. (2018, August 9). Determining the Impact of Mineralogy Composition for Multiphase Flow Through Hydraulically Induced Fractures. Unconventional Resources Technology Conference. doi:10.15530/URTEC-2018-2902986

3. Landry, C. J., Karpyn, Z. T., and Ayala, O. ( 2014), Relative permeability of homogenous‐wet and mixed‐wet porous media as determined by pore‐scale lattice Boltzmann modeling, Water Resour. Res., 50, 3672– 3689, doi:10.1002/2013WR015148.  

4. NeuralNets

5. AWR2

################################################################################
Bibliographic references
################################################################################

1. Unsteady state simulation set-up: Pan, C., Hilpert, M., and Miller, C. T. ( 2004), Lattice‐Boltzmann simulation of two‐phase flow in porous media, Water Resour. Res., 40, W01501, doi:10.1029/2003WR002120 .

2. Contact angle approximation: Huang, H., Thorne, D. T., Schaap, M. G., & Sukop, M. C. (2007). Proposed approximation for contact angles in Shan-and-Chen-type multicomponent multiphase lattice Boltzmann models. Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 76(6), 1–6. https://doi.org/10.1103/PhysRevE.76.066701 .

3. Fluid flow model: Santos, J. E. (2018). Lattice-Boltzmann modeling of multiphase flow through rough heterogeneously wet fractures. University of Texas at Austin (Master thesis). https://repositories.lib.utexas.edu/handle/2152/69246 .


################################################################################
Code references
################################################################################

1. Geometry creation and permeability code modified from Degruyter et al. - http://www.palabos.org/documentation/tutorial/permeability.html
2. Skeletonization code modified from Kollmannsberger Philip- https://www.mathworks.com/matlabcentral/profile/authors/4556277-philip-kollmannsberger
3. Dijkstra's algorithm code modified from Kirk Joseph- https://www.mathworks.com/matlabcentral/fileexchange/12850-dijkstra-s-shortest-path-algorithm


We welcome contributions
