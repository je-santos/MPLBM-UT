================================================================================
Multiphase LBM Toolbox: Permeable media analysis using the Palabos library and in-house codes
================================================================================
.. image:: https://zenodo.org/badge/174389090.svg
   :target: https://zenodo.org/badge/latestdoi/174389090
   

This workflow calculates capillary pressure curves, relative permeability, tortuosity, contact angles, and the percolation pathway of grain packs/image slices. The direct fluid flow simulation is performed using Palabos v2.

This repository was created by Javier E. Santos and Abhishek Bihani in collaboration with Christopher Landry, Hugh Daigle, Masa Prodanovic, Wenhui Song, and Michael Pyrcz

################################################################################
Illustrations
################################################################################

.. figure:: /illustrations/purple_lbm.gif
    :align: right
    :alt: alternate text
    :figclass: align-right

    Non-wetting fluid (purple) traveling through a spherepack.

----------------------------------------------------------------------------

.. figure:: /illustrations/Tropical_simulation.gif
    :align: right
    :alt: alternate text
    :figclass: align-right

    Unsteady flow simulation showing a capillary pressure increment with different colors.

----------------------------------------------------------------------------


.. figure:: /illustrations/percolation.png
    :align: right
    :alt: alternate text
    :figclass: align-right

    Percolating path of a non-wetting fluid (rock and wetting fluid not shown).

----------------------------------------------------------------------------

.. figure:: /illustrations/single_phase.png
    :align: right
    :alt: alternate text
    :figclass: align-right
   

    Velocity magnitude of a single phase flow simulation.


################################################################################
Overview
################################################################################

This library provides tools to interact with the palabos simulation engine. Its aim is to provide ready-to-go tools to compute relevant properties for permeable media. It includes scripts for creating efficient geometries for simulation, single and multiphase models with various boundary conditions, and extensive postprocesing tools. The Palabos backend makes the simulation engine flexible enought to be run on a single core of a personal computer or in hundreds of nodes of a supercomputer.

################################################################################
Requirements
################################################################################

- Matlab or Octave (Python alternative coming soon)

- Unix system (the Windows bash or the terminal in Mac)

- GCC 7.2
- OpenMPI 2.1.1 or MPICH2
- Python 2.7.x (to complie Palabos)


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



A) Pre-processing (Matlab/Octave):

- To create the geometry for simulating with Palabos (.dat file),

    a) If the geometry is a 3D array, use pre-processing/create_geom_edist.m
    b) If the geometry is an image sequence, use pre-processing/createDAT.m

B) Two-Phase LBM Simulation (cpp w/MPI using PALABOS):

    - Update geometry and simulation parameters in input xml file. An explanation of every input (in english and chinese) is provided in   `examples/1_two_phase_template_explanation <examples/1_two_phase_template_explanation/readme.md>`__

C) Post-processing (Matlab/Octave):

    - Read the generated  files using post-porcessing/domains_4_kr.m
    (It will calculate the wetting saturation for all vtk files, will convert the fluid configurations (1 and 2) to .dat files for 1-phase LBM simulation and will find the vtk file where breakthrough occurs and the percolation path/tortuosity. You can choose if you want to generate fluid geometries or calculate percolation path at breakthrough)

D) Single-Phase LBM Simulation (cpp w/MPI using PALABOS):

    - Update geometry and simulation parameters in input xml file

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
- Steady state: an video example is shown in:  `youtube video <https://www.youtube.com/watch?v=wc8ZxwejcHk>`__

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

Note that this process takes a few hours.

-----------------------------------------------------------------------------------------------------------

I am getting a Java Heap Memory error in Matlab?
You need to change the JavaHeapMemory setting in Matlab:

If you're working on a remote system/cluster or supercomputer, the easiest way to do this is to find and change your matlab.prf file diretly.

You will find the matlab.prf in your user storage directory. It will be something like:
```[user home]/.matlab/[Matlab version]/matlab.prf```

Once the file is open add the following line to the file:
```JavaMemHeapMax = [Java Memory Code]```
You will need to find the sepcific Java memory code that corresponds to the amount of RAM you need.

Or you can try this fix (If you are on a remote system, cluster, or supercomputer this may not work): https://www.mathworks.com/matlabcentral/answers/74296-change-java-heap-memory-settings-without-starting-matlab

################################################################################
Citing the workflow
################################################################################

If you use our workflow, please cite as:

Santos J., Bihani A., Landry C., Multi-Phase Porous Media for Palabos. Git code (2019). 
10.5281/zenodo.3842279


################################################################################
Author's Publications
################################################################################
1. Bihani A., Daigle H., Santos J., Landry C., Prodanović M., Milliken K. (2019). H44B-06: Insight into the Sealing Capacity of Mudrocks determined using a Digital Rock Physics Workflow. AGU Fall Meeting, 9-13 December, San Francisco, USA.

2. Santos, J. E., Prodanovic, M., Landry, C. J., & Jo, H. (2018, August 9). Determining the Impact of Mineralogy Composition for Multiphase Flow Through Hydraulically Induced Fractures. Unconventional Resources Technology Conference. doi:10.15530/URTEC-2018-2902986

3. Landry, C. J., Karpyn, Z. T., and Ayala, O. (2014), Relative permeability of homogenous‐wet and mixed‐wet porous media as determined by pore‐scale lattice Boltzmann modeling, Water Resour. Res., 50, 3672– 3689, doi:10.1002/2013WR015148.

4. Santos, J. E., Xu, D., Jo, H., Landry, C. J., Prodanović, M., & Pyrcz, M. J. (2020). PoreFlow-Net: A 3D convolutional neural network to predict fluid flow through porous media. Advances in Water Resources, 138, 103539.


################################################################################
Bibliographic references
################################################################################

1. Unsteady state simulation set-up: Pan, C., Hilpert, M., and Miller, C. T. ( 2004), Lattice‐Boltzmann simulation of two‐phase flow in porous media, Water Resour. Res., 40, W01501, https://doi.org/10.1029/2003WR002120.

2. Contact angle approximation: Huang, H., Thorne, D. T., Schaap, M. G., & Sukop, M. C. (2007). Proposed approximation for contact angles in Shan-and-Chen-type multicomponent multiphase lattice Boltzmann models. Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 76(6), 1–6. https://doi.org/10.1103/PhysRevE.76.066701.

3. Fluid flow model: Santos, J. E. (2018). Lattice-Boltzmann modeling of multiphase flow through rough heterogeneously wet fractures. University of Texas at Austin (Master thesis). https://repositories.lib.utexas.edu/handle/2152/69246.


################################################################################
Code references
################################################################################

1. Absolute permeability code modified from Degruyter et al. - http://www.palabos.org/documentation/tutorial/permeability.html
2. Skeletonization code modified from Kollmannsberger Philip- https://www.mathworks.com/matlabcentral/profile/authors/4556277-philip-kollmannsberger
3. Dijkstra's algorithm code modified from Kirk Joseph- https://www.mathworks.com/matlabcentral/fileexchange/12850-dijkstra-s-shortest-path-algorithm



We welcome contributions
