================================================================================
MPLBM-UT: MultiPhase and singlephase LBM library for permeable media analysis 
================================================================================
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5203295.svg
   :target: https://doi.org/10.5281/zenodo.5203295
   
MPLBM-UT supports the calculation of capillary pressure and relative permeability curves and single phase permeability of 3D binary images. The direct fluid flow simulation is performed using Palabos v2.2.1 and the pre-processing and post-processsing are done in python. We use the Shan-Chen model for the multiphase simulation and the BGK and MRT for the singlephase simulation.

We have tested all the above on a laptop and a supercomputer. If you have any issues, drop by our `FAQ section <https://github.com/je-santos/MPLBM-UT/blob/master/README.rst#faq>`_ . If you cannot find what you are looking for, open an issue.


################################################################################
Illustrations
################################################################################

.. figure:: /illustrations/purple_lbm.gif
    :align: right
    :alt: alternate text
    :figclass: align-right

    Steady-state simulation: Non-wetting fluid (purple) traveling through a spherepack.

----------------------------------------------------------------------------

.. figure:: /illustrations/Tropical_simulation.gif
    :align: right
    :alt: alternate text
    :figclass: align-right

    Unsteady-state flow simulation: Showing capillary pressure increments with different colors.

----------------------------------------------------------------------------


.. figure:: /illustrations/single_phase.png
    :align: right
    :alt: alternate text
    :figclass: align-right
   

    Single-phase simulation: Velocity magnitude in log scale.


################################################################################
Overview
################################################################################

This library provides tools to interact with the powerful palabos simulation engine. Its aim is to provide ready-to-go tools to compute relevant properties for permeable media. It includes scripts for creating efficient geometries for simulation, single and multiphase models with various boundary conditions, and extensive postprocesing tools. The Palabos backend makes the simulation engine flexible enought to be run on a single core of a personal computer or in hundreds of nodes of a supercomputer.

################################################################################
Requirements
################################################################################
- Unix system (`Windows Linux Subsystem <https://docs.microsoft.com/en-us/windows/wsl/>`__ or the terminal in Mac also work)
   - We also use wget in the examples to download 3D images from Digital Rocks Portal. Please make sure it is installed on your system if you are running the examples.
- Python 3.6+
   - The following external modules are required: `Numpy <https://numpy.org/>`__, `PyYAML <https://pypi.org/project/PyYAML/>`__, `Vedo <https://vedo.embl.es/>`__, `PyVista <https://docs.pyvista.org/>`__, `Matplotlib <https://matplotlib.org/>`__, and `Scikit-Image <https://scikit-image.org/>`__. These modules will automatically be installed when running the install script, as they are dependencies for the `mplbm_utils package </src/python/mplbm_utils>`__ included in the repo.
- Modern C++ compiler
   - gcc 7.X, gcc 9.4.0
- MPI
   - OpenMPI 2.1.1, MPICH2
Note: The Palabos documentation is not specific on what C++ compilers and MPI work to run and compile the code. They recommend any "modern" compiler and have tested gcc, Intel, and Portland Group. That being said, as long as you have a recent version of MPI and a recent C++ compiler, it should work. Above are some of the C++ compilers and MPI versions that have worked for us. If you would like the change between different compilers, here are some instructions on `how to install and switch between different gcc and g++ versions <https://linuxconfig.org/how-to-switch-between-multiple-gcc-and-g-compiler-versions-on-ubuntu-20-04-lts-focal-fossa>`_. If you are on a cluster/supercomputer and need to change compilers, you may need to ask your system admins about installation.

################################################################################
Some Prerequisites
################################################################################
Some familiarity with the terminal, unix operating systems, and bash will be very useful. Here are a few things to keep in mind when running simulations:

- You can check the number of CPUs available on your system with the command :code:`nproc`. It is highly recommended to check how many processors are used in the examples (Python examples use the keyword :code:`num procs` in the .yml input files, and check any of the .sh files for the :code:`-np` keyword) before running them becuase it may not be ideal or compatable with your system. If you are running on a supercomputer/cluster, please check your system's documentation for how processors are allocated to users.
- In relation to the above and for those new to LBM, please note that LBM simulations can take quite some time to run. The main factors in simulaiton time are the number of cores used, domain size, and convergence tolerance. In general:

   - The more cores, the faster the simulation (till a point when too much communication between cores slows things down a bit); 
   - The larger to domain, the longer the simulation; 
   - The smaller the convergence tolerance, the longer the simulation. 
  It is a good idea to do some system performance testing (ie find the optimum number of cores for your domain size). Generally on CPUs, each processor will optimally handle a 20^3-50^3 section of a domain, but you will need to adjust accordingly to your system. Please see the `Young-Laplace validation example <examples/young_laplace_validation>`__ to get an idea of how tolerance affects run time. In that example, we used 40 cores on a 50x52x175 domain (~76^3) which gives each processor a ~22^3 section to process.
  
- Having a dedicated `Python virtual environment <https://docs.python.org/3/library/venv.html>`__ or `Anaconda <https://www.anaconda.com/>`__ environment is recommended to avoid compatability issues.
  - A note on running Python: Depending on your system configuration, you may need to use ``python3`` in order to run simulations instead of just ``python`` as shown in the examples. 
- Seeing "nan" appear in simulation or terminal output indicates that something is not quite right. Please check that your simulation parameters are correct. If you still can not resolve the problem, please submit an issue and we can do our best to help or fix the bug!

################################################################################
Installation
################################################################################
The installation script will install/compile `Palabos <src/>`_, `the one-phase simulation code <src/1-phase_LBM>`_, `the two-phase Shan-Chen code <src/2-phase_LBM>`_, and the `mplbm_utils python package </src/python/mplbm_utils>`__. In the parent directory of the repo, run the installation script with the following command:

.. code-block:: bash

  ./Install.sh

################################################################################
Running a simulation
################################################################################
The `examples <examples>`__ are a great place to see how the code can be used in different ways. In general, each example has a Python script, input file, an input folder, and an output folder(s):

- The Python script contains everything needed to run the example; this means that running a simulation is as easy as running ``python 2_phase_sim.py`` or ``python 1_phase_sim.py`` in the terminal. 
- We use YAML format for our input files. Please refer to the `readme <examples/readme.md>`__ in the examples folder to see a description of all the inputs.
- The input folder is where simulation geometries are stored. The output folder(s) are there to store simulation results.  

Please refer to the `unsteady state example  </examples/unsteady_rel_perm>`__ for a complete end-to-end workflow.

An general overview of the main steps of a simulation is given below. These processes are automated by functions in the `mplbm_utils </src/python/mplbm_utils>`__ folder.

----------------------------------------------------------------------------

**1) Parsing inputs**

The input.yml files in each example contain all the input options for a simulation. `parse_input_file.py </src/python/mplbm_utils/parse_input_file.py>`__ parses the input file and stores the entries as a Python dicationary.

**2) Pre-processing**

This is necessary to create an efficient geometry for simulating with Palabos (.dat file) from the initial geometry file. `create_geom_for_palabos.py </src/python/mplbm_utils/create_geom_for_palabos.py>`__ uses the utilities found in `pore_utils.py </src/python/mplbm_utils/pore_utils.py>`__ to create the .dat file.  

**3) Run a simulation**

The examples provide either a python file (``2_phase_sim.py`` or ``1_phase_sim.py``) that provides the details of running a simualtion. Based on the user inputs, `create_palabos_input_file.py </src/python/mplbm_utils/create_palabos_input_file.py>`__ creates an XML file compatible with palabos, and then a bash file is created that contains all the necessary information to run either `ShanChen <src/2-phase_LBM/>`__ for 2-phase or one of the `1-phase permeability options <src/1-phase_LBM>`_.

A relative permeability simulation is also possible after a 2-phase simulation. `create_geom_for_rel_perm.py </src/python/mplbm_utils/create_geom_for_rel_perm.py>`__ processes the resulting 2-phase geometries so realtive permeability can be calculated based on individual single phase simulations.

**4) Post-processing**

`parse_palabos_output.py </src/python/mplbm_utils/parse_palabos_output.py>`__ contains the functions necessary to parse and save palabos outputs as easy-to-use text files. 
`create_geom_for_rel_perm.py </src/python/mplbm_utils/create_geom_for_rel_perm.py>`__ also contains the methods to calculate fluid saturation after a 2-phase simulation. 

**5) Plotting and visualization**

Plotting and visualization can be done using the various utilities provided in the `examples <examples>`__ folder, and the `animation_and_plotting </src/python/animation_and_plotting>`__ folder. The 3D visualization tools create iso-surfaces of the fluid density from the .vti files to visualize fluid interfaces. General plotting utilities are also available to create capillary pressure and realtive permeability curves. You can also view .vti files with `Paraview <https://www.paraview.org/>`_ and perform 2-3-4D visualization of fluid interfaces. 


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

*1. Why am I getting an error like "HYDU_create_process (./utils/launch/launch.c:69): execvp error on file ../../src/2-phase_LBM/ShanChen (No such file or directory)"?*

| A- This is because the simulation code has not compiled correctly. To solve this, first try to again run

.. code-block:: bash

  ./Install.sh
  
If this does not work, you can open the bash terminal from the main folder and type following commands sequentially.  


.. code-block:: bash

   cd src
   unzip palabos-v2.2.1.zip
   cd 2-phase_LBM/build
   cmake ..
   make
   cd ../../1-phase_LBM/build
   cmake ..
   make
   
-------------------------------------------------------------------------------------

*2. I am getting a segmentation error. What to do?* 

| A- Image dimensions are not correct, try switching the dimensions.

-------------------------------------------------------------------------------------

*3. My image is loaded properly but it doesn't look right. What happened?*

| A- This is likely because Palabos engine switches the X and Z coordinates. The inbuilt option to switch X and Z coordinates is available as 

.. code-block:: bash

   geom.swapXZ = true;

in all geometry creation files (`for example line 26 <examples/unsteady_relperm_spherepack/create_geom_4_2phase.m>`_.).

---------------------------------------------------------------------------------------------

*4. I have a SCons compilation error. What to do?*

| A- SCons compilation error: Because of the palabos source code.

  "../MultiphasePorousMediaPalabos-0.1.0/src/palabos-v2.0r0/src/gridRefinement/couplingInterfaceGenerator3D.h" line 145,    "return dataProcessors;" should be "return *dataProcessors;".
  
-----------------------------------------------------------------------------------------------------------

*5. Why am I seeing the same line printed multiple times? / Why is the code so slow?*

| A- Probably MPI is not installed in your system, this could be solved by:

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

*6. I am getting a Java Heap Memory error in Matlab. What to do?*

| A- You need to change the JavaHeapMemory setting in Matlab:

If you're working on a remote system/cluster or supercomputer, the easiest way to do this is to find and change your matlab.prf file diretly.

You will find the matlab.prf in your user storage directory. It will be something like:
```[user home]/.matlab/[Matlab version]/matlab.prf```

Once the file is open add the following line to the file:
```JavaMemHeapMax = [Java Memory Code]```
You will need to find the sepcific Java memory code that corresponds to the amount of RAM you need.

Or you can try this fix (If you are on a remote system, cluster, or supercomputer this may not work): https://www.mathworks.com/matlabcentral/answers/74296-change-java-heap-memory-settings-without-starting-matlab

-------------------------------------------------------------------------------------

*7. I am getting an error related to the class auto_ptr*

| A- This is because you are using a newer cpp compiler. Subsitute the auto_ptr with unique_ptr in the .cpp files


################################################################################
Author's Publications
################################################################################
1. Bihani A., Daigle H., Santos J., Landry C., Prodanović M., Milliken K. (2019). H44B-06: Insight into the Sealing Capacity of Mudrocks determined using a Digital Rock Physics Workflow. AGU Fall Meeting, 9-13 December, San Francisco, USA.

2. Santos, J. E., Prodanovic, M., Landry, C. J., & Jo, H. (2018, August 9). Determining the Impact of Mineralogy Composition for Multiphase Flow Through Hydraulically Induced Fractures. Unconventional Resources Technology Conference. doi:10.15530/URTEC-2018-2902986

3. Landry, C. J., Karpyn, Z. T., and Ayala, O. (2014), Relative permeability of homogenous‐wet and mixed‐wet porous media as determined by pore‐scale lattice Boltzmann modeling, Water Resour. Res., 50, 3672– 3689, doi:10.1002/2013WR015148.

4. Santos, J. E., Xu, D., Jo, H., Landry, C. J., Prodanović, M., & Pyrcz, M. J. (2020). PoreFlow-Net: A 3D convolutional neural network to predict fluid flow through porous media. Advances in Water Resources, 138, 103539.

5. Gigliotti A., Hesse M., Prodanovic M., Pore-Scale Simulation of Two-Phase Melt Percolation During Core Formation in Planetesimals (Mar. 2021). LPSC, https://www.hou.usra.edu/meetings/lpsc2021/pdf/2328.pdf

6. Gigliotti A. (2021, August 13), Two-phase percolation in texturally equilibrated porous media, http://dx.doi.org/10.26153/tsw/21533

7. Santos, J. E., Pyrcz, M. J., & Prodanović, M. (2022). 3D Dataset of binary images: A collection of synthetically created digital rock images of complex media. Data in Brief, 107797.


################################################################################
External Publications
################################################################################

1.  Xu, R. et al (2020). Pore-scale study of water adsorption and subsequent methane transport in clay in the presence of wettability heterogeneity. Water Resources Research

2. Jonathan B. Grunewald et al 2021 J. Electrochem. Soc. 168 024521



################################################################################
Bibliographic References
################################################################################

1. Unsteady state simulation set-up: Pan, C., Hilpert, M., and Miller, C. T. ( 2004), Lattice‐Boltzmann simulation of two‐phase flow in porous media, Water Resour. Res., 40, W01501, https://doi.org/10.1029/2003WR002120.

2. Contact angle approximation: Huang, H., Thorne, D. T., Schaap, M. G., & Sukop, M. C. (2007). Proposed approximation for contact angles in Shan-and-Chen-type multicomponent multiphase lattice Boltzmann models. Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 76(6), 1–6. https://doi.org/10.1103/PhysRevE.76.066701.

3. Fluid flow model: Santos, J. E. (2018). Lattice-Boltzmann modeling of multiphase flow through rough heterogeneously wet fractures. University of Texas at Austin (Master thesis). https://repositories.lib.utexas.edu/handle/2152/69246.


################################################################################
Contributing
################################################################################

We welcome contributions. Find some ideas below:

- [ ] Moving boundary problem: proppant transport or formation damage
- [ ] Fluid with variable density
- [ ] 3D grid refinement
- [ ] better initialization for unsteady state sim
- [ ] matlab to python translation: geometry creation
- [ ] Testing `GPU <https://palabos-forum.unige.ch/t/from-cpu-to-gpu-in-80-days-project-complete/3301>`_ capabilities


