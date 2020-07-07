================================================================================
Relative permeability calculation using an unsteady state approach
================================================================================

In this example, an unsteady state (pressure difference between inlet and outlet) method will be used to simulate flow through a permeable material (in this case, a spherepack).
Additionally, we show how to download and simulate on volumes hosted at the digitalrockportal.

.. figure:: /illustrations/nw_flow.jpg
    :align: center
    :alt: alternate text
    :figclass: align-center

    Percolating path of a non-wetting fluid (rock and wetting fluid not shown).

.. contents::


################################################################################
Instructions
################################################################################

1. *Get geometry:* Fetch the raw image from the digitalrockportal:

.. code-block:: bash

  fetch_image.sh

or manually by downloading it from here_ 

.. _here: https://www.digitalrocksportal.org/projects/47/origin_data/311/

2. *Pre-processing:* Run 

.. code-block:: bash
   
   create_geom_4_2phase.m 
   
to create a computationally efficient (and artifact free) domain for simulation

3. *Two-Phase simulation:* Then, run

.. code-block:: bash

    run2-phase.sh
 
    
To perform the simulation with the inputs specified in the `2-phase .xml-file <examples/unsteady_relperm_spherepack/input_spherepack.xml>`_ at several capillary pressure increments. For more information about this procedure, please refer to the `validation_YoungLaplace example  <examples/validation_YoungLaplace/>`_.

    
4. *Post-processing:* Then, run 
 
.. code-block:: bash
 
   create_geoms_4_kr.m
   
to calculate saturations for the pressure-saturation curve and create geometries for running the relative permeability simulations. The Pc-Saturation curve can then be plot with the generated saturations and the pressure inputs for the simulation.    

5. *One-Phase simulation:* Finally, run

.. code-block:: bash

    run1-phase.sh

to run the simulation for generating the relative permeability values with the inputs specified in the `1-phase .xml-file <examples/unsteady_relperm_spherepack/input_rel_perm.xml>`_.

The capillary Pc-Saturation curve looks like this:

.. figure:: /illustrations/PcSat.PNG
    :align: center
    :alt: alternate text
    :figclass: align-center

    Resulting Pc-Saturation curve
    
################################################################################
Additional Options
################################################################################

- *Handling different input file types*

If you want to run simulations on your projects with different input file types like-

+--------+--------------------------------------------+
| Type 1 | Raw image file (as shown above)            |
+--------+--------------------------------------------+
| Type 2 | MATLAB (.MAT) file                         |
+--------+--------------------------------------------+
| Type 3 | image slices like .tiff/.png/.jpg          |
+--------+--------------------------------------------+
| Type 4 | spherical grain centre coordinates & radii |
+--------+--------------------------------------------+

First, add your initial geometry to the input folder 
Then, update your specific file type in the script and run


.. code-block:: bash

    create_geom_4_2phase_more_file_types.m

-----

- *Running the toolbox on supercomputing resources (eg. TACC)*

1. Copy the toolbox and the created domain for simulation to your workspace on the supercomputer

2. Update details in the example script batchjob_4TACC.sh and make (compile) the (2-phase and 1-phase) codes.

3. Start the LBM simulation by running

.. code-block:: bash

    sbatch batchjob_4TACC.sh
        
Notes-
    
 a) Make sure to test smaller cases first before running a full-fledged simulation on the supercomputer
 b) Observe best practices of making/running jobs from `appropriate nodes`_ 
      
.. _appropriate nodes: https://portal.tacc.utexas.edu/user-guides/stampede2#citizenship-loginnodes

-----


- *Percolation path calculation and visualization*

You can measure the length and visualize the fluid path at first breakthrough after finishing the simulation. To do this, update details and run 

.. code-block:: bash

    percolation_path.m


This can create 3D visualizations like the one shown below. 

.. figure:: /illustrations/percolation.png
    :align: right
    :alt: alternate text
    :figclass: align-right

    Percolating path of a non-wetting fluid (rock and wetting fluid not shown).

################################################################################
References
################################################################################

1. Santos, J. E., Prodanovic, M., Landry, C. J., & Jo, H. (2018, August 9). Determining the Impact of Mineralogy Composition for Multiphase Flow Through Hydraulically Induced Fractures. Unconventional Resources Technology Conference. doi:10.15530/URTEC-2018-2902986

2. Bihani A., Daigle H., Santos J., Landry C., Prodanović M., Milliken K. (2019). Insight into the Sealing Capacity of Mudrocks determined using Digital Rock Physics. TACC Symposium for Texas Researchers (TACCSTER), 26-27 September, Austin, TX, USA. doi:10.26153/tsw/6874

We welcome contributions
