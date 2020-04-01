================================================================================
Relative permeability calculation using an steady state approach (constant saturation)
================================================================================

In this example, an steady state (fluids with constant body forces) method will be used to simulate flow through a permeable material (in this case, a spherepack).
Additionally, we show how to download and simulate on volumes hosted at the digitalrockportal.

.. figure:: /illustrations/steady.JPG
    :align: center
    :alt: alternate text
    :figclass: align-center

    Percolating path of a non-wetting fluid 

.. contents::


################################################################################
Instructions
################################################################################

1. Run:

.. code-block:: bash

  create_geom_4_2phase.m

to build the geometry for simulation. In this example, this geometry is a spherepack with geometrically periodic boundary conditions in all directions.


2. In the input file (input_spherepack_S25.xml), I have set a simulation with periodic boundary conditions of a non-wetting fluid with a saturation of 25% (this can be changed in line 22) with a fluid force driving both fluids in the x-direction.

3. Run:

.. code-block:: bash

  run_rel_perm.sh
  
to run the two-phase simulation. When this is done: 

.. code-block:: bash

  create_geoms_4_kr.m
  
builds geometries to calculate the relative permeability of the sample. These simulations are carried out with:

.. code-block:: bash

  run_rel_perm.sh
  
The final relative permeability curve looks like this:

.. figure:: /illustrations/Steady_kr.PNG
    :align: center
    :alt: alternate text
    :figclass: align-center

    Resulting kr curves


################################################################################
References
################################################################################

1. Santos, J. E., Prodanovic, M., Landry, C. J., & Jo, H. (2018, August 9). Determining the Impact of Mineralogy Composition for Multiphase Flow Through Hydraulically Induced Fractures. Unconventional Resources Technology Conference. doi:10.15530/URTEC-2018-2902986



We welcome contributions
