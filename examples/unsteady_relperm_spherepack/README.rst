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

1. Fetch the raw image from the digitalrockportal:

.. code-block:: bash

  wget "https://www.digitalrocksportal.org/media/projects/47/archive.zip"
  unzip -j archive.zip origin/311/images/spheres_a10_dx0.04_n500_segmented_unsigned_char.raw -d input/.

or manually by downloading it from:

https://www.digitalrocksportal.org/projects/47/origin_data/311/

2. Run create_geom_4sim.m to create a computationally efficient (and artifact free) domain for simulation

3. Run

.. code-block:: bash

    run2-phase.sh
    
 To perform the simulation with the inputs specified in the .xml file



################################################################################
References
################################################################################

1. Santos, J. E., Prodanovic, M., Landry, C. J., & Jo, H. (2018, August 9). Determining the Impact of Mineralogy Composition for Multiphase Flow Through Hydraulically Induced Fractures. Unconventional Resources Technology Conference. doi:10.15530/URTEC-2018-2902986



We welcome contributions
