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

  fetch_image.sh

or manually by downloading it from:

https://www.digitalrocksportal.org/projects/47/origin_data/311/

2. Run 

.. code-block:: bash
   
   create_geom_4_2phase.m 
   
to create a computationally efficient (and artifact free) domain for simulation

3. Then, run

.. code-block:: bash

    run2-phase.sh
    
    
To perform the simulation with the inputs specified in the .xml file. At several capillary pressure increments. For more information about this procedure, please refer to the validation_YoungLaplace in the examples folder.



################################################################################
References
################################################################################

1. Santos, J. E., Prodanovic, M., Landry, C. J., & Jo, H. (2018, August 9). Determining the Impact of Mineralogy Composition for Multiphase Flow Through Hydraulically Induced Fractures. Unconventional Resources Technology Conference. doi:10.15530/URTEC-2018-2902986



We welcome contributions
