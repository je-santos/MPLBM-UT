================================================================================
Simulation using 2D image
================================================================================

In this example, a 2D image will be used for simulation. The 2D simulation is bascily the same as 3D simulation, but with only 1 layer in Z direction.

Additionally, we show how to download and simulate on volumes hosted at the digitalrockportal.

.. contents::


################################################################################
Instructions
################################################################################

1. Fetch the raw image from the digitalrockportal:

.. code-block:: bash

  fetch_image.sh

or manually by downloading it from:

https://www.digitalrocksportal.org/projects/79/origin_data/312/

You need to convert it from .tif to .raw.
We provide a preconverted image.

2. Run 

.. code-block:: bash
   
   create_geom_4_2phase.m 
   
to create a computationally efficient (and artifact free) domain for simulation

3. Then, run

.. code-block:: bash

    run2-phase.sh
    
    
To perform the simulation with the inputs specified in the .xml file. At several capillary pressure increments. For more information about this procedure, please refer to the validation_YoungLaplace in the examples folder.


The capillary pressure-saturation curve looks like this:

.. figure:: /illustrations/PcSat.PNG
    :align: center
    :alt: alternate text
    :figclass: align-center

    Resulting Pc-Sat curve


################################################################################
References
################################################################################

1. Santos, J. E., Prodanovic, M., Landry, C. J., & Jo, H. (2018, August 9). Determining the Impact of Mineralogy Composition for Multiphase Flow Through Hydraulically Induced Fractures. Unconventional Resources Technology Conference. doi:10.15530/URTEC-2018-2902986



We welcome contributions
