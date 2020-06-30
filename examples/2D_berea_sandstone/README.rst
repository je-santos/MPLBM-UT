================================================================================
Simulation using 2D image
================================================================================
Author: JIA Haowei, South West Petroleum University
E-mail contact: nv4dll@outlook.com


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


