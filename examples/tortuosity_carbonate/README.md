================================================================================
Tortuosity calculation in a carbonate rock
================================================================================
In this example, the method to calculate the tortuosity of the flow paths in a carbonate is presented.



################################################################################
Instructions
################################################################################
1. Fetch the raw image from the digitalrockportal:

.. code-block:: bash

    wget "https://www.digitalrocksportal.org/media/projects/57/archive.zip" &&
    unzip -j archive.zip 57/origin/287/images/grid_25.bin -d input/. &&
    rm archive.zip

or manually, from:

https://www.digitalrocksportal.org/projects/57/origin_data/287/

2. Run create_geom_4sim.m to create a computationally efficient (and artifact free) domain for simulation.

3. Run

.. code-block:: bash

    run1-phase.sh


To perform the simulation with the inputs specified in the .xml file
