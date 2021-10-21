================================================================================
Pythonic tools
================================================================================

Future python release


################################################################################
Instructions
################################################################################

1. Fetch the raw image from the digitalrockportal:

.. code-block:: bash

    ./get_geom_fromDRP.sh

or manually, from:

https://www.digitalrocksportal.org/projects/57/origin_data/287/

2. Run:

.. code-block:: bash

	python create_geom_4sim.m

to create a computationally efficient (and artifact free) domain for simulation.

3. Run

.. code-block:: bash

    run1-phase.sh


To perform the simulation with the inputs specified in the .xml file

