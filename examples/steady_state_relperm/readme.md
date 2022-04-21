# Steady State Relative Permeability
The idea of this example is to (1) initialize non-wetting fluid in the domain with varying levels of saturation from the output of a PoreSpy drainage simulation; and then (2) equilibrate the fluid configuration with applied forces and periodic boundary conditions; finally, (3) run a relative permeability simulation using the equilibrated fluid configurations. The example as set up here took about 4.5 hours to run on 60 cores.

To run the example, run this command in your terminal:\
```python 2_phase_sim.py```

For visualization, ```2_phase_viz.py``` and ```1_phase_viz.py``` are also available.

Some notes:
1) Depending on your goals, the convergence tolerance can be set a bit lower than the recommended 1e-4 for steady state simulations. Keep in mind, the lower the tolerance, the longer the simulation.
2) You may also need to use some trial and error to find the right forces. In order to figure out the right forces, watch how the capillary number (Ca) changes at the beginning of your simualtion. Ideally for pore-scale flow, Ca should be ~10^-4 - 10^-6. 
3) For this example, the forces for both fluids are set to 1e-4 and Ca was ~10^-4 for all runs. The convergence tolerance was 1e-4.

<img src=/illustrations/finney_pack_porespy_init_2.png width="350" height="280"> <img src=/illustrations/finney_pack_lbm_animation_2.gif width="350" height="280">

<img src=/illustrations/finney_pack_porespy_init_3.png width="350" height="280"> <img src=/illustrations/finney_pack_lbm_animation_3.gif width="350" height="280">

Non-wetting fluid shown in orange traveling through a [Finney Pack](https://www.digitalrocksportal.org/projects/47). Wetting fluid is shown in light blue, and grains are light gray.

----------------------------------------------------------------------------

<img src=/illustrations/finney_pack_relperm_curve.png width="500">
