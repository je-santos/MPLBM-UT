# Steady State Relative Permeability
The idea of this example is to (1) initialize non-wetting fluid in the domain with varying levels of saturation; and then (2) equilibrate the fluid configuration with applied forces and periodic boundary conditions; finally, (3) run a relative permeability simulation using the equilibrated fluid configurations. The example as set up here took about 3 hours to run on 50 cores.

To run the example, run this command in your terminal:\
```python 2_phase_sim.py```

For visualization, ```2_phase_viz.py``` and ```1_phase_viz.py``` are also available.

A small note: the convergence tolerance should be set a bit lower than the recommended 1e-4 for steady state simulations. You may also need to use some trial and error to find the right forces and convergence tolerance. In this example the convergence tolerance is 2e-5, and the force applied to both fluids is 3e-4.

<img src=./finney_pack_lbm_animation.gif width="500" height="400">

Non-wetting fluid shown in orange traveling through a [Finney Pack](https://www.digitalrocksportal.org/projects/47). Wetting fluid is shown in light blue, and grains are light gray.

----------------------------------------------------------------------------

<img src=./relperm_curve.png width="500">
