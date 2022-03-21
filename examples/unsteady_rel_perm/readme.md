# Two-Phase Simulation Python Workflow (a.k.a. laptop-pore)
This is an example of an unsteady state relative permeability simulation. The main outputs will be capillary pressure and relative permeability curves. The number of points (number of pressure steps, under the ```num pressure steps``` key) on the curves and the highest pressure value (```minimum radius``` key) can be set in the input.yml file. On 30 cores, this example took about 45 minutes. This example can also be run on a laptop since the geometry is small. If you allocate 4 cores, it should take on the order of a few hours.

Run the following line in your terminal:\
```python 2_phase_sim.py```

You can run ```python 2_phase_viz.py``` followed by ```python create_gif_and_mp4.py``` to visualize outputs and create an animation.

<img src=/illustrations/laptop_pore_lbm_animation.gif width="500" height="400">

Non-wetting fluid shown in orange traveling through a [texturally equilibrated pore](https://www.digitalrocksportal.org/projects/65). Wetting fluid is shown in blue.

----------------------------------------------------------------------------

<img src=/illustrations/laptop_pore_pc_and_relperm_curve.png width="300" height="500">
