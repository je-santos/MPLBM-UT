# Two-Phase Micromodel of a Texturally Equilbrated Porous Medium (TEPM)
This example includes utilities that can be used to create a micromodel from a binary image. An unsteady state simulation is setup here, but it can be modified for steady state or single phase as well. This example took about 4.5 hours to run on 40 cores.

To start the simulation, run this in your terminal:\
```python 2_phase_sim.py```

To create the animation, run the following two commands:\
```python 2_phase_micromodel_viz.py```\
```python create_gif_and_mp4.py```
  
----------------------------------------------------------------------------

<img src=/illustrations/rg_theta30_phi30_lbm_animation.gif width="500" height="400">

Non-wetting fluid shown in orange traveling through a [TEPM](https://www.digitalrocksportal.org/projects/65) micromodel. Wetting fluid is shown in blue.
