# Two-Phase Simulation Python Workflow (a.k.a. laptop-pore)
This is an example of an unsteady state relative permeability simulation. The main outputs will be capillary pressure and relative permeability curves. The number of points (number of pressure steps, under the ```num pressure steps``` key) on the curves and the highest pressure value (```minimum radius``` key) can be set in the input.yml file. On 30 cores, this example took about 45 minutes. This example can also be run on a laptop since the geometry is small. If you allocate 4 cores, it should take on the order of a few hours.

There are currently two options to run this example:

### Option 1: Using python
Run the following line in your terminal:\
```python 2_phase_sim.py```

You can run ```python 2_phase_viz.py``` followed by ```python create_gif_and_mp4.py``` to visualize outputs and create an animation.

A couple notes:
  1) Depending on your system configuration, you may need to use ```python3``` instead of just ```python```
  2) We test on as many systems as we can, but sometimes Python and bash do not always play well together. If you can't run this using Python, please try option 2.

### Option 2: From the terminal
  1) Run ```pwd``` in the terminal once in your simulation directory. Copy and paste the output into ```input.yml``` under the ```simulation directory``` key. 
  2) Run ```run_two_phase_sim.sh``` to start a two-phase Shan-Chen LBM simulation.
  3) Run ```run_rel_perm.sh``` to start the relative permeability simulation.

<img src=./lbm_animation.gif width="500" height="400">

Non-wetting fluid shown in orange traveling through a [texturally equilibrated pore](https://www.digitalrocksportal.org/projects/65). Wetting fluid is shown in blue.

----------------------------------------------------------------------------

<img src=./pc_and_relperm_curve.png width="300" height="500">
