# Two-Phase Simulation Python Workflow (a.k.a. laptop-pore)
There are currently two options to run this example:

### Option 1: Using python
Run the following line in your terminal:\
```python 2_phase_sim.py```

A couple notes:
  1) Depending on your system configuration, you may need to use ```python3``` instead of just ```python```
  2) We test on as many systems as we can, but sometimes Python and bash do not always play well together. If you can't run this using Python, please try option 2.

### Option 2: From the terminal
  1) Run ```pwd``` in the terminal once in your simulation directory. Copy and paste the output into ```input.yml``` under the ```simulation directory``` key. 
  2) Run ```run_two_phase_sim.sh``` to start a two-phase Shan-Chen LBM simulation.
  3) Run ```run_rel_perm.sh``` to start the relative permeability simulation.

More details are to come on all the inputs. For now, there are some comments in ```input.yml``` explaining some of the new things.

<img src=./lbm_animation.gif width="500" height="400">

Non-wetting fluid shown in orange traveling through a [texturally equilibrated pore](https://www.digitalrocksportal.org/projects/65). Wetting fluid is shown in blue.

----------------------------------------------------------------------------

<img src=./pc_and_relperm_curve.png width="300" height="500">
