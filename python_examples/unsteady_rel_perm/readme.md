# Two-Phase Simulation Python Workflow
All you have to do is the following:
  1) Run ```pwd``` in the terminal once in your simulation directory. Copy and paste the output into ```input.yml``` under the ```simulation directory``` key. 
  2) Run ```run_two_phase_sim.sh``` to start a two-phase Shan-Chen LBM simulation.

More details are to come on all the inputs. For now, there are some comments in ```input.yml``` explaining some of the new things.


.. figure:: lbm_animation.gif
    :align: right
    :alt: alternate text
    :figclass: align-right

    Non-wetting fluid shown in orange traveling through a texturally equilibrated pore. Wetting fluid is shown in blue.

----------------------------------------------------------------------------

.. figure:: pc_and_relperm_curve.png 
    :align: right
    :alt: alternate text
    :figclass: align-right

    Capillary pressure and relative permeability curves.
