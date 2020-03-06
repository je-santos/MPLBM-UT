================================================================================
Numerical validation using a bundle of tubes
================================================================================

In this example, I show the trade-off between convergence threshold and accuracy (compared to the Young-Laplace solution for a capillary cube with a circular cross-section).

.. figure:: /illustrations/tubes_nw.gif
    :align: center
    :alt: alternate text
    :figclass: align-center

    Drainage process in a bundle of tubes

.. contents::


################################################################################
Instructions
################################################################################

1. The domain is created using the create_cap_tube_domain.m file. In there a different number of individual tubes with varying radius can be specified.
2. The convergence threshold (mass balance can be changed in the input_tubes.xml file)
3. The simulation can be carried-out by running the run2-phase.sh file

################################################################################
Results
################################################################################

.. figure:: /illustrations/cap_tubes.png
    :align: center
    :alt: alternate text
    :figclass: align-center

    The simulations show that an epsilon of 1e-2 provides an accurate solution. The curves diverge a negligible ammount when the tube gets really small (it gets increasingly difficult to define what is the actual radious of a cross-section with a few voxels in length).


################################################################################
References
################################################################################

1. Santos, J. E., Prodanovic, M., Landry, C. J., & Jo, H. (2018, August 9). Determining the Impact of Mineralogy Composition for Multiphase Flow Through Hydraulically Induced Fractures. Unconventional Resources Technology Conference. doi:10.15530/URTEC-2018-2902986



We welcome contributions
