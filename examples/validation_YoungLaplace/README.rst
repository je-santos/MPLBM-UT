================================================================================
Numerical validation using a bundle of tubes
================================================================================

In this example, I show the trade-off between convergence threshold and accuracy (compared to the Young-Laplace, or in this specific case, the Washburn equation)

.. figure:: /illustrations/nw_flow.jpg
    :align: center
    :alt: alternate text
    :figclass: align-center

    Percolating path of a non-wetting fluid (rock and wetting fluid not shown).

.. contents::


################################################################################
Instructions
################################################################################

1. The domain is created using the create_cap_tube_domain.m file. In there a different number of individual tubes with varying radius can be specified.
2. The convergence threshold (mass balance can be changed in the input_tubes.xml file)
3. The simulation can be carried-out by running the run2-pahse.sh file

################################################################################
References
################################################################################

1. Santos, J. E., Prodanovic, M., Landry, C. J., & Jo, H. (2018, August 9). Determining the Impact of Mineralogy Composition for Multiphase Flow Through Hydraulically Induced Fractures. Unconventional Resources Technology Conference. doi:10.15530/URTEC-2018-2902986



We welcome contributions
