# Young-Laplace Validation
<img src=./capillary_tube_viz.png width="600"> 
This example has three different purposes:

1. It is an additional validaiton of MPLBM-UT
2. It provides a tolerance test
3. It provides some performance data. 


To run the test, enter the following line in your terminal:\
```python 2_phase_sim.py```

You can visualize the outputs using ```2_phase_viz.py```.

## Validation
Tolerances of 1e-4 and 1e-5 overlap each other and agree very well with the Young-Laplace equation (the difference comes down to spatial discritization error i.e. tubes/cylinders in discrete space are not perfect circles). This means that for accurate and consistent results we should choose a tolerance of 1e-4 or below.

<img src=./young_laplace_validation.png width="500"> 

## Performance
As tolerance decreases, simulation time rapidly increases. We recommend using 1e-4 as a tolerance since it provides a good combination between speed and accuracy. 

<img src=./young_laplace_validation_performance.png width="500"> 
