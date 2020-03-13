# English explanation of inputs for 2 phase flow XML file (input_sc.xml)
## followed by the explaination in Chinese

<img src="https://github.com/je-santos/MultiphasePorousMediaPalabos/blob/master/illustrations/LBM%20geometry%203D.png" align="middle" width="600" height="400" alt="geometry inputs explanation">

Figure 1- Geometry setup for 2-phase flow simulations


**The inputs are explained in the same sequence as the xml file**

**load_savedstated:** This input gives the user the option of loading a previous incomplete simulation (True), or starting a new simulation (False).

**geometry:** There are multiple inputs required under this heading.

**file_geom:** This input asks for the name of the geometry for simulation. PALABOS requires the input geometry be in .dat file format which can be created using the [pre-processing steps](https://github.com/je-santos/MultiphasePorousMediaPalabos/tree/master/pre-processing) The geometry file should be placed in the same folder as the [2-phase simulation code](https://github.com/je-santos/MultiphasePorousMediaPalabos/tree/master/src/2-phase_LBM) or if placed elsewhere, the path should be modified to point to the geometry.

**size:** This input requires the size (in voxels) in the X, Y, and Z directions of the geometry (Figure 1).

*Note:*

1) PALABOS conducts simulations in X-direction, so please double-check X and Z directions of the geometry.
2) If the blank slices and a mesh are added in the pre-proccesing step, the original size in the X-direction would be larger.

**init:** There are multiple inputs required under this heading.

**fluid1:** This input requires the initial positions of fluid 1 (usually invading fluid). As shown in Figure 1, these inputs are x1 to x2, y1 to y2, and z1 to z2 given in orange color. Fluid 1 has one edge as the YZ plane at x=x1 at inlet side of the geometry and the second edge at x=x2 (which is also the left edge of fluid interface at t=0). As the entire fluid in the Y-Z space is filled with fluid 1 between x=x1 and x=x2, the Y and Z direction limits are the geometry limits.

**fluid2:** This input requires the initial positions of fluid 2 (usually defending fluid). As shown in Figure 1, these inputs are x1 to x2, y1 to y2, and z1 to z2 given in blue color. Fluid 2 has one edge as the YZ plane at x=x1 (which is the right edge of fluid interface at t=0) and the second edge at x=x2 at outlet side of the geometry. As the entire fluid in the Y-Z space is filled with fluid 1 between x=x1 and x=x2, the Y and Z direction limits are again the geometry limits.

**fluids:** There are multiple inputs required under this heading.

**Gc:** Interparticle (cohesion) force. This input controls the fluid-fluid interfacial tension. This value assures phase separation. A stable separation is reached by Gc > 1/(rho_f1+rho_f2). A value of 0.9 is recommended.

**omega_f1 and omega_f2:** This parameter is used to calculate the kinematic viscosity of fluid 1 and 2, using:  v=( 1 / omega_fi - 0.5 ) / c^2.

**force_f1 and force_f2** If this term is different than zero, a driving force will be added to each fluid (i.e. gravitational) in the x-direction. The pressure boundary conditions are suggested to be turned off and periodicity should be enabled (to reach a steady-state).

**Wetting forces:**
(G_ads_f1_s1, G_ads_f1_s2, G_ads_f1_s3, G_ads_f1_s4): These terms refer to the interaction force between the fluids and the solid walls. This code has the option to add 4 different wetting conditions ( 4 different solid surfaces ), but more could be added with ease. In the 3D image, the voxels labeled with 1, 3, 5, 6 are assigned G_ads_f1_s1, G_ads_f1_s2, G_ads_f1_s3, G_ads_f1_s4, respectively (2 is reserved for inside solids, 4 for the neutral-wet mesh and 0 for the fluids). The contact angle is calculated as: cos(theta) = 4*G_ads_f1_si/( Gc*( rho_f1-rho_d ) )

<img src="https://latex.codecogs.com/svg.latex?\Large&space;cos(\theta)=\frac{4G_{ads_{f1,si}}}{G_c(rho_{f1}-rho_d)}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

**rho_f1:** This input takes the initial density of fluid 1 throughout the geometry.

**rho_f2:** This input takes the initial density of fluid 2 throughout the geometry.

*Note:* A stable value for both densities is 2. High density ratios tend to be numerically unstable.

**pressure_bc:** This input asks if a pressure gradient will be applied in the geometry. Please use True for un-steady state flow simulations and False for steady state calculations.

**rho_f1_i:** This input takes the initial density of fluid 1 at the inlet pressure boundary and is kept constant.

**rho_f2_i:** This input takes the initial density of fluid 2 at the outlet pressure boundary at the beginning of the simulation.

**rho_f2_f:** This input takes the final density of fluid 2 at the outlet pressure boundary at the end of the simulation. The difference between the inlet and the outlet pressure boundaries decides the capillary pressure.

**rho_d:** This input takes the dissolved density of one phase in the other (both fluid1 and fluid2). The default value may be kept 0.06.

**drho_f2:** This input takes the decrement in the pressure of fluid 2 at the outlet pressure boundary at every step (capillary pressure change) in the simulation. A range of 0.01 to 0.1 may be input depending on balance between sensitivity / computational time, as smaller decrement will require a longer time but will have greater sensitivity to measure change in fluid movement.

The pressure in the Shan-Chen model is calculated as:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;P(x)=\frac{\rho_1(x)+\rho_2(x)}{3}+G_c\frac{\rho_1(x)\rho_2(x)}{3}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

often, the second term can be neglected because its many orders of magnitude smaller than the first.

To calculate the capillary pressure of the system we use:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;P_c=P_{nw}-P_{w}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

if we substitute the expression to calculate the pressure, we get:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;P_c=\frac{2}{3}-\frac{(2-\Delta\rho)}{3}=\frac{\Delta\rho}{3}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

where the first term represents the pressure at the inlet and consequently, the second term is the pressure at the outlet.

Comparing this expression with the Young-Laplace equation for a capillary tube (circular cross section), we finally get:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;\Delta\rho=\frac{6\sigma}{r}cos(\theta)" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

For the example shown above, we performed a buble test where the interfacial tension (sigma) showed a value of 0.15. Substituting that, for a non-wetting condition of G_ads_f1 = -0.4 (156.4 degrees), we get:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;\Delta\rho\approx\frac{0.825}{r}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />



**output:** There are multiple inputs required under this heading.

**out_folder:** This input takes the name of the folder where all the output files will be stored.

**save_sim:** This input asks if user wants to save the simulation lattices after every capillary pressure decrement for both fluid 1 and fluid 2. The saved files are large (> 1 GB) but are overwritten after every decrement and may be used to restart the simulation from that step.

**convergence:** This input takes the value of the convergence criterion for the simulation. The convergence value is inversely proportional to the computational time and the accuracy. See the examples for best practices.

**it_max:** This input takes the value of maximum iterations allowed at  if the convergence is not reached. Roughly, in our supercomputer (LoneStar5) we can get 1M its ruuning in 96 cores for two days for a 200x200x200 domain. This is useful in case one would like to save the simulation dats files before the max run time (in our case, 2 days).

**it_conv:** This input takes the value of number of iterations after which to check if convergence criterion is satisfied. This takes a certain computational overhead, since general statistics have to be computed.

**it_gif:** This input takes the value of number of iterations after which 2D images of the geometry crossection (.gif) showing the density and current fluid configuration are to be saved.

**rho_vtk:** This input asks if 3D geometry files (.vtk) for both fluids 1 and 2 are to be saved: True, or only for for Fluid 1: False.

**it_vtk:** This input takes the  number of iterations after which the 3D geometry files (.vtk) showing the density and current fluid configuration are to be saved. If this number is greater than 100000 no vtk files will be output (which keeps the output lighter and the run faster).

**print_geom:** This input asks if a 3D geometry file (.vtk) is to be saved at the beginning of the simulation.

**print_stl:** This input asks if a 3D geometry file (.stl) is to be saved at the beginning of the simulation. This file may be used for 3D printing the geometry or for viewing.

*Note:* The .vtk files may be viewed in [PARAVIEW](https://www.paraview.org/).


# Chinese explanation of inputs for 2 phase flow XML file (input_sc.xml)

**两相流模型输入参数解释（input_sc.xml）

<img src="https://github.com/je-santos/MultiphasePorousMediaPalabos/blob/master/illustrations/LBM%20geometry%203D.png" align="middle" width="600" height="400" alt="geometry inputs explanation">


图1-两相流模型几何模型示意图
**安装xml文件里头的参数顺序进行解释**

**load_savedstated:** 该选项用于确定是否进行断点模拟（True）或者重新开始新模拟（False）

**geometry:** 几何选项下有多个输入参数

**file_geom:** 该选项要求输入几何文件名字，文件名以.dat格式命名，可以在前处理过程阶段创建 (https://github.com/je-santos/MultiphasePorousMediaPalabos/tree/master/pre-processing) 几何文件需放在该文件夹[2-phase simulation code](https://github.com/je-santos/MultiphasePorousMediaPalabos/tree/master/src/2-phase_LBM) 或者也可以放置于其他位置，但需要声明文件指向路径

**size:** 几何模型X，Y，Z方向上体素大小（图1）

*注意点:*

1) PALABOS 模拟X方向上的流动，所以需要注意模拟原始模型文件X方向和Z方向是都设置正确
2) 如果前处理过程增加了网格量和单相流体占据网格，原始X方向尺寸会变大

**init:** 初始化选项下有多个输入参数

**fluid1:** 该选项要求输入侵入相流体的初始三维空间位置，图1所示， x1 -x2, y1 - y2,  z1 - z2 （黄色标记）.

**fluid2:** 该选项要求输入原始赋存相，也就说被侵入相流体的初始三维空间位置，图1所示，  x1 - x2, y1 -y2, z1 - z2 （蓝色标记）.

**fluids:** 流体性质选项下有多个输入参数

**Gc:** 粒间作用力，控制流体与流体之间界面张力，推荐0.9取值，需满足 Gc > 1/(rho_f1+rho_f2).

**omega_f1 and omega_f2:** 用于计算流体1和流体2的运动粘度  v=( 1 / omega_fi - 0.5 ) / c^2.

**force_f1 and force_f2** 设置驱动力选项，如果该项取值不为0，驱动力会设置到每一相流体，例如重力，但模型边界需要设置为周期性边界

**Wetting forces:**
(G_ads_f1_s1, G_ads_f1_s2, G_ads_f1_s3, G_ads_f1_s4): 这些参数用于声明流体与壁面之间的相互作用力，可设置4种不同润湿性条件，体素值 1, 3, 5, 6 分别设置为 G_ads_f1_s1, G_ads_f1_s2, G_ads_f1_s3, G_ads_f1_s4, 体素值2 被设置为内部骨架, 4 被设置为中性润湿，0被设置为流体. 润湿角采用下式进行计算: cos(theta) = 4*G_ads_f1_si/( Gc*( rho_f1-rho_d ) )

<img src="https://latex.codecogs.com/svg.latex?\Large&space;cos(\theta)=\frac{4G_{ads_{f1,si}}-b\pm\sqrt{b^2-4ac}}{G_c(rho_{f1}-rho_d)}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

**rho_f1:** fluid 1的初始密度.

**rho_f2:** fluid 2 的初始密度.

*Note:* 推荐fluid1和fluid2的稳定值2，高密度比会变得数值不稳定

**pressure_bc:** 该参数用来确认是否采用压力梯度边界，True为非稳态流动， False 为稳态流动.

**rho_f1_i:** 入口压力边界fluid 1的固定密度.

**rho_f2_i:** 模拟开始时出口压力边界fluid 2的密度.

**rho_f2_f:** 模拟结束后出口fluid 2的密度，入口和出口压力边界决定毛管力。

**rho_d:** 两相流体在各相中的密度，推荐取值 0.06.

**drho_f2:** 低压力边界fluid2密度每一步（毛管力变化）降低的密度值，推荐取值0.01 到 0.1，减小值越小计算时间越长。

Shan-Chen模型计算压力:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;P(x)=\frac{\rho_1(x)+\rho_2(x)}{3}+G_c\frac{\rho_1(x)\rho_2(x)}{3}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

式中，第二项较小，可忽略不计

采用下式计算毛管力:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;P_c=P_{nw}-P_{w}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

可重新整理为:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;P_c=\frac{2}{3}-\frac{(2-\Delta\rho)}{3}=\frac{\Delta\rho}{3}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

式中第一项表示入口压力，第二项表示出口压力.

通过与圆形孔隙拉普拉斯方程对比, 可得到:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;\Delta\rho=\frac{6\sigma}{r}cos(\theta)" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

采用气泡算例，界面张力给定为0.15， G_ads_f1 = -0.4 (156.4度), 得到:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;\Delta\rho\approx\frac{0.825}{r}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />



**output:** 输出模块下有多个输入参数.

**out_folder:** 设置输出文件夹选项

**save_sim:** 选择是否将每一步毛管力变化下模型输出结果进行存储，可用于断点模拟

**convergence:** 收敛性判别标准，取值需参考示例

**it_max:** 达不到收敛性条件下，每个毛管力下可进行的最长迭代计算次数

**it_conv:** 每隔声明的迭代步判断是否收敛

**it_gif:** 声明迭代步下保存gif格式的二维两相密度分布

**rho_vtk:**是否保存vtk格式的fluid1和fluid2文件，true fluid1和fluid2都保存，false只保存fluid1

**it_vtk:** vtk格式文件输出的迭代步间隔，100000 以上不输出vtk格式文件

**print_geom:** 模拟开始情况下，是否输出vtk格式几何模型

**print_stl:** 模拟开始情况下，是否输出stl格式几何模型.

*Note:* The .vtk 文件可在 [PARAVIEW](https://www.paraview.org/)可视化
