## Explanation of inputs for 2 phase flow XML file (input_sc.xml)
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
