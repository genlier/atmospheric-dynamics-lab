# 大气动力学实验室

**[英文主页](https://genlier.github.io/atmospheric-dynamics-lab/)** · **[中文网页](https://genlier.github.io/atmospheric-dynamics-lab/?lang=zh)** · [English README](README.md)

用傅里叶谱方法探索涡旋演化、地转适应与浮力驱动流动。英文为网页默认语言，可随时切换中文。

## 自由实验

支持正、负涡度画笔、橡皮擦、数学表达式和 JSON 网格导入。可设置雷诺数、beta、积分时长、最大时间步长和保存帧数，随后在浏览器里计算并播放动画。

周期区域为 $[-\pi,\pi)\times[-\pi,\pi)$，提供 32×32 和 64×64 网格，保存 61、121 或 241 个计算帧。播放支持暂停、时间轴拖动、速度调整、色标切换和格点取值，并显示动能与涡度拟能。初值和计算结果可导出为 JSON，数据保留在当前浏览器。

自由实验计算正压涡度方程。九个参考算例还包括浅水地转适应、两层准地转和热泡上升，这些补充模式提供已有计算结果。

参考算例的平滑播放对已有快照进行插值。自由实验由数值积分生成较密集的计算帧，再对相邻帧插值播放。五个正压算例可通过“将此初值载入自由实验”重新计算。自由实验采用页面显示的网格和参数；二维湍流参考结果由 256×256 抽样至 64×64 显示。

## 运行

在仓库目录运行：

~~~bash
python -m http.server 8000 --directory docs
~~~

浏览 http://localhost:8000。全部脚本随仓库提供，计算通过 HTTP 或 HTTPS 网页启动。

MATLAB 正压模式通过文件开头的 p.model_case 选择 1–5 号测试：

~~~matlab
run('models/matlab/bve.m')
~~~

浅水模式位于 models/matlab/shallow_water.m，通过 p.exp_case 选择初始高度或涡度扰动。

Python：

~~~bash
python -m pip install -r requirements.txt
python models/python/bve.py --case 1 --output outputs/case1.npz
python models/python/experiments.py
python models/python/supplementary.py
~~~

前者保存数值数组及衰减涡格的数值、解析动能序列；后两者涵盖非对称双涡旋、浅水、两层准地转与热泡，参数见源码。

## 数学模型

$$
\frac{\partial\zeta}{\partial t}+J(\psi,\zeta)+\beta\frac{\partial\psi}{\partial x}
=\frac{1}{\mathrm{Re}}\nabla^2\zeta,\qquad
\nabla^2\psi=\zeta.
$$

速度为 $u=-\partial_y\psi$、$v=\partial_x\psi$。傅里叶微分用波数因子表示空间导数，泊松反演采用零均值涡度，并令流函数零模态为零。非线性乘积采用 2/3 截断，平流项用 RK4 推进，黏性扩散用指数因子处理。浏览器时间步长由设置上限、速度 CFL 估计和 beta 波频率共同约束。

初值 JSON 包含 n 和 values。数值按 x 优先、随后 y 递增排列，坐标从 −π 开始。结果文件还包含初值、参数、扣除的均值、完成状态，以及各计算帧的时刻、涡度、动能和涡度拟能。

## 验证

~~~bash
node tests/test_browser_solver.js
python -m unittest discover -s tests
python tools/audit_publication.py
~~~

检查涵盖 FFT 反演、解析衰减、beta 波相位、非线性动能与涡度拟能守恒，以及表达式解析。网页交互验证涵盖中英文切换、绘图、计算、导入导出、取消计算和手机布局。

仓库只收录精选源码及数值快照。个人报告、讲义文件、凭据、私人路径和原始大数据均未纳入，提交使用账号的 GitHub 隐私邮箱。

