# 大气动力学实验室

**[英文主页](https://genlier.github.io/atmospheric-dynamics-lab/)** · **[中文网页](https://genlier.github.io/atmospheric-dynamics-lab/?lang=zh)** · [English README](README.md)

用傅里叶谱方法探索涡旋演化、地转适应与浮力驱动流动。英文为网页默认语言，可随时切换中文。

## 自由实验

支持正、负涡度画笔、橡皮擦、数学表达式和 JSON 网格导入。可设置雷诺数、beta、积分时长、最大时间步长和保存帧数，随后在浏览器里计算并播放动画。

周期区域为 $[-\pi,\pi)\times[-\pi,\pi)$，提供 32×32 和 64×64 网格，保存 61、121、601 或 1201 个计算帧。播放支持暂停、时间轴拖动、速度调整、色标切换和格点取值，并显示动能与涡度拟能。初值和计算结果可导出为 JSON，数据保留在当前浏览器。

自由实验计算正压涡度方程。九个参考算例还包括浅水地转适应、两层准地转和热泡上升，这些补充模式提供已有计算结果。

每个参考动画均从初值完整积分到终止时刻，保存 601 张实际计算 PNG 帧，再按 30 fps 顺序编码。二维湍流动画使用原来的 256×256 网格。自由实验同样直接播放计算保存帧，帧数设置决定记录密度。用于初值编辑的小型预览保留 64×64 网格。

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



## 复现录制动画

安装 requirements-animation.txt，并将 MATLAB、FFmpeg 加入命令行路径后运行：

~~~bash
matlab -batch "addpath('tools'); run_dense_bve"
python tools/run_dense_supplements.py
python tools/render_dense_animations.py case1 case2 case3 case4 case5 asymmetric shallow-water qg thermal
python tools/update_reference_data.py
python tools/verify_animations.py
~~~

记录脚本调用现有算法并保存密集时次，bve.m 原文件保持不变。PNG 序列保存在 outputs/dense/frames，视频及时间元数据保存在 docs/animations。网格、时间步长与终止时刻沿用原模型设置。
