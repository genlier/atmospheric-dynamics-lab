# Atmospheric Dynamics Lab

用傅里叶伪谱方法探索涡旋演化、地转适应与浮力驱动流动。包含 MATLAB / NumPy 数值模式，以及可直接在浏览器打开的交互展示。

**[打开交互实验室](https://genlier.github.io/atmospheric-dynamics-lab/)**

## 内容

| 算例 | 观察内容 |
| --- | --- |
| 衰减涡格 | 涡格黏性衰减；总动能与精确解对比 |
| 涡旋相互作用 | 同号涡旋互旋、合并与轴对称化 |
| 剪切不稳定 | 涡度带弯曲、卷起与涡核形成 |
| 涡旋漂移 | 行星涡度梯度、漂移与 Rossby 波 |
| 二维湍流 | 涡旋并合与动能谱演化 |
| 非对称双涡旋 | 幅值和半径差异对演化的影响 |
| 浅水地转适应 | 比较高度扰动、涡度扰动的调整 |
| 两层准地转 | 垂直切变与层间位涡耦合 |
| 热泡上升 | 均匀温度背景中的 Boussinesq 浮力过程 |

## 快速开始

### 网页

无需构建、登录或外部服务。直接打开 `docs/index.html`，或运行：

```bash
python -m http.server 8000 --directory docs
```

浏览 `http://localhost:8000`。支持算例切换、保存时次播放、物理量选择、色标切换、格点取值及动能/能谱查看。时间滑块浏览离散的数值快照，网页不重新积分方程。二维湍流的显示网格由 256×256 抽样至 64×64，诊断谱仍来自原网格；两层准地转算例仅提供末时刻。

### MATLAB

```matlab
run('models/matlab/bve.m')
```

在文件开头修改 `p.model_case`（1–5）选择测试。`models/matlab/shallow_water.m` 提供浅水模式，通过 `p.exp_case` 选择初始高度或涡度扰动。

### Python

```bash
python -m pip install -r requirements.txt
python models/python/bve.py --case 1 --output outputs/case1.npz
python models/python/experiments.py
python models/python/supplementary.py
```

`bve.py` 保存末时刻数值数组；衰减涡格同时保存数值与解析动能时序。`experiments.py` 包含参数探索和浅水算例；`supplementary.py` 包含两层准地转及热泡。辅助探索脚本的部分配置与 MATLAB 基准不同，查看各脚本开头的参数后运行，勿将不同配置的结果混为一次计算。

## 方法与参数

周期区域内采用傅里叶展开，空间微分转为谱系数与波数因子的乘积。正压模式以涡度为预报量，通过泊松反演求流函数；非线性乘积在物理空间计算，施加 2/3 截断退混淆。使用 RK4 推进平流，并以指数因子处理黏性扩散。

| MATLAB 测试 | 网格 | Re | beta | dt | T |
| --- | --- | --- | --- | --- | --- |
| 衰减涡格 | 64×64 | 1 | 0 | 0.000375 | 6 |
| 双涡旋 | 64×64 | 50 | 0 | 0.0000833 | 6 |
| 剪切层 | 64×64 | 100 | 0 | 0.0000833 | 6 |
| 漂移 | 64×64 | 50 | 1 | 0.0000833 | 6 |
| 湍流 | 256×256 | 40000 | 0 | 0.05 | 300 |

快照来自已保存的数值计算，未从论文或讲义图片提取。网页里的辅助模式与仓库脚本采用相同控制关系和参数。热泡采用双周期边界和均匀温度背景；它与上下平板加热的 Rayleigh–Bénard 算例具有不同边界和背景条件。

## 验证

```bash
python -m unittest discover -s tests
python tools/audit_publication.py
```

解析涡格测试检查动能衰减与精确解的一致性。静态检查检查个人路径、邮箱、常见密钥模式及不应发布的文件类型。

## 目录

```text
models/matlab/   MATLAB 求解器
models/python/   NumPy 求解器与辅助算例
docs/           GitHub Pages 页面和精简快照
tests/          解析解验证
tools/          发布内容检查
```

发布目录仅包含精选代码、说明和数值快照。个人报告、姓名学号、讲义文件、配置凭据、下载资料、原始大数据及历史缓存均不在仓库中。GitHub 提交使用通用项目署名。网页没有访问统计或第三方资源请求。
