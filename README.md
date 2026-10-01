# 无透镜成像：标定与重建

在 MATLAB 中打开项目根目录，先运行 `setup_project`。

| 用途 | 过程可视化版 | 快速版 |
|---|---|---|
| 仿真：球面波、标定与重建比较 | `demo_sim` | `demo_sim_fast` |
| 实验：读取数据并重建 | `demo_exp` | `demo_exp_fast` |
| 彩色：三通道重建与融合 | `reconstruct_color` | `reconstruct_color_fast` |

每对入口使用相同的数值设置。快速版关闭过程窗口，仍保存结果和对比图；仿真两版使用相同求解器，并非两种算法。

## 仿真：只修改主入口

打开 **`main/demo_sim.m`**，顶部依次设置光学参数、采集位置、标定参数、样品和迭代次数。无需另建配置文件；`demo_sim_fast` 直接复用这些参数。

```matlab
setup_project
[result, folder] = demo_sim;  % 显示过程和最终对比
% [result, folder] = demo_sim_fast;  % 同样计算，关闭窗口
```

### 选择素材

修改主入口的 `options.sample.type`：

| 值 | 内容 | 用途 |
|---|---|---|
| `'mixed'` | 吸收条纹、吸收斑点和光滑相位 | 默认幅相混合样品 |
| `'resolution'` | 不同宽度的横纵条纹 | 纯吸收细节对比；非标准 USAF 标定靶 |
| `'cells'` | 三个细胞状相位轮廓 | 纯相位样品 |
| `'widefield'` | 分布在较大范围的吸收点与相位斑点 | 检查新增视场的恢复质量 |
| `'image'` | 用户提供的图片 | 自定义幅度和相位 |

内置样品由代码生成，不依赖外部素材。自定义图片模式设置：

```matlab
options.sample.type = 'image';
options.sample.amplitudeFile = 'D:\my_data\amplitude.png';
options.sample.phaseFile = 'D:\my_data\phase.png';
```

幅度图的灰度值映射到幅度 `[0,1]`，白色透明、黑色不透光；相位图映射到 `[0, phaseScale]` rad。图片按比例缩放到物面，图片范围之外的复场为零；不做自动对比度拉伸。幅度路径为空表示图片范围内幅度为 1，相位路径为空表示相位为 0。两张图片应事先配准、具有相同宽高比。可使用自己的 BMP、PNG、TIFF 或 JPEG，例如本地已有的 CCTV 仿真素材。

`options.sample.sizePixels=512` 独立设置图片最长边；`options.imageSize` 只设置相机窗口边长，两者不再绑定。旧的 `sizeFraction` 参数已移除。对 256×256 的 cameraman 图，这是 2 倍上采样；插值放大用于构造大样品，不会凭空增加原始纹理信息。

例如保持 512 像素样品，把窗口改成 300×300：

```matlab
options.sample.sizePixels = 512;
options.imageSize = 300;
options.padding = 320; % 计算网格边长为 300+2*320=940
```

改为 128×128 或 200×200 窗口时，只需修改 `imageSize`。样品超过计算网格时，报错会列出实际尺寸；不要为了改变相机窗口而放大样品。

相机窗口超过重采样后的图片尺寸时会自动缩小。当前相机为方形，非正方形图片按短边限制；幅度图决定图片范围，仅提供相位图时以相位图为准。实际值和请求值分别保存为 `result.options.imageSize`、`requestedImageSize`。

倾斜照明导致某个相机像素的几何采样位置移出图片时，该像素在测量输出和采集示意中置零，`result.cameraMasks` 同时将其标为无效。噪声不会填入这些越界位置。重建不对无效位置施加零强度约束，避免把缺失数据误当成真实测量；所有模型共用同一采样掩膜。简化求解器每轮将图片支撑之外的物面场固定为零；实验核心的差别见下文。

这是显式的有限图片采样规则：真实衍射可能越过几何边缘，不能用补零推断那里的物理光强必然为零。图像模式现在直接传播有限支撑内的复场，标定点阵保留独立的透明背景模型。评价只计算图片范围内的覆盖区域，`fov_summary.csv` 同时给出几何面积和样品范围内的有效面积，避免把零填充算作新增信息。

`options.sample.offsetPixels=[x,y]` 可独立移动图片在物面上的位置（整数像素），不改变光源位置。

### 大样品与有限相机：移动后看到其它区域

仿真先在较大物面网格上生成完整样品，再传播到各个探测器高度，最后截取固定大小的相机窗口。每张测量图仍是 `imageSize × imageSize`，不能把完整物面图直接当相机图。

主入口中可以独立调整样品和相机。运行中央 128×128 窗口示例：

```matlab
setup_project
[result, folder] = demo_sim(struct('imageSize',128));
```

查看 `04_camera_acquisition.png`：上排为完整样品与逐帧物面覆盖红框；中排是窗口对应的清晰幅度真值，只作示意；下排才是传播后相机实际记录的强度。静态图选取 6 个代表位置，`camera_acquisition.gif` 展示全部 13 帧。`options.fov.saveAcquisitionGif=false` 可关闭动画保存。

128 窗口示例：样品到相机距离从 0.5 mm 增至 12.5 mm，光源到样品 30 mm，光源偏移 `[300,-180]` 像素。样品中心对准首帧物面覆盖中心，后续帧逐渐看到右上侧的其它部分。传播网格为 768×768，未覆盖的样品仍参与前向衍射；每帧仅用中央 128×128 强度作重建约束。独立点阵标定使用完整 192×192 图，避免点阵被小窗口截断过多。

这里相机保持横向位置、沿轴向改变距离。离轴球面照明使物面覆盖同时缩放和偏移。红框是几何范围，真实强度由波传播计算，包含边缘衍射；当前相机使用网格点强度采样，尚未模拟像素面积积分、满阱饱和或读出噪声。

### 改大倾斜角为什么会失败

当前角谱实现直接在间距为 `pixelSize` 的网格上采样复场。其局部载频必须满足每轴 `abs(sin(theta))/wavelength < 1/(2*pixelSize)`，这只是该数值实现的必要条件，并不是同轴强度成像系统的普适倾角上限；样品频谱也会占用带宽。

例如 `pixelSize=3 um`、`wavelength=514 nm`、`sourceToSample=15 mm`、`M=1000` 时，横向光源距离为 3 mm，中心倾角约 11.31 度，复场载频约 1.14 周期/像素，超过 0.5 的采样界限。直接离散化会混叠。要可靠模拟这种倾角，需要去载频/移位频谱传播，或更细的计算网格并单独建模相机采样；仅增加补边不能解决。

`checkSimulationSampling` 现在会检查复场载频、标定中心点是否移出相机、物面覆盖是否超出计算网格。仿真还利用已知真值检查标定是否严重失败，避免错误几何继续生成貌似合理的视场图。参数不满足时会报出原因，不会自动替换用户设置。采样背景可参考 [off-axis diffraction 的频谱采样研究](https://www.eee.hku.hk/optima/pub/journal/2307_OPT.pdf)。

### 客观评价与视场扩展

主对比图标出联合视场的幅度 PSNR、SSIM；`03_quality_metrics.png` 展示各区域评价。`quality_by_region.csv` 保存各组重建的完整指标：

- 幅度 RMSE、PSNR（固定峰值 1）、SSIM、复场 NRMSE。
- 圆周相位 RMSE、相位 PSNR（固定范围 `2*pi`）、相位 phasor SSIM。
- 评价像素数、有效 SSIM 窗口数，以及最终重建场的测量幅度残差。

评价分别覆盖 `first_fov`（首帧）、`union_fov`（多帧并集）、`added_fov`（并集减首帧）、`overlap_fov`（全部帧共同覆盖）。各组使用同一真值几何掩膜，避免不同裁剪范围影响公平性。指标在原始数值上计算，不做显示归一化或幅度增益拟合。

相位先按整个联合视场消除一次全局相位偏置，各区域共用该偏置；`rawFields` 保留未对齐的重建。真值幅度不超过 `phaseAmplitudeThreshold` 的位置不评价相位。相位 SSIM 比较正弦/余弦表示，避免 `-pi/pi` 跳变，并非直接对相位灰度图计算的 SSIM。幅度与相位 SSIM 都只使用完全落在评价区域内的 11×11 窗口；无新增区域或区域太窄时相应指标为 `NaN`。

SSIM 使用 MATLAB 的 [ssim](https://www.mathworks.com/help/images/ref/ssim.html)，固定 `DynamicRange=1`、`Radius=1.5`。完美重建的 PSNR 为 `Inf`。

视场计算直接复用实验的 `sampleGeometrySupport`。原始仿真图没有做平移配准，因此配准剪切设为零。光源投影位于传感器内部时，各高度的物面视场可能互相包含，扩展面积为零；增加 `padding` 只扩大计算网格，不增加测量信息。

运行一个确实产生几何扩展的例子（不修改主入口当前设置）：

```matlab
[result, folder] = demo_sim(struct('iterations',60)); % 或 demo_sim_fast(...)
```

`02_expanded_fov.png` 对比同一多帧重建的首帧裁剪与联合视场，显示覆盖帧数和新增区域误差；这不是单帧求解与多帧求解的对比。灰色区域没有几何覆盖，不计入有效重建。`fov_summary.csv` 同时记录真值、标定估计的首帧与联合面积。`options.fov.displayUnion` 决定主图显示联合视场还是首帧。

扩展功能复用实验的几何覆盖与显示规则。前三组使用简化幅度投影，第四组直接调用实验核心 `reconstructMultiPlane`。几何覆盖增加不保证边缘恢复质量与中心相同，应结合新增区域指标判断。

### 仿真具体比较什么

```text
demo_sim 顶部参数
  ├─ 独立点阵 + 离轴球面波 → 衍射图 → calibratePointSource → 估计 M/N/Z
  └─ 被测样品 + 同一真实球面波 → 多帧强度数据
                                    ├─ 已知真值球面波重建（参考）
                                    ├─ 标定后球面波重建（估计 M/N/Z）
                                    ├─ 平面波 + 平移重建（对照）
                                    └─ 实验核心重建（同一标定 M/N/Z，所选 TV / 自适应约束）
```

四组共用测量数据和迭代次数，前三组共用简化求解器。`06_core_vs_calibrated.png` 并排显示真值、标定球面波简化重建、实验核心重建；`amplitude_phase_comparison.png` 显示全部四组。

`sourceToSample` 是光源到样品距离；标定 `Z` 是光源到首个探测器平面的距离。`sampleDistances` 和 `calibrationDistances` 是绝对传播距离，两组首个距离须相同。

前三组用简化幅度投影隔离传播模型和标定的影响，不包含 TV 或自适应约束。第四组保留实验核心的反馈迭代、覆盖权重融合及所选约束，TV 按核心原有流程在输出阶段应用。所有组固定使用已知样品距离，本次比较不混入自动聚焦误差。相位恢复仍可能有误差，应结合相位指标判断。

在 `main/demo_sim.m` 顶部设置第四组，无需修改实验配置：

```matlab
options.core.enabled = true;
options.core.tvEnabled = true;
options.core.adaptiveConstraintEnabled = true;
options.core.adaptiveConstraintStrength = 0.05;
options.core.adaptiveConstraintEdgeWidth = 32;
```

仿真帧直接嵌入计算网格，不重复做实验平移配准。输入适配保留标定 M/N/Z，将绝对采集距离转换为逐次位移；输出除以估计球面波的径向幅度，得到与其它组一致的透射率单位。核心使用仿真照明单位设置白场幅度，不拟合样品亮度。完整设置与迭代残差保存在 `experimental_core/` 下。

图像样品的支撑处理有一项差别：简化算法每轮将图片范围外的物场置零；实验核心保留原有迭代约束，只在最终展示时将该范围外置零。因而这是两套完整流程的比较，不能把结果差异全部归因于某个单独约束。四组指标使用相同真值区域，最终测量残差使用各自未经展示裁剪的重建场计算。

`convergence.png` 保留各求解器的内部残差：核心曲线位于参考探测器面、在输出 TV 之前，简化算法曲线由物面场正向计算。比较最终数据拟合时请使用 `metrics.csv` 中按同一评分流程计算的 `amplitude_residual`，不要直接比较两种内部曲线的高低。

### 本地调试与批量实验

无需复制主入口，只传入要修改的字段；未指定的参数沿用主入口。拼错参数名会报错。

```matlab
o = struct('iterations', 120);
o.sample = struct('type', 'cells', 'phaseScale', 0.8);
[result, folder] = demo_sim_fast(o);
```

私人批量脚本放在 Git 忽略的 `local/experiments/`；使用 `setup_project('local')` 加入路径。主入口适合保留可公开复现的参数。旧的仿真配置已移到 `local/archive/simulation_settings/`，当前流程不读取它们。

## 实验：输入适配与核心计算分开

```text
个人原始文件 → 个人读取/预处理 → 标准 inputData → 核心重建 → 结果
```

本地采集使用根目录 **`reconstruction_config.m`** 设置路径、波长、像素尺寸、位移、聚焦方法和重建选项，然后运行 `demo_exp` 或 `demo_exp_fast`。`prepareInput=true` 重新准备输入，`false` 复用已有 MAT。`distanceSteps` 表示逐次移动量，代码累加成位置。

公开工程的使用者按 [输入格式](docs/INPUT_FORMAT.md) 自行读取数据并保存标准 MAT，再调用：

```matlab
cfg = loadProjectConfig('reconstruct', false);
cfg.inputFile = 'my_prepared_input.mat';
cfg.options.focus.prior = 1.68e-3;
[fields, folder] = demo_exp_fast(cfg);
```

`view_reconstruction(folder)` 查看结果，`reconstruct_color` 处理彩色重建。原始文件名、设备目录和个人分析留在适配层，不进入核心算法。

### 彩色与灰度共用核心

本机彩色入口优先读取根目录 `color_config.m`：图片根目录、B/G/R 子目录、标定根目录、波长、像素尺寸和相邻位移都在这里设置。`prepareInput=true` 自动生成通道 MAT；只改算法参数时可设为 `false` 复用输入，然后运行 `reconstruct_color`。

`folder = reconstruct_color_fast;` 使用同一配置和算法，关闭过程窗口及聚焦曲线 PNG 导出，仍保存彩色图、各通道灰度图、复场和聚焦数值。迭代次数和记录间隔保持配置值；若选择手动聚焦且未填写 ROI，仍需手动选区。

彩色的迭代、聚焦、TV、自适应约束、灰度显示和融合参数全部在 `color_config.m` 独立设置，不读取灰度配置。彩色和灰度共用 `defaultReconstructionOptions` 与 `reconstructMultiPlane`，三个通道分别重建。测量数据和标定保持每通道独立；额外逐帧亮度校正默认关闭。没有本机配置时，公开入口使用 `config/color_example.m` 和标准通道 MAT。

结果位于 `outputs/color/`：`channel_B/`、`channel_G/`、`channel_R/` 保留完整重建过程与供融合使用的画布；`grayscale/B/`、`grayscale/G/`、`grayscale/R/` 集中保存各通道最后一次记录的灰度结果。每个灰度目录都有 `amplitude.png`、`phase_heatmap.png`、联合视场和首帧视场图片，以及数值 `reconstruction.mat`。这些图使用单独灰度重建的裁剪、掩膜、幅度归一化和相位色轴规则，在 RGB 配准和颜色调整之前导出，不重复传播或重建。MAT 中 `sourceResult` 指向原始完整结果。

相同输入和参数才能复用通道缓存，改变约束模式会重新计算。`color_config.m` 被 Git 忽略；旧的 `local/config/color_local.m`、`color_input_local.m` 不参与新的根目录配置流程。

彩色各通道每次记录迭代时，在其 `iteration_XXXX/results/` 中额外保存 `amplitude_cropped.png` 与 `phase_heatmap_cropped.png`。使用与单独灰度重建一致的有效视场裁剪和显示规则，去掉计算画布的外围大黑边；原有完整画布图片和复场仍用于颜色融合。

**给别人使用：**公开包不包含你的本地配置和输入适配器。使用者只需在 `config/color_example.m` 设置标准通道清单路径、重建参数和融合参数，然后运行 `reconstruct_color` 或 `reconstruct_color_fast`。清单和通道 MAT 的格式见 [输入格式](docs/INPUT_FORMAT.md)；各通道的波长、像素尺寸和标定来自输入 MAT。用户可按自己的相机格式编写输入适配器，核心无需修改。

配置结构表示一次彩色重建任务：`inputFile` 描述输入，`options` 描述重建与灰度导出，`colorBalance` 和 `physicalIntensity` 描述融合。入口接收完整配置，显式传参优先于本机默认文件，适合多个样品的批处理。各任务的参数独立，避免彩色设置随灰度配置变化。

## 标定

本地修改 `calibration_config.m`，运行 `demo_calibration`。公开示例见 `config/calibration_example.m`。

标定基于原 `circle_search_2.0/mainV5` 的背景归一化、投影寻峰、网格拟合和几何回归。无效帧被跳过时保留真实采集位置。M/N 的像素原点采用 `floor(size/2)+1`，与原 mainV5 的 `size/2` 有一个像素的约定差异。输出 `calibration.mat`、兼容输入的 `MNZ_result.mat` 和拟合图。

## 目录与输出

| 目录 | 职责 |
|---|---|
| `main/` | 用户入口；仿真参数也在这里 |
| `src/reconstruction/`、`src/optics/` | 重建与传播 |
| `src/calibration/`、`src/focus/` | 标定与聚焦 |
| `src/geometry/`、`src/constraints/` | 有效视场与约束 |
| `src/simulation/` | 样品生成与模型比较 |
| `src/io/`、`src/utils/` | 标准输入、结果保存与通用辅助 |
| `src/compatibility/` | 旧接口兼容入口 |
| `config/`、`examples/` | 公开配置示例和标准输入示例 |
| `local/` | 私人输入适配、调试、分析和旧文件；Git 忽略 |
| `outputs/` | 运行结果；Git 忽略 |

每次运行创建独立目录：

```text
outputs/simulation/model_comparison_时间戳/
  01_calibrated_vs_plane.png       主要幅相对比
  02_expanded_fov.png              首帧、联合视场与新增区域误差
  03_quality_metrics.png           分区域客观评价
  04_camera_acquisition.png        完整样品、逐帧覆盖及有限相机强度
  camera_acquisition.gif           相机移动采集动画（可关闭）
  simulation.mat                  真值、重建、实际参数和测量数据
  06_core_vs_calibrated.png        实验核心与标定球面波的幅相对比
  experimental_core/              实验核心运行设置、最终复场与诊断
  metrics.csv                     各组误差
  quality_by_region.csv            分区域 PSNR、SSIM、相位误差等
  fov_summary.csv                  真值及标定几何的视场面积
  calibration_parameters.csv      M/N/Z 真值与估计
  amplitude_phase_comparison.png  含已知真值参考组的总览
  simulated_measurements.png      测量数据
  calibration_fit.png             标定拟合
  convergence.png                 收敛曲线
outputs/reconstruction/reconstruction_时间戳_framesN_iterN_recordN/
  run_config.mat / run_manifest.json
  iteration_0020/results/ / diagnostics/
outputs/calibration/calibration_时间戳/
```

MAT 文件保留数值结果；图片仅用于查看。默认对比图的幅度色轴为 `[0,1.2]`、相位为 `[-1,1]` rad，改变相位强度后应结合 MAT 数据查看。

## 依赖、检查与公开导出

在 MATLAB R2025b 验证；图像处理需要 Image Processing Toolbox，标定寻峰需要 Signal Processing Toolbox。实验核心 GPU 运行需要 Parallel Computing Toolbox；仿真前三组使用 CPU，第四组沿用实验核心的 GPU/CPU 自动选择。

```matlab
setup_project
addpath('tests')
run_core_tests
verify_workflow_profiles
```

`local/tools/export_core.ps1` 按白名单导出公开工程。个人配置、数据、分析、结果和归档不随核心导出。更多细节见 [重建说明](docs/RECONSTRUCTION.md) 和 [配置示例](config/README.md)。
