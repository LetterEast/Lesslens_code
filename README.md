# Lensless APRW reconstruction

项目按照 CCTV-phase-retrieval 的层级重新整理，数据创建和重建完全分离。

## 目录

```text
Lesslens code/
├─ main/
│  ├─ create_input_data.m   # 创建统一重建输入
│  ├─ reconstruct.m         # 标准重建
│  └─ reconstruct_fast.m    # 无窗口快速重建
├─ src/
│  ├─ APRW.m                # 重建、自动聚焦和结果保存
│  ├─ prepareMeasurements.m # 图像读取与配准
│  └─ propagate.m           # 同轴角谱传播
├─ data/
│  ├─ calibration/          # MNZ 标定数据
│  └─ reconstruction_input.mat
├─ tools/                   # 与主重建无关的实验脚本
└─ ResultFolder/            # 重建结果
```

## 使用方法

### 1. 创建输入数据

编辑 `main/create_input_data.m` 顶部的图像目录、标定文件、波长、像素尺寸和距离步长，然后运行：

```matlab
cd main
create_input_data
```

程序生成 `data/reconstruction_input.mat`。输入图像应当已经完成项目外的预处理。

### 2. 重建

带实时显示：

```matlab
reconstruct
```

无窗口快速版本：

```matlab
reconstruct_fast
```

两个重建入口只读取 `data/reconstruction_input.mat`，不会访问原始图像或重新配准。

融合权重按每个像素的覆盖数 C 分配：覆盖它的每个平面占 1/C，未覆盖的占 0。已取消四边余弦衰减和基于软权重阈值的测量支持截断；边缘及单平面覆盖区仍完整参与。权重从硬掩膜重新计算，已有输入 MAT 无需重新生成。局部覆盖自适应不改变现有振幅上限、TV 或输出裁剪流程。

## 输出

### 球面相位移除对照测试

在 `main` 中运行 `test_spherical_removal`，读取最近保存的参考面复场，不重复迭代。也可传入指定的 `Iter_XXXX` 路径。对照分支直接回传，测试分支先去除参考面球面相位再按同一聚焦距离回传；两路均不做 TV、振幅截断或支持掩膜清零。输出在该迭代的 `sphere_removal_test/时间戳/`，两路完整画布及相同范围裁剪图分别保存，不拼图，使用统一显示范围，并保存复场 MAT。移除球面相位后不再去除样品面球面相位。本测试改变了传播条件，不能仅凭清晰度变化判断放大率；未调整测试分支的等效传播距离或重新聚焦。

在 `main` 中运行 `test_geometry_crop` 可单独测试几何裁剪：保持球面波正常回传，按保存的几何参数计算每个平面的放大率 `(Z+zi)/(Z-focusDistance)`，将相机像素边界（含配准平移）映射到物面后取并集。配准坐标中光源横向位置使用保存的 `M/N`，避免重复加入已由配准处理的偏轴位移。输出在 `geometry_crop_test/时间戳/`，原裁剪和几何裁剪各自保存，均为 TV 前振幅且显示范围相同。并集内没有几何支持的角落直接零填充为黑色；另存支持掩膜、逐平面边界、覆盖图和放大率。测试与正式流程共用 `src/sampleGeometrySupport.m`，采用点光源、平行平面和现有配准坐标模型；几何支持并非严格的衍射信息边界。

结果保存在 `ResultFolder/APRW_*/Iter_XXXX/`，并分为：

```text
results/
├─ reconstruction.mat
├─ amplitude.png
├─ phase_heatmap.png
├─ originalFOV_amplitude.png
└─ originalFOV_phase_heatmap.png

diagnostics/
├─ autofocus.mat
├─ adaptive_tv.mat
├─ sample_geometry.mat
├─ coverage_confidence.png
├─ tv_risk.png
├─ tv_lambda.png
└─ meta.mat
```

相位热力图使用中心视场相位的 1%–99% 分位范围抑制边缘离群值，再使用 `hot` 色表增强物体相位特征。实际显示范围保存在 `reconstruction.mat` 的 `phaseDisplayLimits` 和 `originalPhaseDisplayLimits` 中。重建默认启用样品面自适应复数 TV，最大重叠区域受到保护，低覆盖边缘约束更强。

APRW 在测量区域内使用局部覆盖权重融合，在测量区域外保留传播预测场。输出采用相机边界映射到样品面后的几何联合视场，不按固定像素数内缩。振幅和相位 PNG 按同一有效掩膜将无效区域置零，不保存透明通道，显示分位统计排除无效区域。现有 TV 参数及权重算法保持不变。可通过 `options.output.trimPropagationBoundary` 额外内缩传播边界，宽度记录在 `propagationMargin` 中。`originalFOV_*` 已改为第一相机映射到物面的几何视场，不再保持相机原始像素宽高，且不做额外传播边界内缩。

`results/reconstruction.mat` 同时保存 APRW 得到的参考面复场 `field`、TV 后样品面结果 `object` 和第一相机物面视场结果 `objectOriginalFOV`。因此可以在不重新运行 APRW 的情况下重新自动聚焦或改变反向传播距离。

需要 MATLAB Image Processing Toolbox；有可用 GPU 时会自动加速，否则使用 CPU。

正式重建默认启用几何输出，无需新增开关或重新创建输入 MAT。`reconstruction.mat` 新增与 `object/validMask` 同尺寸的 `sampleCoverage` 以及 `originalValidMask`；完整几何信息位于 `diagnostics/sample_geometry.mat`。旧单侧修剪逻辑及其 `samplePlaneTrimPixels/objectPlaneOffsetPixels` 输出已移除，`bounds/originalBounds` 仍为完整计算画布中的像素坐标。

彩色入口 `tools/color.m` 为每个波长保留完整物面画布用于配准，读取对应几何 `validMask/sampleCoverage` 并同步平移，最后取三通道共同支持范围。RGB 与单通道 PNG 均将无效区域零填充为黑色，不保存透明通道，不拼图。缓存版本已更新，旧裁剪缓存不会被复用。曝光校正、颜色平衡和 TV 参数不受本次修改影响。

可运行 `tools/verify_geometry_output` 检查零距离几何、三种 APRW 输出模式、振幅及相位无效区域零填充，以及已有几何测试结果的逐像素复现。
