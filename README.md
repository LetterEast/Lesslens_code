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

## 输出

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
├─ coverage_confidence.png
├─ tv_risk.png
├─ tv_lambda.png
└─ meta.mat
```

相位热力图使用中心视场相位的 1%–99% 分位范围抑制边缘离群值，再使用 `hot` 色表增强物体相位特征。实际显示范围保存在 `reconstruction.mat` 的 `phaseDisplayLimits` 和 `originalPhaseDisplayLimits` 中。重建默认启用样品面自适应复数 TV，最大重叠区域受到保护，低覆盖边缘约束更强。

`results/reconstruction.mat` 同时保存 APRW 得到的参考面复场 `field`、TV 后样品面结果 `object` 和严格原始视场结果 `objectOriginalFOV`。因此可以在不重新运行 APRW 的情况下重新自动聚焦或改变反向传播距离。

需要 MATLAB Image Processing Toolbox；有可用 GPU 时会自动加速，否则使用 CPU。
