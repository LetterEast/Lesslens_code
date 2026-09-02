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

结果保存在 `ResultFolder/APRW_*/Iter_XXXX/`：

- `Object_amp_*`、`Object_phs_*`：扩大视场结果；
- `Object_originalFOV_*`：严格裁回第一张输入图尺寸的结果；
- `Object_trusted_*`：至少两个测量面覆盖的可靠区域；
- `autofocus.mat`：自动聚焦距离与评价曲线；
- `meta.mat`：重建参数。

需要 MATLAB Image Processing Toolbox；有可用 GPU 时会自动加速，否则使用 CPU。
