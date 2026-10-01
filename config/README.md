# 公共配置

| 名称 | 入口 | 行为 |
| --- | --- | --- |
| `reconstruct` | `reconstruct` | 显示重建过程和对焦结果 |
| `reconstruct_fast` | `reconstruct_fast` | 无窗口，默认只记录最后一次迭代或提前收敛时的结果 |
| `color` | `reconstruct_color` | 读取通道 MAT 清单并执行 RGB 重建与融合 |

`src/io/defaultReconstructionOptions.m` 定义公共重建默认值，灰度、彩色与仿真核心共用。

```matlab
cfg = loadProjectConfig('reconstruct');
cfg.options.focus.prior = 1.68e-3;
cfg.options.focus.halfRange = 0.6e-3;
cfg.options.focus.step = 0.01e-3;
cfg.options.output.rootFolder = fullfile(pwd, 'outputs', 'reconstruction', 'sample01');
reconstruct(cfg);
```

无第二参数时，先加载公共 `*_example.m`，再加载可选的 `local/config/*_local.m`。本地脚本只需覆盖字段，例如：

```matlab
% local/config/reconstruct_local.m
config.inputFile = fullfile(projectRoot, 'data', 'sample01.mat');
config.options.focus.halfRange = 0.4e-3;
```

`projectRoot` 由加载器提供。脚本不要修改它，也不要执行 `clear` 或 `cd`。推荐使用绝对路径或 `fullfile(projectRoot, ...)`。

使用 `loadProjectConfig(name, false)` 仅加载公共默认值，适合测试和可复现示例。入口接收完整配置，推荐加载后再修改，而非手工省略字段。

## 聚焦控制

- `focus.method`：`'fft'`（默认，原高通振幅 FFT 指标）或 `'adfrft'`（分数域强度指标）。旧配置未指定时使用 `'fft'`。
- `focus.frftOrder`：ADFrFT 阶数，默认 `0.8`，范围 `0 < p <= 1`。
- `focus.cropSize`：仅对 ADFrFT 生效，默认中央 `256×256` 像素；实际区域受图像尺寸限制，`0` 表示全图。先传播完整复场，再截取评价区域。
- `focus.roiMode`：`'center'` 使用上述中央区域；`'full'` 使用全图；`'manual'` 弹出预览供手动框选样品，忽略 `cropSize`。仅对 ADFrFT 生效。
- `focus.roiBounds`：手动模式留空会弹窗；也可填写已知的 `[xmin ymin xmax ymax]` 整数像素坐标，直接复用，不弹窗。
- `focus.prior/halfRange/step`：单位均为米；搜索范围为中心 ± 半宽。请按实际几何选择范围。
- `focus.halfRange = 0`：固定距离，仅有一个采样点，不代表进行了焦距搜索。
- `showFigures`：控制重建与聚焦窗口。
- `output.saveFocusPlot`：控制 `autofocus.png` 导出，默认开启；聚焦数据 MAT 始终保存。
- 聚焦在记录迭代时执行，`recordEvery` 越小，搜索次数和输出量越多。

例如在 `local/config/reconstruct_local.m` 中添加以下字段，可让无参数的 `reconstruct` 使用分数域指标：

```matlab
config.options.focus.method = 'adfrft';
config.options.focus.frftOrder = 0.8;
config.options.focus.cropSize = 256;
```

快速和彩色入口使用同样字段。两种指标的数值尺度不同，不能直接比较分数大小来判断哪种方法更好。默认原指标评价全图，ADFrFT 默认评价中央区域；比较结果时需同时考虑评价区域的差异。

## 彩色控制

本机直接修改根目录 `color_config.m`，包含输入、重建、聚焦、约束、输出和融合设置；它不读取 `reconstruction_config.m`。公开包中只需修改 `config/color_example.m`，再运行 `reconstruct_color` 或 `reconstruct_color_fast`。彩色入口在没有根目录配置时仅加载公共示例，不自动叠加旧 `local/config/color_local.m`。

批处理可以显式传入完整配置，不受本机彩色配置影响：

```matlab
cfg = loadProjectConfig('color', false);
cfg.inputFile = 'D:\my_data\color_input.mat';
cfg.options.iterations = 100;
cfg.options.adaptiveConstraint.enabled = true;
folder = reconstruct_color_fast(cfg);
```

`inputFile` 指向通道清单；`channelNames` 标识清单内通道，默认 `{'B','G','R'}`；`rgbOrder = [3,2,1]` 指定合成 RGB 时的索引。波长、像素尺寸、距离和几何信息来自各通道 MAT。

`reuseExistingReconstructions` 控制缓存复用；公共默认也读取 `COLOR_FORCE_RECONSTRUCT=1` 环境变量。修改输入 MAT 或 APRW 参数后会重新重建。原始图像变化后应由输入适配器重新生成 MAT。
