# 重建输入接口

核心接收**已配准到共同画布的强度测量**。目录扫描、设备标定文件解析、前景提取、图像配准等由输入适配器完成。

## 灰度输入

MAT 文件内保存名为 `inputData` 的标量结构体：

| 字段 | 格式与含义 |
| --- | --- |
| `images` | 至少两个二维 `single/double` 数组组成的 cell；尺寸相同，有限、非负的强度值，不是场振幅 |
| `distanceSteps` | 每张图像对应一个相邻轴向间隔，单位 m；平面绝对位置由 `cumsum` 得到 |
| `wavelength` | 波长，单位 m，正标量 |
| `pixelSize` | 像素尺寸，单位 m/pixel，正标量 |
| `geometry` | 下表所示几何结构体 |

例如间隔 `[0, 0.2, 0.2] * 1e-3` 对应平面位置 `0、0.2、0.4 mm`。通常第一个间隔为零。

### geometry

| 字段 | 含义 |
| --- | --- |
| `M / N` | 当前已配准坐标中点光源的横向位置，单位 pixel；现有剪切配准后通常均为 0 |
| `Z` | 点光源到参考探测面的距离，单位 m，必须为正 |
| `orig_M / orig_N` | 配准前点光源的横向位置，单位 pixel；用于传感器视场映射 |
| `orig_size` | 原始相机 `[height,width]` |
| `padSize` | 对称填充 `[padY,padX]`，非负整数；`orig_size + 2*padSize` 必须等于画布尺寸 |
| `ValidMaskHard` | 每张图像的二值有效测量掩膜 cell，尺寸与画布一致，且不能全空 |

每幅图像与其掩膜必须做相同的几何配准；填充区不能标记成真实测量。核心根据硬掩膜重新计算覆盖权重，不要求适配器提供软权重。`Z - focusDistance` 和 `Z + cumsum(distanceSteps)` 必须为正。

强度缩放应在输入适配阶段统一约定。当前 APRW 含振幅上限约束，不能将任意位深的整数像素值未经归一化直接当作强度输入。

```matlab
setup_project
inputData = syntheticInput();  % 可运行的格式示例
validateReconstructionInput(inputData);
save(fullfile(pwd, 'data', 'reconstruction_input.mat'), 'inputData', '-v7.3');
```

若 `data/` 尚不存在，应先创建。也可直接调用 `[fields, folder] = reconstructMultiPlane(inputData, options)`；公共文件入口会先进行输入检查。

## 彩色输入

先为三个通道分别准备灰度格式的 MAT，再保存通道清单：

```matlab
colorInput.channelFiles = {'channel_B.mat', 'channel_G.mat', 'channel_R.mat'};
save('color_input.mat', 'colorInput');
```

通道路径可为绝对路径，或相对于清单文件所在目录的路径。核心逐通道加载，避免同时保留所有通道测量堆栈。默认顺序 B/G/R 与 `channelNames/rgbOrder` 对应，通道名称应使用不同的简单目录名。

各通道的波长、间隔和标定从自身 `inputData` 读取。当前 RGB 配准使用平移模型，通道必须使用相同的像素尺寸及一致的横向坐标约定，不包含跨像素尺度的重采样。

彩色流程可在通道内部执行标量曝光均衡。若输入适配器已完成此处理，可设置 `config.inputBrightness.enabled = false`。

## 接入自己的采集流程

实现自己的导入函数，返回上述结构体或保存对应 MAT，即可使用公共入口。原始数据文件名、采集日期、网络路径和设备专用脚本不需要进入核心仓库。`examples/syntheticInput.m` 给出了完整的小规模结构示例。
