# 核心回归测试

```matlab
setup_project
addpath(fullfile(pwd, 'tests'))
run_core_tests
```

所有公共测试只使用合成数据，不读取个人输入、不调用 `local/` 工具。输出写入被 Git 忽略的 `outputs/`。

| 测试 | 检查内容 |
| --- | --- |
| `verify_frft` | 整数阶的精确变换关系、非整数阶的高斯幅值近似不变性 |
| `verify_autofocus` | 固定距离、平坦曲线、非法搜索步长、ADFrFT 公式与 ROI、方法选择 |
| `verify_geometry_output` | 几何边界、三种 APRW 输出模式、无效像素零填充 |
| `verify_project_entrypoints` | 标准/快速灰度、聚焦数据与 PNG、绘图坐标、彩色融合及缓存复用、无窗口运行后的图窗清理 |

这些测试验证实现和接口的一致性，不替代真实光学实验的重建质量评估。

支持图形窗口的 MATLAB 会话还可运行 `verify_manual_focus_roi`。测试自动操作选区确认和取消，检查跨迭代复用、坐标保存以及结果回看中的选区轮廓；该交互测试不包含在默认无交互测试集中。
