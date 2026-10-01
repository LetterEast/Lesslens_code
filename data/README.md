# 本地输入

此目录除本说明外全部由 Git 忽略。公共代码读取已经准备好的 MAT 文件，格式见 [输入接口](../docs/INPUT_FORMAT.md)。也可用配置指向项目外的文件。

- 灰度：默认 `reconstruction_input.mat`，包含 `inputData`。
- 彩色：默认 `color_input.mat`，包含 `colorInput.channelFiles`，列出三个通道 MAT。
- 原始图像、标定和历史文件可在这里分目录存放，不属于公共示例。

数据准备工具由使用者自行实现，或作为个人扩展放在 `local/`。
