function outputRoot = reconstruct_color_fast(config)
%RECONSTRUCT_COLOR_FAST 彩色快速入口：相同算法与参数，关闭过程窗口。
% 修改根目录 color_config.m 后运行：folder = reconstruct_color_fast;
% 自动准备输入的行为由 prepareInput 决定；彩色和各通道灰度结果照常保存。
% 不导出聚焦曲线图片，仍保存聚焦数值和各次记录。手动聚焦选区仍需交互。
% 可传入完整 config 用于批处理；不修改迭代次数或约束。
projectRoot=fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);setup_project;
if nargin<1,config=[];end
outputRoot=reconstruct_color(config,false);
end
