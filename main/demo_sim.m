function [result,folder]=demo_sim(overrides,showFigures)
%DEMO_SIM 仿真主入口：修改下面参数，运行 demo_sim 或 demo_sim_fast。
% 批量调试可传入部分参数，例如 demo_sim(struct('iterations',120))。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setup_project;

%% 1. 光学真值：用于生成数据及已知真值参考组；长度单位为米
options.wavelength=514e-9;           % 波长，514 nm
options.pixelSize=3e-6;              % 探测器像素尺寸，3 um
options.sourceToSample=30e-3;        % 光源到样品距离，不是标定输出 Z
options.sourceOffsetPixels=[0,-0]; % 光源位置 [x,y]，像素；倾斜照明使物面覆盖移动
options.imageSize=1000;              % 相机窗口边长；超过图片短边时自动限制到图片短边
options.padding=320;                % 每侧补边；网格边长=imageSize+2*padding
% 当前求解器直接采样复场；过大倾角会触发采样检查，不应通过增大 M 强行扩展。
% 样品大小在 sample.sizePixels 中独立设置。
% 光源投影在传感器外时，轴向移动可使物面覆盖产生并集扩展；补边本身不是新增视场。

%% 2. 采集位置：样品到探测器的绝对距离，不是相邻帧移动量
options.sampleDistances=(0.5:1:12.5)*1e-3;      % 相机轴向移动，13 帧逐渐覆盖其它区域
options.calibrationDistances=(0.5:0.25:1.5)*1e-3; % 点阵短距离采集，避免离焦过大
% 两组首个位置须相同。标定 Z = sourceToSample + sampleDistances(1)。

%% 3. 标定：从独立点阵衍射图估计 M/N/Z
options.calibration.imageSize=192;        % 独立标定窗口，不随样品采集窗口改变
options.calibration.gridPitchPixels=28;   % 生成点阵的物面间距，像素
options.calibration.approxPitchPixels=29; % 寻峰用探测器近似间距，像素
options.calibration.backgroundSigma=30;  % 背景估计宽度，像素
options.calibration.projectionSmooth=5;  % 投影平滑宽度，像素

%% 4. 样品选择：mixed / resolution / cells / image / widefield
% mixed：原混合幅相样品；resolution：条纹靶；cells：细胞状纯相位样品。
% image：自定义图片。亮度映射为幅度或相位，具体规则见 createSimulationSample。
options.sample.type='image';
options.sample.amplitudeFile='cell.tif';  % image 模式：幅度图片；空表示幅度全为 1
options.sample.phaseFile='peppers.png';      % image 模式：相位图片；空表示相位为 0
options.sample.phaseScale=1;      % 相位乘数；image 模式对应最大相位 [rad]
options.sample.absorptionScale=1; % 吸收乘数；0 为纯相位
options.sample.sizePixels=512;     % 图片最长边，像素；与 imageSize 独立
% 图片范围之外为零；移动后越界的相机像素补零，并用 cameraMasks 标为无效采样。
options.sample.offsetPixels=[5,-3]; % 图片中心对齐首帧物面覆盖中心（按当前几何取整）
% widefield 为覆盖较大范围的幅相靶，便于观察首帧视场外的细节。

%% 5. 视场显示与客观评价
options.fov.displayUnion=true;     % 显示标定几何得到的联合视场；false 显示首帧视场
options.fov.saveAcquisitionGif=true; % 保存相机逐帧覆盖大样品的动画
options.evaluation.phaseAmplitudeThreshold=0.05; % 真值幅度低于此值，不评价相位
% 第四组：直接调用实验核心 reconstructMultiPlane，与标定球面波简化算法比较。
options.core.enabled=true;
options.core.tvEnabled=true;
options.core.adaptiveConstraintEnabled=true;
options.core.adaptiveConstraintStrength=0.05;
options.core.adaptiveConstraintEdgeWidth=32;
% 共用 iterations；聚焦固定在仿真首个距离，避免混入自动聚焦误差。
% 幅度 PSNR 峰值固定为 1；相位 PSNR 使用圆周误差及固定范围 2*pi。
% SSIM 仅统计窗口完全位于评价区域内的像素；窄小区域无有效窗口时为 NaN。

%% 6. 求解与输出：四组共用数据及迭代次数，第四组调用实验核心
options.iterations=20;
options.recordEvery=options.iterations;           % 显示和误差记录间隔
options.noiseStd=0;               % 样品强度图加性高斯噪声标准差
options.seed=7;
options.outputRoot=fullfile(root,'outputs','simulation');

%% 执行：日常使用只修改上面的参数
if nargin>=1 && ~isempty(overrides)
    options=applyOverrides(options,overrides,'options');
end
if nargin<2,showFigures=true;end
validateattributes(showFigures,{'logical'},{'scalar'});
options.showFigures=showFigures;
[result,folder]=simulateModelComparison(options);
end

function base=applyOverrides(base,changes,path)
% 递归覆盖指定字段；参数名拼错时立即报错。
assert(isstruct(changes)&&isscalar(changes),'Simulation overrides must be a scalar struct.');
names=fieldnames(changes);
for k=1:numel(names)
    name=names{k};
    assert(isfield(base,name),'Unknown simulation option: %s.%s',path,name);
    if isstruct(base.(name))
        base.(name)=applyOverrides(base.(name),changes.(name),[path '.' name]);
    else
        base.(name)=changes.(name);
    end
end
end
