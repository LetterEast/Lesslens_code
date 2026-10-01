
function [sub_pks,sub_locs] = findProjectionPeaks(proj,Pitch)
% 1. 基础寻峰
[pks,locs] = findpeaks(proj, 'MinPeakDistance', Pitch*0.8, 'MinPeakHeight', 0.3*max(proj));

if ~isempty(locs)
    % --- 修正点 1：计算整体波幅作为参考基准 ---
    amplitude = max(proj) - min(proj);
    threshold = amplitude * 0.2;

    % --- 修正点 2：左边缘过滤 ---
    diff_left = pks(1) - proj(1);
    if locs(1) < Pitch*0.6 && diff_left < threshold
        pks(1) = [];
        locs(1) = [];
    end

    % --- 修正点 3：右边缘过滤 ---
    if ~isempty(locs)
        diff_right = pks(end) - proj(end);
        if locs(end) > (numel(proj) - Pitch*0.6) && diff_right < threshold
            pks(end) = [];
            locs(end) = [];
        end
    end
end

% 2. 边界安全检查
valid_idx = locs > 1 & locs < numel(proj);
locs = locs(valid_idx);
pks = proj(locs);

if isempty(locs)
    sub_pks = []; sub_locs = []; return;
end

% 3. 三点抛物线插值 (初步亚像素定位)
p_m1 = proj(locs - 1);
p_0  = pks;
p_p1 = proj(locs + 1);
denominator = (p_m1 - 2.*p_0 + p_p1);

denominator(abs(denominator)<eps) = eps;
delta_locs = max(-0.5,min(0.5,0.5 .* (p_m1 - p_p1) ./ denominator));
temp_sub_locs = locs + delta_locs; % 临时存储初步位置

% ============================ 最高精度改进：局部质心精细化 ============================
% 在抛物线插值的基础上，利用局部灰度信息进行重心校正，消除干扰线产生的偏移
sub_locs = zeros(size(temp_sub_locs));
sub_pks  = zeros(size(temp_sub_locs));

half_win = round(Pitch / 4); % 定义局部质心窗口大小
for i = 1:length(temp_sub_locs)
    % 确定当前点周围的局部索引 ROI
    cur_loc = temp_sub_locs(i);
    roi_idx = round(cur_loc - half_win) : round(cur_loc + half_win);
    roi_idx = roi_idx(roi_idx > 0 & roi_idx <= numel(proj)); % 边界检查

    % 提取局部投影数据并【动态去背景】
    % 减去局部最小值能极大程度消除你图中那两条"干扰线"带来的基数抬升
    local_proj = proj(roi_idx);
    local_proj = local_proj - min(local_proj);

    % 计算灰度加权质心 (精度远高于单纯的抛物线拟合)
    if sum(local_proj) > 0
        sub_locs(i) = sum(roi_idx(:) .* local_proj(:)) / sum(local_proj);
        % 对应的峰值高度沿用插值结果或取局部最大值
        sub_pks(i) = max(proj(roi_idx));
    else
        sub_locs(i) = temp_sub_locs(i);
        sub_pks(i) = p_0(i);
    end
end
% ====================================================================================
end
