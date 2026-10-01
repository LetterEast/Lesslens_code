function [result,folder]=demo_sim_fast(overrides)
%DEMO_SIM_FAST 复用 demo_sim 顶部参数；关闭窗口，仍保存结果和对比图。
if nargin<1,overrides=[];end
[result,folder]=demo_sim(overrides,false);
end
