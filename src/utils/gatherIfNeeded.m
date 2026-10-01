function data = gatherIfNeeded(data)
%GATHERIFNEEDED Bring a GPU array to the CPU; leave ordinary arrays unchanged.
% The CPU path does not require Parallel Computing Toolbox functions.
if isa(data, 'gpuArray'), data = gather(data); end
end
