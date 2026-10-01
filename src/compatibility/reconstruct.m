function varargout=reconstruct(varargin)
%RECONSTRUCT Legacy alias; new code calls reconstructPreparedData.
[varargout{1:nargout}]=reconstructPreparedData(varargin{:});
end
