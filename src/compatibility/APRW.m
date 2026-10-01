function varargout = APRW(varargin)
%APRW Compatibility entry; new code calls reconstructMultiPlane.
[varargout{1:nargout}] = reconstructMultiPlane(varargin{:});
end
