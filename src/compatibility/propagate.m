function field = propagate(varargin)
%PROPAGATE Compatibility entry; new code calls propagateAngularSpectrum.
field = propagateAngularSpectrum(varargin{:});
end
