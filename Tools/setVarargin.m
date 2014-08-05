function opt = setVarargin( allfields, varargin )
% setVarargin( allfields, varargin )

opt   = struct();

nVarargs = length(varargin);
if mod(nVarargs,2)~=0
    error('Optional arguments should be pairs: "Property name", Value')
end

for k = 1:2:nVarargs
    opt.(varargin{k}) = varargin{k+1};
end

% Default values in options
for k = 1:size(allfields,1)
    if ~isfield(opt,allfields{k,1})
        opt.(allfields{k,1}) = allfields{k,2};
    end
end