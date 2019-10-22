function data = prep(obj,data)
% Initial pre-processing of data cell array prior to model fitting
%
% Syntax:
%   data = obj.prep(data)
%
% Description:
%   Implements an initial sanity cleaning of the input cell array of data
%   to remove "bad" voxels.
%
% Inputs:
%   data                  - A cell array of [v t]  matricies. The fMRI
%                           time-series data across t TRs, for v vertices /
%                           voxels.
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   data                  - A cell array of [v t]  matricies. The fMRI
%                           time-series data across t TRs, for v vertices /
%                           voxels.
%


% Obj variables
nAcqs = obj.nAcqs;
dimdata = obj.dimdata;
dimtime = obj.dimtime;

% Set "bad" voxels to have uniform zero values
bad = cellfun(@(x) any(~isfinite(x),dimtime) | all(x==0,dimtime),data,'UniformOutput',0);  % if non-finite or all 0
bad = any(cat(dimtime,bad{:}),dimtime);  % badness in ANY run
for ii=1:nAcqs
    data{ii}(repmat(bad,[ones(1,dimdata) size(data{ii},dimtime)])) = 0;
end


end
