function seeds = seeds(obj,data,vxs)
% Generate parameter seeds for the non-linear search
%
% Syntax:
%   seeds = obj.seeds(data,vxs)
%
% Description:
%   Generates a set of seed parameters for each voxel/vertex in vxs.
%
% Inputs:
%   data                  - A matrix [v t] or cell array of such
%                           matricies. The fMRI time-series data across t
%                           TRs, for v vertices / voxels. The data should
%                           have bassed through the prep stage.
%   vxs                   - Vector. A list of vertices/voxels to be
%                           processed.
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   seeds                 - Cell array. Each cell contains a matrix of
%                           [v nParams] and is one of the seed sets.
%

% Derived vars
totalVxs = size(data{1},1);

% Generate seeds, which is just the initial values
x0 = obj.initial;
seeds{1} = repmat(x0,totalVxs,1);

end


