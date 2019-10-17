function [lb, ub] = modelBounds(res,fixed)
% Bounds for the forward model
%
% Inputs:
%   res                   - Scalar. the resolution of the stimulus in
%                           pixels (typically ~100).
%   fixed                 - Vector, optional. Indicates the indices in the
%                           lower bound to be set to NaN, thus indicating
%                           to the calling function that these parameters
%                           are to be fixed at their x0 value.
%                           

% Define outputs
nParams = 7;
lb = nan(1,nParams);
ub = nan(1,nParams);

% Handle incomplete arguments
if nargin==0
    return
end
if nargin==1
    fixed = [];
end

% The lower bounds
lb(1) = 1-res(1)+1;     % xPosition
lb(2) = 1-res(2)+1;     % yPosition
lb(3) = 0;              % sigma
lb(4) = 0;              % gain (amplitude) of response
lb(5) = 0;              % compressive exponent
lb(6) = -2;             % HRF temporal shift (units of TRs)
lb(7) = NaN;            % Unused, so always fixed

% The upper bounds
ub(1) = 2*res(1)-1;     % xPosition
ub(2) = 2*res(2)-1;     % yPosition
ub(3) = Inf;            % sigma
ub(4) = Inf;            % gain (amplitude) of response
ub(5) = Inf;            % compressive exponent
ub(6) = 2;              % HRF temporal shift (units of TRs)
ub(7) = NaN;            % Unused, so always fixed

% Fix the specified params
lb(fixed) = NaN;

end

