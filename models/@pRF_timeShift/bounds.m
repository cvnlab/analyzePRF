function [lb, ub] = bounds(obj)
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

res = obj.res;
nParams = obj.nParams;

% Define outputs
lb = nan(1,nParams);
ub = nan(1,nParams);

% The lower bounds
lb(1) = -res(1)/2;      % xPosition
lb(2) = -res(2)/2;      % yPosition
lb(3) = 0;              % sigma
lb(4) = 0;              % gain (amplitude) of response
lb(5) = 0;              % compressive exponent
lb(6) = -3;             % HRF temporal shift (units of TRs)

% The upper bounds
ub(1) = 1.5*res(1);     % xPosition
ub(2) = 1.5*res(2);     % yPosition
ub(3) = max(res);            % sigma
ub(4) = 1000;           % gain (amplitude) of response
ub(5) = 3;              % compressive exponent
ub(6) = 3;              % HRF temporal shift (units of TRs)

end

