function [lb, ub] = bounds(obj)
% Returns bounds on the model parameters
%
% Syntax:
%   [lb, ub] = obj.bounds()
%
% Description:
%   Bounds for the deiveHRF model.
%
% Inputs:
%   none
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   lb, ub                - 1 x nParams vectors.
%


% Obj variables
nParams = obj.nParams;

% Define outputs
lb = nan(1,nParams);
ub = nan(1,nParams);

% The lower bounds
lb(1) = 2;              % gamma1
lb(2) = 1.1;            % gamma2 = gamma1 * x(2)
lb(3) = 1;              % ratio of gamma 1 / gamma 2 amplitudes
lb(4) = 0;              % gain (amplitude) of response

% The upper bounds
ub(1) = 10;              % gamma1
ub(2) = 3;              % gamma2 = gamma1 * x(2)
ub(3) = 200;             % ratio of gamma 1 / gamma 2 amplitudes
ub(4) = Inf;            % gain (amplitude) of response

end

