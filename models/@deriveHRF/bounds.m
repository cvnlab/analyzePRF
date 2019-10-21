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
lb(2) = 3;              % gamma2
lb(3) = 1;              % gammaScale
lb(4) = -Inf;           % gain (amplitude) of response

% The upper bounds
ub(1) = 8;              % gamma1
ub(2) = 10;             % gamma2
ub(3) = 20;             % gammaScale
ub(4) = Inf;            % gain (amplitude) of response

end

