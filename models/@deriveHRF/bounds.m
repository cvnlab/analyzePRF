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
lb(1) = -1;              % Aratio
lb(2) = 2;              % alpha1
lb(3) = 0.5;            % beta1
lb(4) = 6;              % alpha2
lb(5) = 0;              % beta2                    
lb(6) = 0;              % gain (amplitude) of response

% The upper bounds
ub(1) = 1;              % Aratio
ub(2) = 10;             % alpha1
ub(3) = 2;              % beta1
ub(4) = 25;             % alpha2
ub(5) = 1.5;            % beta2
ub(6) = Inf;            % gain (amplitude) of response

end

