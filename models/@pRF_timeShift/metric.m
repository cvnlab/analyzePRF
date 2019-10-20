function metric = metric(obj, signal, x)
% Evaluates the match between a signal and a model fit
%
% Syntax:
%   metric = obj.metric(signal, x)
%
% Description:
%   Given a time series signal and the parameters of the forward model,
%   returns a metric that describes how well the two match.
%
% Inputs:
%   signal                - 1 x time vector. The data to be fit.
%   x                     - 1 x nParams vector of parameter values.
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   metric                - Scalar.
%

% Implement an R^2 metric
metric = calccorrelation(signal, obj.forward(x))^2;

end

