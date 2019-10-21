function signal = clean(obj,signal)
% Cleans a passed time-series in preparation for non-linear fitting
%
% Syntax:
%   signal = obj.clean(signal)
%
% Description:
%   The passed signal vector is cleaned and returned. In this
%   implementation, nuisance effects as modeled in the projection matrix T
%   are removed from the data. This includes removal of the zeroeth
%   frequency and thus resulting vector is mean centered.
%
% Inputs:
%   signal                - 1 x time vector.
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   signal                - 1 x time vector.
%

signal = double(obj.T*signal);

end
