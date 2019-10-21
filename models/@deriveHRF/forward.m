function [fit, hrf] = forward(obj, pp)
% Forward model for the derive hrf search
%
% Syntax:
%   fit = obj.forward(pp)
%
% Description:
%   Returns a time-series vector that is a square-wave vector of neural
%   activity as defined by the stimulus, subject to convolution by an HRF
%   that is defined by the params.
%
% Inputs:
%   pp                    - [1 nParams] vector.
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   fit                   - [nTRs 1] vector.
%   hrf                   - [duration 1] vector.
%


% Force the parameters within bounds
if obj.forceBounds
    [lb, ub] = obj.bounds;
    idx = pp < lb;
    pp(idx) = lb(idx);
    idx = pp > ub;
    pp(idx) = ub(idx);
end

% Unpack the params
gamma1 = posrect(pp(1));
gammaTimeRatio = posrect(pp(2));
gamma2 = gamma1 * gammaTimeRatio;
gammaScale = posrect(pp(3));
gain = posrect(pp(4));

% Obj variables
stimulus = obj.stimulus;
acqGroups = obj.acqGroups;
tr = obj.tr;
duration = obj.duration;

% THe neural signal is the stimulus, scaled by the gain.
neuralSignal =  posrect(gain) * stimulus;

% Define the timebase in TRs
timebase = 0:tr:ceil(duration/tr);

% Create the double gamma function
hrf = gampdf(timebase, gamma1, 1) - ...
    gampdf(timebase, gamma2, 1)/gammaScale;

% Set to zero at onset
hrf = hrf - hrf(1);

% Normalize the kernel to have unit area. Transpose the vector so that it
% is time x 1
hrf = (hrf/sum(abs(hrf)))';

% Convolve the neural signal by the passed hrf, respecting the boundaries
% betwen the acquisitions
fit = conv2run(neuralSignal,hrf,acqGroups);

% Partial the data to remove the effects that are represented in the
% regression matrix T
fit = obj.T*fit;


end

