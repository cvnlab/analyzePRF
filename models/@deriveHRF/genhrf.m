function genhrf(obj)
% Creates and stores an HRF using a double-gamma model
%
% Syntax:
%   obj.genhrf
%
% Description:
%   Creates an HRF kernel for convolution. The stored kernel has unit area
%   so that it preserves signal area after convolution. The kernel is
%   specified in a time x 1 vector orientation.
%
%   Typical values for the HRF parameters (which are in units of seconds)
%   are [6 12 10 20].
%
% Inputs:
%   none
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   none
%


% Obj variables
tr = obj.tr;
hrfParams = obj.hrfParams;

% Unpack hrfParams
gamma1 = hrfParams(1);
gamma2 = hrfParams(2);
gammaScale = hrfParams(3);
duration = hrfParams(4);

% Defin the timebase
timebase = 0:tr:ceil(duration/tr);

% Create the double gamma function
hrf = gampdf(timebase,gamma1, 1) - ...
    gampdf(timebase, gamma2, 1)/gammaScale;

% Set to zero at onset
hrf = hrf - hrf(1);

% Normalize the kernel to have unit area
hrf = hrf/sum(abs(hrf));

% Store the hrf in the object. Transpose the vector so that it is time x 1
obj.hrf = hrf';

end

