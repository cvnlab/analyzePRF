function fit = forward(obj, pp)
% Forward model for the pRF search
%
% Syntax:
%   fit = obj.forward(pp)
%
% Description:
%   Returns a time-series vector that is the predicted response to a 2D
%   stimulus, based upon the parameters provided in pp.
%
% Inputs:
%   pp                    - 1xnParams vector.
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   fit                   - 1xtime vector.
%


% Obj variables
stimulus = obj.stimulus;
acqGroups = obj.acqGroups;
res = obj.res;
hrf = obj.hrf;
tr = obj.tr;
xx = obj.xx;
yy = obj.yy;
resmx=max(res);

% Gaussian at [x, y] pp(1), pp(2), of sigma size pp(3)
gaussWindow = makegaussian2d(resmx,pp(1),pp(2),pp(3),pp(3),xx,yy,0,0);

% Normalization scalar
gaussNorm = (2*pi*abs(pp(3))^2);

% Gaussian window normalized, cropped to <res>, and vectorized
gaussVector =  vflatten(placematrix(zeros(res), gaussWindow / gaussNorm));

% Dot product of the stimulus by the Gaussian window (the neural signal),
% subjected to a compressive non-linearity by raising to the pp(5)
% exponent.
neuralSignal =  pp(4) * (stimulus*gaussVector).^ pp(5);

% Shift the hrf by the number of seconds specified in pp(6)
hrf = fshift(hrf,pp(6)/tr);

% Convolve the neural signal by the passed hrf, respecting the boundaries
% betwen the acquisitions
fit = conv2run(neuralSignal,hrf,acqGroups);

% Partial the data to remove the effects that are represented in the
% regression matrix T
fit = obj.T*fit;


end

