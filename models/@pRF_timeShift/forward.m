function fit = forward(obj, pp)
% Forward model for the pRF search

% Obj variables
stimulus = obj.stimulus;
res = obj.res;
hrf = obj.hrf;
xx = obj.xx;
yy = obj.yy;
resmx=max(res);

% Gaussian at [x, y] pp(1), pp(2), of sigma size pp(3)
gaussWindow = makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0);

% Normalization scalar
gaussNorm = (2*pi*abs(pp(3))^2);

% Gaussian window normalized, cropped to <res>, and vectorized
gaussVector =  [vflatten(placematrix(zeros(res), gaussWindow / gaussNorm)); 0];

% Dot product of the stimulus by the Gaussian window (the neural signal),
% subjected to a compressive non-linearity by raising to the pp(5)
% exponent.
neuralSignal = (stimulus*gaussVector).^ posrect(pp(5));

% Scale the neural signal to have unit amplitude and then apply the gain
neuralSignal = neuralSignal - min(neuralSignal);
maxSignal = max(neuralSignal);
if maxSignal~=0
    neuralSignal = (neuralSignal / maxSignal) * posrect(pp(4));
end

% Define the acquisition groupings
acqGroups = stimulus(:,prod(res)+1);

% Shift the hrf by the number of TRs specified in pp(6)
hrf = fshift(hrf,pp(6));

% Convolve the neural signal by the passed hrf, respecting the boundaries
% betwen the acquisitions
fit = conv2run(neuralSignal,hrf,acqGroups);

% Partial the data to remove the effects that are represented in the
% regression matrix T
fit = obj.T*fit;


end

