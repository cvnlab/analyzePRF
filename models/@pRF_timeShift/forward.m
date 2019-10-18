function signal = forward(obj, pp)
% Forward model for the pRF search

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

% Combine the stimulus cells into a single vector
stimulus = catcell(1,stimulus);

% Dot product of the stimulus by the Gaussian window (the neural signal),
% subjected to a compressive non-linearity by raising to the pp(5)
% exponent, and scaled in amplitude by pp(4).
neuralSignal = posrect(pp(4)) * (stimulus*gaussVector) .^ posrect(pp(5));

% Define the acquisition groupings
acqGroups = stimulus(:,prod(res)+1);

% Shift the hrf by the number of TRs specified in pp(6)
hrf = fshift(hrf,pp(6));

% Convolve the neural signal by the passed hrf, respecting the boundaries
% betwen the acquisitions
signal = conv2run(neuralSignal,hrf,acqGroups);

end

