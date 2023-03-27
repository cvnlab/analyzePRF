%% Example 2: Inspect the model fits obtained by analyzePRF

%% Download dataset (if necessary) and add analyzePRF to the MATLAB path

setup;

%% Load and analyze the data (this follows the example1.m script)

% Start parallel MATLAB to speed up execution
if isempty(gcp)
  parpool;
end

% Load in the data
load('exampledataset.mat');

% Upsample the data to match the stimulus rate
data = tseriesinterp(data,2,1,2);

% Analyze the data
results = analyzePRF(stimulus,data,1,struct('seedmode',[0 1],'display','off'));
%%

%% Perform some setup

% Define some variables
res = [100 100];                    % row x column resolution of the stimuli
resmx = 100;                        % maximum resolution (along any dimension)
hrf = results.options.hrf;          % HRF that was used in the model
degs = results.options.maxpolydeg;  % vector of maximum polynomial degrees used in the model

% Pre-compute cache for faster execution
[d,xx,yy] = makegaussian2d(resmx,2,2,2,2);

% Prepare the stimuli for use in the model
stimulusPP = {};
for p=1:length(stimulus)
  stimulusPP{p} = squish(stimulus{p},2)';  % this flattens the image so that the dimensionality is now frames x pixels
  stimulusPP{p} = [stimulusPP{p} p*ones(size(stimulusPP{p},1),1)];  % this adds a dummy column to indicate run breaks
end

% Define the model function.  This function takes parameters and stimuli as input and
% returns a predicted time-series as output.  Specifically, the variable <pp> is a vector
% of parameter values (1 x 5) and the variable <dd> is a matrix with the stimuli (frames x pixels).
% Although it looks complex, what the function does is pretty straightforward: construct a
% 2D Gaussian, crop it to <res>, compute the dot-product between the stimuli and the
% Gaussian, raise the result to an exponent, and then convolve the result with the HRF,
% taking care to not bleed over run boundaries.
modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));

% %%%%% BRIEF ASIDE:
% 
% Note that the above line of code is highly compressed and is optimized for speed,
% not necessarily readability. Here we show how the line can be broken up, omitting
% some of the less critical steps:
% 
% % here, we have a set of parameters, like those present in the
% % output from analyzePRF. pp is just a 5-element vector with
% % [ROW COL SIGMA GAIN EXPONENT].
% pp = results.params(1,:,1);
% 
% % make a Gaussian and give it a specific scale
% mm = makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0);
% mm = mm / (2*pi*abs(pp(3))^2);
% 
% % prepare the stimulus (here, we take just the first run of stimulus)
% dd = stimulus{1};  % e.g. 100 x 100 x 300
% dd = squish(dd,2);  % e.g. 10000 x 300
% 
% % compute the neural drive (dot product of stimulus with Gaussian,
% % raise to an exponent, apply gain factor)
% neural = posrect(pp(4))  *  (  (dd' * mm(:)) .^ posrect(pp(5)) );  % time x 1
% 
% % convolve with HRF
% ts = conv(neural,hrf);
% ts = ts(1:length(neural));  % crop to original length
% 
% %%%%%

% Construct projection matrices that fit and remove the polynomials.
% Note that a separate projection matrix is constructed for each run.
polymatrix = {};
for p=1:length(degs)
  polymatrix{p} = projectionmatrix(constructpolynomialmatrix(size(data{p},2),0:degs(p)));
end

%% Inspect the data and the model fit

% Which voxel should we inspect?  Let's inspect the second voxel.
vx = 2;

% For each run, collect the data and the model fit.  We project out polynomials
% from both the data and the model fit.  This deals with the problem of
% slow trends in the data.
datats = {};
modelts = {};
for p=1:length(data)
  datats{p} =  polymatrix{p}*data{p}(vx,:)';
  modelts{p} = polymatrix{p}*modelfun(results.params(1,:,vx),stimulusPP{p});
end

% Visualize the results
figure; hold on;
set(gcf,'Units','points','Position',[100 100 1000 100]);
plot(cat(1,datats{:}),'r-');
plot(cat(1,modelts{:}),'b-');
straightline(300*(1:4)+.5,'v','g-');
xlabel('Time (s)');
ylabel('BOLD signal');
ax = axis;
axis([.5 1200+.5 ax(3:4)]);
title('Time-series data');
