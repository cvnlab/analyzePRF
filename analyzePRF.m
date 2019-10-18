function results = analyzePRF(stimulus,data,varargin)


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('stimulus',@(x)(iscell(x) || ismatrix(x)));
p.addRequired('data',@(x)(iscell(x) || ismatrix(x)));

p.addParameter('modelClass','prf_timeShift',@ischar);
p.addParameter('tr',0.8,@isscalar);
p.addParameter('vxs',0.8,@isvector);
p.addParameter('hrf',[],@isvector);
p.addParameter('maxiter',500,@isscalar);
p.addParameter('typicalGain',30,@isscalar);


% parse
p.parse(stimulus,data, varargin{:})


% Place the stimulus and data in a cell if not already so
if ~iscell(stimulus)
    stimulus = {stimulus};
end
if ~iscell(data)
    data = {data};
end
vxs = p.Results.vxs;


dimdata = 1;
dimtime = 2;
xyzsize = size(data{1},1);

numvxs = prod(xyzsize);

% calc
res = sizefull(stimulus{1},2);
numruns = length(data);


% Vectorize the stimuli. Also, add a dummy column to indicate run breaks
for ii=1:length(stimulus)
    stimulus{ii} = squish(stimulus{ii},2)';
    stimulus{ii} = [stimulus{ii} ii*ones(size(stimulus{ii},1),1)]; 
end

% deal with data badness (set bad voxels to be always all 0)
bad = cellfun(@(x) any(~isfinite(x),dimtime) | all(x==0,dimtime),data,'UniformOutput',0);  % if non-finite or all 0
bad = any(cat(dimtime,bad{:}),dimtime);  % badness in ANY run
for ii=1:numruns
    data{ii}(repmat(bad,[ones(1,dimdata) size(data{ii},dimtime)])) = 0;
end

% what HRF should we use?
if isempty(p.Results.hrf)
    hrf = getcanonicalhrf(p.Results.tr,p.Results.tr)';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE MODEL

% Create the model object
modelObj = pRF_timeShift(stimulus,res,hrf);

% Assign the typical gain parameter
modelObj.gain = p.Results.typicalGain;

% define the model with a call out to the external function modelCore
modelfun = @(pp) modelObj.forward(pp);

% Model bounds. The NaN in the lb indicates that these parameters are fixed
% Parameters are:
%   xPos, yPos, sigma, amplitude, gain, hrfShift
modelObj.fixed = [5 6];
[lb,ub] = modelObj.bounds();

% Initial model, with the compressive non-linearity fixed
M1 = {[] [lb; ub] modelfun};

% Second model, which takes the params from the first stage to generate a
% seed for model fitting, allowing the compressive non-linearity parameter
% to vary
modelObj.fixed = [];
[lb,ub] = modelObj.bounds();
M2 = {@(ss)ss [lb; ub] @(ss)modelfun};

% Chain the sub-models into the full model variable
model = {M1 M2};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE SEEDS

% init
seeds = [];

% generic large seed
    modelObj.seedScale = 'large';
    seeds = [ seeds; modelObj.initial() ];

% generic small seed
    modelObj.seedScale = 'small';
    seeds = [ seeds; modelObj.initial() ];

% super-grid seed
    [supergridseeds,rvalues] = analyzePRFcomputesupergridseeds(res,stimulus,data,modelfun, ...
        1,dimdata,dimtime, ...
        p.Results.typicalGain,[],modelObj);

% make a function that individualizes the seeds
if exist('supergridseeds','var')
    seedfun = @(vx) [[seeds];
        [subscript(squish(supergridseeds,dimdata),{vx ':'})]];
else
    seedfun = @(vx) [seeds];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE STIMULUS AND DATA




data = @(vxs) cellfun(@(x) subscript(squish(x,dimdata),{vxs ':'})',data,'UniformOutput',0);


options = optimset('Display','off','FunValCheck','on', ...
    'MaxFunEvals',Inf,'MaxIter',pp.Results.maxiter, ...
    'TolFun',1e-6,'TolX',1e-6, ...
    'OutputFcn',@(a,b,c) outputfcnsanitycheck(a,b,c,1e-6,10));
options.Algorithm = 'levenberg-marquardt';



for ii=1:vnum
               params0 = ...
                lsqcurvefit(@(x,y) double(modelObj.forward(x) / datastd),seed,[],double(data),[],[],options);

    
end


