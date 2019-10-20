function results = analyzePRF(stimulus,data,varargin)


%% input parser
p = inputParser; p.KeepUnmatched = false;

% Required
p.addRequired('stimulus',@(x)(iscell(x) || ismatrix(x)));
p.addRequired('data',@(x)(iscell(x) || ismatrix(x)));

p.addParameter('modelClass','pRF_timeShift',@ischar);
p.addParameter('modelOpts',{'typicalGain',30},@iscell);
p.addParameter('modelPayload',{},@iscell);
p.addParameter('tr',0.8,@isscalar);
p.addParameter('vxs',[],@isvector);
p.addParameter('hrf',[],@isvector);
p.addParameter('maxIter',500,@isscalar);
p.addParameter('verbose',true,@islogical);

% parse
p.parse(stimulus,data, varargin{:})

verbose = p.Results.verbose;

%% Alert the user
if verbose
    fprintf(['Fitting the ' p.Results.modelClass ' model.\n\n']);
end


%% Massage inputs and set constants

% Place the stimulus and data in a cell if not already so
if ~iscell(stimulus)
    stimulus = {stimulus};
end
if ~iscell(data)
    data = {data};
end

% Identify the row and columns of the data matrix
dimdata = 1;
totalVxs = size(data{1},dimdata);

% Define vxs (the voxel/vertex set to process)
if isempty(p.Results.vxs)
    vxs = 1:size(data{1},totalVxs);
else
    vxs = p.Results.vxs;
end


%% Stimulus prep

% Obtain the dimensions of the stimulus frames
res = [size(stimulus{1},1) size(stimulus{1},2)];

% Vectorize the stimuli. Also, add a dummy column to indicate run breaks
for ii=1:length(stimulus)
    stimulus{ii} = squish(stimulus{ii},2)';
    stimulus{ii} = [stimulus{ii} ii*ones(size(stimulus{ii},1),1)];
end


%% HRF prep
if isempty(p.Results.hrf)
    hrf = getcanonicalhrf(p.Results.tr,p.Results.tr)';
else
    hrf = p.Results.hrf;
end


%% Set up model

% Create the model object
model = feval(p.Results.modelClass,data,stimulus,res,hrf,p.Results.tr,...
    'payload',p.Results.modelPayload, ...
    p.Results.modelOpts{:});

% Set model verbosity
model.verbose = verbose;

% Prep the raw data
data = model.prep(data);

% Generate seeds
seeds = model.seeds(data,vxs);


%% Fit the data

% Options for lsqcurvefit
options = optimset('Display','off','FunValCheck','on', ...
    'MaxFunEvals',Inf,'MaxIter',p.Results.maxIter, ...
    'TolFun',1e-6,'TolX',1e-6, ...
    'OutputFcn',@(a,b,c) outputfcnsanitycheck(a,b,c,1e-6,10));
options.Algorithm = 'levenberg-marquardt';

% Alert the user
if verbose
    tic
    fprintf(['Fitting non-linear model over ' num2str(length(vxs)) ' vertices:\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% Loop through the voxels/vertices in vxs
parfor ii=1:length(vxs)
    
    % Update progress bar
    if verbose && mod(ii,round(length(vxs)/50))==0
        fprintf('\b.\n');
    end
    
    data2 = cellfun(@(x) x(vxs(ii),:),data,'UniformOutput',0);
    data3 = @(vxs) cellfun(@(x) subscript(squish(x,dimdata),{vxs ':'})',data2,'UniformOutput',0);
    data4 = catcell(1,data3(1));
    data5 = model.clean(data4);
    
    % Loop over seed sets
    seedParams = nan(length(seeds),model.nParams);
    seedMetric = nan(length(seeds),1);
    
    for ss = 1:length(seeds)
        seed = seeds{ss}(vxs(ii),:);
        x0 = seed;
        
        % Loop over model stages
        for bb = 1:model.nStages
            x0 = lsqcurvefit(...
                @(x,y) model.forward([x x0(model.fixSet{bb})]),...
                x0(model.floatSet{bb}),...
                [],...
                data5,[],[],options);
            x0 = [x0 seed(model.fixSet{bb})];
        end
        
        seedParams(ss,:) = x0;
        seedMetric(ss) = model.metric(data5,x0);
    end
    
    [~,bestSeedIdx]=max(seedMetric);
    
    parParams(ii,:) = seedParams(bestSeedIdx,:);
    parMetric(ii) = seedMetric(bestSeedIdx);
    
end

% report completion of loop
if verbose
    toc
    fprintf('\n');
end

% Map the par variables into full variables
params = nan(totalVxs,model.nParams);
params(vxs,:) = parParams;
clear parParams

metric = nan(totalVxs,1);
metric(vxs) = parMetric;
clear parMetric;


%% Obtain the results
results = model.results(params, metric);

% Add the model information
results.model.class = p.Results.modelClass;
results.model.inputs = {stimulus, res, hrf p.Results.tr};
results.model.opts =  p.Results.modelOpts;
results.model.payload =  p.Results.modelPayload;

% Store the calling options
results.meta.vxs = p.Results.vxs;
results.meta.tr = p.Results.tr;
results.meta.maxIter = p.Results.maxIter;

end
