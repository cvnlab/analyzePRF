function results = analyzePRF(stimulus,data,tr,varargin)
% Non-linear model fitting of fMRI time-series data
%
% Syntax:
%  results = analyzePRF(stimulus,data)
%
% Description:
%   Lorem ipsum
%
% Inputs:
%   stimulus              - A matrix [x y t] or cell array of such
%                           matrices. 
%   data                  - A matrix [v t] or cell array of such
%                           matricies. The fMRI time-series data across t
%                           TRs, for v vertices / voxels.
%   tr                    - Scalar. The TR of the fMRI data in seconds.
%
% Optional key/value pairs:
%  'modelClass'           - Char vector. The name of one of the available
%                           model objects. Choices include:
%                             {'pRF','pRF_timeShift'}
%  'modelOpts'            - A cell array of key-value pairs that are passed
%                           to the model object at that time of object
%                           creation.
%  'modelPayload'         - A cell array of additional inputs that is
%                           passed to the model object. The form of the
%                           payload is defined by the model object.
%  'vxs'                  - Vector. A list of vertices/voxels to be
%                           processed.
%  'maxIter'              - Scalar. The maximum number of iterations
%                           conducted by lsqcurvefit in model fitting.
%  'verbose'              - Logical.
%
% Outputs:
%   results               - Structure
%


%% input parser
p = inputParser; p.KeepUnmatched = false;

% Required
p.addRequired('stimulus',@(x)(iscell(x) || ismatrix(x)));
p.addRequired('data',@(x)(iscell(x) || ismatrix(x)));
p.addRequired('tr',@isscalar);

p.addParameter('modelClass','pRF_timeShift',@ischar);
p.addParameter('modelOpts',{'typicalGain',30},@iscell);
p.addParameter('modelPayload',{},@iscell);
p.addParameter('vxs',[],@isvector);
p.addParameter('maxIter',500,@isscalar);
p.addParameter('verbose',true,@islogical);

% parse
p.parse(stimulus,data,tr, varargin{:})
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


%% Set up model
% Create the model object
model = feval(p.Results.modelClass,data,stimulus,p.Results.tr,...
    'payload',p.Results.modelPayload, ...
    p.Results.modelOpts{:});

% Set model verbosity
model.verbose = verbose;

% Prep the raw data
data = model.prep(data);

% Generate seeds
seeds = model.seeds(data,vxs);


%% Fit the data
% Options for fmincon
options = optimset('Display','off', ...
    'MaxFunEvals',Inf,'MaxIter',p.Results.maxIter, ...
    'TolFun',1e-6,'TolX',1e-6);

% Pre-compute functions that will asemble the parameters in the different
% model stages
for bb = 1:model.nStages
	order = [model.floatSet{bb} model.fixSet{bb}];
	[~,sortOrder]=sort(order);
	xSort{bb} = @(x) x(sortOrder);
end

% Obtain the model bounds
[lb, ub] = model.bounds;

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

    % Squeeze the data from a cell array into a single concatenated time
    % series for the selected voxel/vertex
    datats = cell2mat(cellfun(@(x) x(vxs(ii),:),data,'UniformOutput',0))';
    
    % Apply the model cleaning step, which may include regression of
    % nuisance components.
    datats = model.clean(datats);
    
    % Pre-allocate the seed search result variables to keep par happy
    seedParams = nan(length(seeds),model.nParams);
    seedMetric = nan(length(seeds),1);
    
    % Loop over seed sets
    for ss = 1:length(seeds)
        seed = seeds{ss}(vxs(ii),:);
        x0 = seed;

        % Loop over model stages
        for bb = 1:model.nStages
            
            % Get the params in the fix and float set for this stage
            fixSet = model.fixSet{bb};
            floatSet = model.floatSet{bb};
            
            % Call the non-linear fit function
            myObj = @(x) norm(datats - model.forward(xSort{bb}([x x0(fixSet)])));
            x = fmincon(myObj,x0(floatSet),[],[],[],[], ...
                lb(floatSet),ub(floatSet), ...
                [],options);
            % Update the x0 guess with the searched params
            x0(model.floatSet{bb}) = x;
        end
        
        % Store the final params
        seedParams(ss,:) = x0;
        
        % Evaluate the model metric
        seedMetric(ss) = model.metric(datats,x0);
    end
    
    % Save the best result across seeds
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


%% Prepare the results variable
results = model.results(params, metric);

% Add the model information
results.model.class = p.Results.modelClass;
results.model.inputs = {stimulus, p.Results.tr};
results.model.opts =  p.Results.modelOpts;
results.model.payload =  p.Results.modelPayload;

% Store the calling options
results.meta.vxs = p.Results.vxs;
results.meta.tr = p.Results.tr;
results.meta.maxIter = p.Results.maxIter;

end
