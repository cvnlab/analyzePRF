function results = analyzePRF(stimulus,data,varargin)


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('stimulus',@(x)(iscell(x) || ismatrix(x)));
p.addRequired('data',@(x)(iscell(x) || ismatrix(x)));

p.addParameter('modelClass','prf_timeShift',@ischar);
p.addParameter('tr',0.8,@isscalar);
p.addParameter('vxs',[],@isvector);
p.addParameter('hrf',[],@isvector);
p.addParameter('maxIter',500,@isscalar);
p.addParameter('typicalGain',30,@isscalar);
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
dimtime = 2;
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



%% Data prep
nAcqs = length(data);

% deal with data badness (set bad voxels to be always all 0)
bad = cellfun(@(x) any(~isfinite(x),dimtime) | all(x==0,dimtime),data,'UniformOutput',0);  % if non-finite or all 0
bad = any(cat(dimtime,bad{:}),dimtime);  % badness in ANY run
for ii=1:nAcqs
    data{ii}(repmat(bad,[ones(1,dimdata) size(data{ii},dimtime)])) = 0;
end


%% HRF prep
if isempty(p.Results.hrf)
    hrf = getcanonicalhrf(p.Results.tr,p.Results.tr)';
else
    hrf = p.Results.hrf;
end


%% Set up model

% Create the model object
model = pRF_timeShift(stimulus,res,hrf);

% Assign the typical gain parameter
model.gain = p.Results.typicalGain;

% Set model verbosity
model.verbose = verbose;

% Generate seeds
seeds = model.seeds(data,vxs);

% % init
% seeds = [];
%
% % generic large seed
% modelObj.seedScale = 'large';
% seeds = [ seeds; modelObj.initial() ];
%
% % generic small seed
% modelObj.seedScale = 'small';
% seeds = [ seeds; modelObj.initial() ];
%
%
% % make a function that individualizes the seeds
% if exist('supergridseeds','var')
%     seedfun = @(vx) [[seeds];
%         [subscript(squish(supergridseeds,dimdata),{vx ':'})]];
% else
%     seedfun = @(vx) [seeds];
% end


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
    datastd = std(data4);
    data5 = data4 / datastd;
    data6 = data5 - mean(data5);
    
    %     seeds = seedfun(ii);
    %     for ss = 1:size(seeds,1)
    %
    
    seed = seeds(vxs(ii),:);
    x0 = seed;
    
    for bb = 1:model.nStages
        x0 = lsqcurvefit(...
            @(x,y) double(model.forward([x x0(model.fixSet{bb})]) / datastd),...
            x0(model.floatSet{bb}),...
            [],...
            double(data6),[],[],options);
        x0 = [x0 seed(model.fixSet{bb})];
    end
    
    parParams(ii,:) = x0;
    parMetric(ii) = model.metric(data6,x0);
    
    %     end
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
results = model.results(obj, params, metric);
results.meta.p = p;

end
