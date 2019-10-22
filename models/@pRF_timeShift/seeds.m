function seeds = seeds(obj,data,vxs)
% Generate parameter seeds for the non-linear search
%
% Syntax:
%   seeds = obj.seeds(data,vxs)
%
% Description:
%   Generates a set of seed parameters for each voxel/vertex in vxs. The
%   seeds consist of a seed positioned at the stimulus center, as well as a
%   seed that is determined by grid search. A set of time-series
%   predictions are created for a plausible set of parameters that vary
%   across location within the stimulus and sigma size. The prediction that
%   is closest to the time-series data for a given voxel is found, and the
%   parameters of that prediction are assigned as the third seed.
%
% Inputs:
%   data                  - A matrix [v t] or cell array of such
%                           matricies. The fMRI time-series data across t
%                           TRs, for v vertices / voxels. The data should
%                           have bassed through the prep stage.
%   vxs                   - Vector. A list of vertices/voxels to be
%                           processed.
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   seeds                 - 3x1 cell array. Each cell contains a matrix of
%                           [v nParams] and is one of the three seed sets.
%


% Obj variables
res = obj.res;
resmx=max(res);
nParams=obj.nParams;
verbose=obj.verbose;

% Derived vars
totalVxs = size(data{1},1);

% Store the seedScale
seedScale = obj.seedScale;

% Medium scale seeds
obj.seedScale = 'medium';
x0 = initial(obj);
seeds{1} = repmat(x0,totalVxs,1);

% restore the seedScale
obj.seedScale = seedScale;


%% Internal constants
% Log-spaced eccentricities
eccs = [0 0.00551 0.014 0.0269 0.0459 0.0731 0.112 0.166 0.242 0.348 0.498 0.707 1];
% Linear spaced angles
angs = linspacecircular(0,2*pi,16);
% Plausible sigma sizes
maxn = floor(log2(resmx));
ssindices = 2.^(0:maxn);
% Plausible exponents
expts = [0.05];


%% Generate simulated time series
% Evaluate the forward model for a set of parameters that vary in x, y,
% sigma, and compressive exponent.

% construct full list of seeds (seeds x params) [R C S G N]
allseeds = zeros(length(eccs)*length(angs)*length(ssindices)*length(expts),nParams);
cnt = 1;
for p=1:length(eccs)
    for q=1:length(angs)
        if p==1 && q>1  % for the center-of-gaze, only do the first angle
            continue;
        end
        for s=1:length(ssindices)
            for r=1:length(expts)
                thisSeed = obj.initial;
                % Update the center x, y, sigma, and expt. Leave the x0
                % gain value (in position 4) unchanged.
                thisSeed([1 2 3 5]) = [...
                    (1+res(1))/2 - sin(angs(q)) * (eccs(p)*resmx), ...
                    (1+res(2))/2 + cos(angs(q)) * (eccs(p)*resmx), ...
                    ssindices(s)*sqrt(expts(r)), ...
                    expts(r) ...
                    ];
                allseeds(cnt,:) = thisSeed;
                cnt = cnt + 1;
            end
        end
    end
end


%% Time series from the forward model
% Alert the user
if verbose
    tic
    fprintf(['Generating ' num2str(size(allseeds,1)) ' seeds:\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% Create a function handle to avoid broadcasting the model obj
forward = @(p) obj.forward(p);

% Loop over all seeds and generate the predicted time series
parfor ii = 1:size(allseeds,1)

    % Update progress bar
    if verbose && mod(ii,round(size(allseeds,1)/50))==0
        fprintf('\b.\n');
    end

    % Evaluate the forward model
    predts(:,ii) = forward(allseeds(ii,:));
    
end

% report completion of loop
if verbose
    toc
    fprintf('\n');
end

% Not sure what this is for
predts = unitlength(predts,1,[],0);

% Clean the predictt
predts = obj.clean(predts);


%% Match model predictions to data
% Break the computation into chunks of 100 voxels to reduce memory
% footprint
chunks = chunking(1:length(vxs),100);
bestseedix = {};

% Set up a function to avoid broadcasting the obj to the parpool
cleanFun = @(x) obj.clean(x);

% Loop through chunks
parfor p=1:length(chunks)
    
    % time x voxels
    datats = unitlength(catcell(2,cellfun(@(x) subscript(squish(x,1),{vxs(chunks{p}) ':'}),data,'UniformOutput',0))',1,[],0);
    
    % Clean the time series
    datats = cleanFun(datats);
    
    % voxels x 1 with index of the best seed (max corr)
    [~,bestseedix{p}] = max(datats'*predts,[],2);  % voxels x seeds -> max corr along dim 2 [NaN is ok]
    
end
bestseedix = catcell(1,bestseedix);  % voxels x 1

% prepare output
vxsSeeds = allseeds(bestseedix,:);  % voxels x parameters
gridSeeds = nan(totalVxs,nParams);
gridSeeds(vxs,:)=vxsSeeds;

% Add the gridSeeds to the seed set cell array
seeds{2} = gridSeeds;

end


