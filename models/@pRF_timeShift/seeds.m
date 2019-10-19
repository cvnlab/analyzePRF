function seeds = seeds(obj,data,vxs)


%% Obj variables
stimulus = obj.stimulus;
res = obj.res;
resmx=max(res);
nParams=obj.nParams;
verbose=obj.verbose;

%% Internal constants
% Log-spaced eccentricities
eccs = [0 0.00551 0.014 0.0269 0.0459 0.0731 0.112 0.166 0.242 0.348 0.498 0.707 1];
% Linear spaced angles
angs = linspacecircular(0,2*pi,16);
% Plausible sigma sizes
maxn = floor(log2(resmx));
ssindices = 2.^(0:maxn);
% Plausible exponents
expts = [0.5 0.25 0.125];


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
% generate predicted time-series [note that some time-series are all 0]
predts = zeros(sum(cellfun(@(x) size(x,1),stimulus)),size(allseeds,1),'single');  % time x seeds

forward = @(p) obj.forward(p);


% Alert the user
if verbose
    tic
    fprintf(['Generating ' num2str(size(allseeds,1)) ' seeds:\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

parfor p=1:size(allseeds,1)

    % Update progress bar
    if verbose && mod(p,round(size(allseeds,1)/50))==0
        fprintf('\b.\n');
    end

    % Evaluate the forward model
    predts(:,p) = forward(allseeds(p,:));
    
end

% report completion of loop
if verbose
    toc
    fprintf('\n');
end

predts = unitlength(predts,1,[],0);

% Clean the time series
predts = obj.clean(predts);


%% Match model predictions to data

% Break the computation into chunks of voxels to reduce memory footprint
chunks = chunking(1:length(vxs),100);
bestseedix = {};

% Set up a function to avoid broadcasting the obj to the parpool
cleanFun = @(x) obj.clean(x);

% Loop through chunks
parfor p=1:length(chunks)
    
    % time x voxels
    datats = unitlength(catcell(2,cellfun(@(x) subscript(squish(x,1),{vxs(chunks{p}) ':'}),data,'UniformOutput',0))',1,[],0);
    
    % Implement the time series clean stage
    datats = cleanFun(datats);
    
    % voxels x 1 with index of the best seed (max corr)
    [~,bestseedix{p}] = max(datats'*predts,[],2);  % voxels x seeds -> max corr along dim 2 [NaN is ok]
    
end
bestseedix = catcell(1,bestseedix);  % voxels x 1

% prepare output
vxsSeeds = allseeds(bestseedix,:);  % voxels x parameters
totalVxs = size(data{1},1);
seeds = nan(totalVxs,nParams);
seeds(vxs,:)=vxsSeeds;

seeds = {seeds};

end


