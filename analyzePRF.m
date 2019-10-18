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


% init
seeds = [];

% generic large seed
modelObj.seedScale = 'large';
seeds = [ seeds; modelObj.initial() ];

% generic small seed
modelObj.seedScale = 'small';
seeds = [ seeds; modelObj.initial() ];

% super-grid seed
maxPolyDeg = cellfun(@(x) round(size(x,dimtime)*p.Results.tr/60/2),data);
[supergridseeds,~] = analyzePRFcomputesupergridseeds(modelObj,res,stimulus,data, ...
    maxPolyDeg,dimdata,dimtime, ...
    p.Results.typicalGain,[]);

% make a function that individualizes the seeds
if exist('supergridseeds','var')
    seedfun = @(vx) [[seeds];
        [subscript(squish(supergridseeds,dimdata),{vx ':'})]];
else
    seedfun = @(vx) [seeds];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE STIMULUS AND DATA





options = optimset('Display','off','FunValCheck','on', ...
    'MaxFunEvals',Inf,'MaxIter',p.Results.maxiter, ...
    'TolFun',1e-6,'TolX',1e-6, ...
    'OutputFcn',@(a,b,c) outputfcnsanitycheck(a,b,c,1e-6,10));
options.Algorithm = 'levenberg-marquardt';

vnum = length(vxs);
for ii=1:vnum
    data2 = cellfun(@(x) x(ii,:),data,'UniformOutput',0);
    data3 = @(vxs) cellfun(@(x) subscript(squish(x,dimdata),{vxs ':'})',data2,'UniformOutput',0);
    data4 = catcell(1,data3(1));
    datastd = std(data4);
    data5 = data4 / datastd;
    data6 = data5 - mean(data5);
    
    seeds = seedfun(ii);
    for ss = 1:size(seeds,1)

        
        seed = seeds(ss,:);
        x0 = seed;

        for bb = 1:modelObj.nStages
            x0 = lsqcurvefit(...
                @(x,y) double(modelObj.forward([x x0(modelObj.fixSet{bb})]) / datastd),...
                x0(modelObj.floatSet{bb}),...
                [],...
                double(data6),[],[],options);
            x0 = [x0 seed(modelObj.fixSet{bb})];
        end
        
        params(ss,:) = x0;
        metric(ss) = calccorrelation(data6,modelObj.forward(x));
        
    end
    
end


