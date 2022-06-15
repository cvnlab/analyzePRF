function [seeds,rvalues] = analyzePRF_computesupergridseeds(res,stimulus,data,modelfun,maxpolydeg,dimdata,dimtime,typicalgain,noisereg)

% function [seeds,rvalues] = analyzePRF_computesupergridseeds(res,stimulus,data,modelfun,maxpolydeg,dimdata,dimtime,typicalgain,noisereg)
%
% <res> is [R C] with the resolution of the stimuli
% <stimulus> is a cell vector of time x (pixels+1)
% <data> is a cell vector of X x Y x Z x time (or XYZ x time)
% <modelfun> is a function that accepts parameters (pp) and stimuli (dd) and outputs predicted time-series (time x 1)
% <maxpolydeg> is a vector of degrees (one element for each run). elements can be NaN (meaning don't use any polynomials).
% <dimdata> is number of dimensions that pertain to voxels
% <dimtime> is the dimension that is the time dimension
% <typicalgain> is a typical value for the gain in each time-series.
%   can be NaN, in which case we attempt to calculate a more proper gain seeding.
% <noisereg> is [] or a set of noise regressors (cell vector of matrices)
%
% this is an internal function called by analyzePRF.m.  this function returns <seeds>
% as a matrix of dimensions X x Y x Z x parameters (or XYZ x parameters)
% with the best seed from the super-grid.  also, returns <rvalues> as X x Y x Z
% (or XYZ x 1) with the corresponding correlation (r) values.
%
% history:
% - 2022/06/14 - the case of non-square stimulus preparations was being handled
%                incorrectly. it is now fixed.
% - 2021/11/11 - remove the forcing of single format and instead inherit from <data>.
%                also, explicitly handle some pernicious corner cases.
% - 2015/02/07 - make less memory intensive

% internal notes:
% - note that the gain seed is fake (it is not set the correct value but instead
%   to the <typicalgain>).  however, if <typicalgain> is passed as NaN, we attempt
%   to perform a more proper gain seeding.

% internal constants
eccs = [0 0.00551 0.014 0.0269 0.0459 0.0731 0.112 0.166 0.242 0.348 0.498 0.707 1];
angs = linspacecircular(0,2*pi,16);
expts = [0.5 0.25 0.125];

% calc
numvxs = prod(sizefull(data{1},dimdata));  % total number of voxels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculate a long list of potential seeds

% calc
resmx = max(res);

% calculate sigma gridding (pRF size is 1 px, 2 px, 4 px, ..., up to resmx)
maxn = floor(log2(resmx));
ssindices = 2.^(0:maxn);

% construct full list of seeds (seeds x params) [R C S G N]
fprintf('constructing seeds.\n');
allseeds = zeros(length(eccs)*length(angs)*length(ssindices)*length(expts),5);
cnt = 1;
for p=1:length(eccs)
  for q=1:length(angs)
    if p==1 && q>1  % for the center-of-gaze, only do the first angle
      continue;
    end
    for s=1:length(ssindices)
      for r=1:length(expts)
        allseeds(cnt,:) = [(1+resmx)/2 - sin(angs(q)) * (eccs(p)*resmx) ...
                           (1+resmx)/2 + cos(angs(q)) * (eccs(p)*resmx) ...
                           ssindices(s)*sqrt(expts(r)) 1 expts(r)];
        cnt = cnt + 1;
      end
    end
  end
end
allseeds(cnt:end,:) = [];  % chop because of the omission above

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% generate the predicted time-series for each seed

% generate predicted time-series [note that some time-series are all 0]
predts = zeros(sum(cellfun(@(x) size(x,1),stimulus)),size(allseeds,1),class(data{1}));  % time x seeds
temp = catcell(1,stimulus);
fprintf('generating super-grid time-series...'); tic
parfor p=1:size(allseeds,1)
  predts(:,p) = modelfun(allseeds(p,:),temp);
end
fprintf('done.'); toc
clear temp;

% % inspect for sanity on range and fineness [OBSOLETE]
% figure;
% for r=1:4:length(rrindices)
%   for c=1:4:length(ccindices)
%     for s=1:length(ssindices)
%       cnt = (r-1)*(length(ccindices)*length(ssindices)) + (c-1)*length(ssindices) + s;
%       clf;
%       plot(predts(:,cnt)');
%       title(sprintf('r=%d, c=%d, s=%d',r,c,s));
%       pause;
%     end
%   end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% prepare data and model predictions

% construct polynomials, noise regressors, and projection matrix
pregressors = {};
for p=1:length(maxpolydeg)
  if isnan(maxpolydeg(p))
    pregressors{p} = zeros(size(data{p},dimtime),0);
  else
    pregressors{p} = constructpolynomialmatrix(size(data{p},dimtime),0:maxpolydeg(p));
  end
  if ~isempty(noisereg)
    pregressors{p} = cat(2,pregressors{p},noisereg{p});
  end
end
pmatrix = projectionmatrix(blkdiag(pregressors{:}));

% project out and scale to unit length
predts = pmatrix*predts;  % time x seeds
predtslen = vectorlength(predts,1);  % 1 x seeds
predts = bsxfun(@rdivide,predts,predtslen);  % time x seeds
  % really ugly (deal with time series that have no length)
predts(repmat(predtslen==0,[size(predts,1) 1])) = 0;
predtslen(predtslen==0) = 1;
  % OLD: predts = unitlength(pmatrix*predts,1,[],0);  % time x seeds   [NOTE: some are all NaN]
  % OLD: datats = unitlength(pmatrix*squish(catcell(dimtime,data),dimdata)',1,[],0);  % time x voxels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% find the seed with the max correlation

% compute correlation and find maximum for each voxel
chunks = chunking(1:numvxs,100);
rvalues = {};
bestseedix = {};
propergain = {};
fprintf('finding best seed for each voxel.\n');
parfor p=1:length(chunks)

  % project out and scale to unit length
  datats = pmatrix*catcell(2,cellfun(@(x) subscript(squish(x,dimdata),{chunks{p} ':'}),data,'UniformOutput',0))';
  datatslen = vectorlength(datats,1);  % 1 x voxels
  datats = bsxfun(@rdivide,datats,datatslen);  % time x voxels
    % really ugly (deal with time series that have no length)
  datats(repmat(datatslen==0,[size(datats,1) 1])) = 0;
  datatslen(datatslen==0) = 1;
    % OLD: datats = unitlength(temp,1,[],0);  % time x voxels

  % voxels x 1 with index of the best seed (max corr)
  [rvalues{p},bestseedix{p}] = max(datats'*predts,[],2);  % voxels x seeds -> max corr along dim 2 [NaN is ok]
  
  % voxels x 1 with best fitting beta weight
  propergain{p} = rvalues{p} .* (vflatten(datatslen) ./ vflatten(predtslen(bestseedix{p})));
     % best fitting beta weight is something like corr(a(:),b(:)) * bb/aa

end
rvalues = catcell(1,rvalues);        % voxels x 1
bestseedix = catcell(1,bestseedix);  % voxels x 1
propergain = catcell(1,propergain);  % voxels x 1

% prepare output
rvalues = reshape(rvalues,[sizefull(data{1},dimdata) 1]);
seeds = allseeds(bestseedix,:);  % voxels x parameters

% Note: This modification is similar to the special implementation for the HCP 7T RETINOTOPY DATASET.
% Here, we set the gain seed to be 0.75 of the gain found in the grid fit.
if isnan(typicalgain)
  propergain(~isfinite(propergain)) = 0;  % SOMEWHAT HACKY BUT FOR SAFETY!!!
  seeds(:,4) = .75*max(propergain,0);  % NOTE: if the best seed in a corr sense is negative weight, just set to 0
else
  seeds(:,4) = typicalgain;
end

% prepare output
seeds = reshape(seeds,[sizefull(data{1},dimdata) size(allseeds,2)]);





%  predts(:,p) = modelfun([allseeds(p,:) flatten(hrf)],temp);
% <hrf> is T x 1 with the HRF
