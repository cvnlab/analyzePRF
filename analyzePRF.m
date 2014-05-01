function [results,supergridseeds] = ...
  analyzePRF(stimulus,data,tr,vxs,wantglmdenoise,seedmode,wantfithrf,xvalmode,numperjob,maxiter,display,typicalgain)

% function [results,supergridseeds] = ...
%   analyzePRF(stimulus,data,tr,vxs,wantglmdenoise,seedmode,wantfithrf,xvalmode,numperjob,maxiter,display,typicalgain)
%
% <stimulus> provides the apertures as a cell vector of R x C x time.
%   values should be in [0,1].  the number of time points can differ across runs.
% <data> provides the data as a cell vector of voxels x time.  can also be
%   X x Y x Z x time.  the number of time points should match the number of 
%   time points in <stimulus>.
% <tr> is the TR in seconds (e.g. 1.5)
% <vxs> (optional) is a vector of voxel indices to analyze.  (we automatically
%   sort the voxel indices and ensure uniqueness.)  default to 1:V where 
%   V is total number of voxels.  to reduce computational time, you may want 
%   to create a binary brain mask and perform find() on it.
% <wantglmdenoise> (optional) is whether to use GLMdenoise to determine
%   nuisance regressors to add into the PRF model.  note that in order to use
%   this feature, there must be at least two runs (and conditions must repeat
%   across runs).  we automatically determine the GLM design matrix based on
%   the contents of <stimulus>.  special case is to pass in the noise regressors 
%   directly (e.g. from a previous call).  default: 0.
% <seedmode> (optional) is
%   0 means use default initial seeds
%   1 means use best seed based on super-grid
%   2 means use random seed
%   3 means 0 and 1
%   M where M is X x Y x Z x parameters with the seed to use for each voxel.
%   default: 0.  sketch out this strategy. maybe make a figure.  helps you think about bias.
% <wantfithrf> (optional) is whether to fit the HRF on a voxel-by-voxel basis.
%   (INTERACTION WITH GLMDENOISE>..??.)
%   default: 0.
% <xvalmode> (optional) is
%   0 means just fit all the data
%   1 means two-fold cross-validation (first half of runs; second half of runs)
%   2 means two-fold cross-validation (first half of each run; second half of each run)
%   default: 0.  (note that we round when halving.)
% <numperjob> (optional) is
%   [] means to run locally (not on the cluster)
%   N where N is a positive integer indicating the number of voxels to
%     analyze in each cluster job.
%   default: [].
% <maxiter> (optional) is maximum number of iterations.
%   default: 500.
% <display> (optional)
% <typicalgain> (optional)
%
% analyze PRF data and return the results.
% <supergridseeds> is X x Y x Z x parameters with the super-grid seed selected.
%
% notes:
% - pRF gain is restricted to be positive
%
% history:
% 2014/04/27 - gain seed is now set to 0; add gain to the output
% 2014/04/29 - use typicalgain now (default 10). allow display input.
%
% example:

% internal notes:
% - convert inputs to options struct!
% - for cluster mode, need to make sure fitnonlinearmodel is compiled (compilemcc.m)
% - to check whether local minima are a problem, can look at results.resnorms
% - having clear input and clear output will be good for integration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REPORT

fprintf('*** analyzePRF: started at %s. ***\n',datestr(now));
stime = clock;  % start time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTERNAL CONSTANTS

% define
remotedir = '/scratch/knk/input/';
remotedir2 = '/scratch/knk/output/';
remotelogin = 'knk@login1.chpc.wustl.edu';
remoteuser = 'knk';
corrthresh = .99;  % used in determining which apertures are the same

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP AND PREPARATION

% massage
if ~iscell(stimulus)
  stimulus = {stimulus};
end
if ~iscell(data)
  data = {data};
end

% calc
is3d = size(data{1},4) > 1;
if is3d
  dimdata = 3;
  dimtime = 4;
  xyzsize = sizefull(data{1},3);
else
  dimdata = 1;
  dimtime = 2;
  xyzsize = size(data{1},1);
end
numvxs = prod(xyzsize);

% calc
res = sizefull(stimulus{1},2);
resmx = max(res);
numruns = length(data);

% input
if ~exist('vxs','var') || isempty(vxs)
  vxs = 1:numvxs;
end
if ~exist('wantglmdenoise','var') || isempty(wantglmdenoise)
  wantglmdenoise = 0;
end
if ~exist('seedmode','var') || isempty(seedmode)
  seedmode = 0;
end
if ~exist('wantfithrf','var') || isempty(wantfithrf)
  wantfithrf = 0;
end
if ~exist('xvalmode','var') || isempty(xvalmode)
  xvalmode = 0;
end
if ~exist('numperjob','var') || isempty(numperjob)
  numperjob = [];
end
if ~exist('maxiter','var') || isempty(maxiter)
  maxiter = 500;
end
if ~exist('display','var') || isempty(display)
  display = 'final';
end
if ~exist('typicalgain','var') || isempty(typicalgain)
  typicalgain = 10;
end

% calc
usecluster = ~isempty(numperjob);

% prepare stimuli
for p=1:length(stimulus)
  stimulus{p} = squish(stimulus{p},2)';  % frames x pixels
  stimulus{p} = [stimulus{p} p*ones(size(stimulus{p},1),1)];  % add a dummy column to indicate run breaks
  stimulus{p} = single(stimulus{p});  % make single to save memory
end

% deal with data badness (set bad voxels to always all 0)
bad = cellfun(@(x) any(~isfinite(x),dimtime) | all(x==0,dimtime),data,'UniformOutput',0);  % if non-finite or all 0
bad = any(cat(dimtime,bad{:}),dimtime);  % badness in ANY run
for q=1:numruns
  data{q}(repmat(bad,[ones(1,dimdata) size(data{q},dimtime)])) = 0;
end

% calc mean volume
meanvol = mean(catcell(dimtime,data),dimtime);

% what HRF should we start from?
hrf = getcanonicalhrf(tr,tr)';
numinhrf = length(hrf);

% what polynomials should we use?
maxpolydeg = cellfun(@(x) round(size(x,dimtime)*tr/60/2),data);

% initialize cluster stuff
if usecluster
  localfilestodelete = {};
  remotefilestodelete = {};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% APPLY GLMDENOISE [perhaps make this a subfunction?]

if isequal(wantglmdenoise,1)

  %%%%% figure out GLM design matrix
  
  % concatenate everything and drop the dummy column (X is pixels x frames)
  X = catcell(1,stimulus)';
  X(end,:) = [];
  
  % figure out where blanks are (logical vector, 1 x frames)
  blanks = all(X==0,1);

  % normalize each frame (in preparation for computing correlations)
  X = unitlength(X,1,[],0);
  X(:,blanks) = 0;

  % initialize the result (this will grow in size on-the-fly)
  glmdesign = zeros(size(X,2),0);

  % do the loop
  wh = find(~blanks);
  cnt = 1;
  while ~isempty(wh)
    ix = wh(1);                        % pick the first one to process
    corrs = X(:,ix)' * X;              % compute correlation with all frames
    spots = find(corrs > corrthresh);  % any frame with r > .99 counts as the same
    glmdesign(spots,cnt) = 1;          % add to design matrix
    X(:,spots) = 0;                    % blank out (since we're done with those columns)
    blanks(spots) = 1;                 % update the list of blanks
    wh = find(~blanks);
    cnt = cnt + 1;
  end

  % finally, un-concatenate the results
  glmdesign = splitmatrix(glmdesign,1,cellfun(@(x) size(x,1),stimulus));
  
  % clean up
  clear X;
  
  %%%%% run GLMdenoise to get the noise regressors

  % what directory to save results to?
  glmdenoisedir = [tempname];
  assert(mkdir(glmdenoisedir));

  % call GLMdenoise
  fprintf('using GLMdenoise figure directory %s\n',[glmdenoisedir '/GLMdenoisefigures']);
  if wantfithrf
    hrfmodel = [];
    hrfknobs = [];
  else
    hrfmodel = 'assume';
    hrfknobs = hrf;
  end
  results = GLMdenoisedata(glmdesign,data,tr,tr,hrfmodel,hrfknobs, ...
                           struct('numboots',0), ...
                           [glmdenoisedir '/GLMdenoisefigures']);

  % get the noise regressors
  noisereg = cellfun(@(x) x(:,1:results.pcnum),results.pcregressors,'UniformOutput',0);
  
  % get new HRF
  if wantfithrf
    hrf = results.modelmd{1};
    numinhrf = length(hrf);
  end

  % save 'results' to a file
  file0 = [glmdenoisedir '/GLMdenoise.mat'];
  fprintf('saving GLMdenoise results to %s (in case you want them).\n',file0);
  results = rmfield(results,{'pcweights' 'models' 'modelse'});  % remove some boring fields
  save(file0,'results','noisereg');
  
  % clean up
  clear results;

elseif isequal(wantglmdenoise,0)

  noisereg = [];

else

  noisereg = wantglmdenoise;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE MODEL

% pre-compute some cache
[d,xx,yy] = makegaussian2d(resmx,2,2,2,2);

% define the model (parameters are R C S G N [HRF])
modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),pp(5+(1:numinhrf))',dd(:,prod(res)+1));
model = {{[] [1-res(1)+1 1-res(2)+1 0    0   NaN repmat(NaN,[1 numinhrf]);
              2*res(1)-1 2*res(2)-1 Inf  Inf Inf repmat(Inf,[1 numinhrf])] modelfun} ...
         {@(ss)ss [1-res(1)+1 1-res(2)+1 0    0   0   repmat(NaN,[1 numinhrf]);
                   2*res(1)-1 2*res(2)-1 Inf  Inf Inf repmat(Inf,[1 numinhrf])] @(ss)modelfun} ...
         {@(ss)ss [1-res(1)+1 1-res(2)+1 0    0   0   repmat(-Inf,[1 numinhrf]);
                   2*res(1)-1 2*res(2)-1 Inf  Inf Inf repmat(Inf,[1 numinhrf])] @(ss)modelfun}};

% if don't want to fit the HRF, exclude the last model step
if ~wantfithrf
  model = model(1:2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE SEEDS

% default output for supergridseeds is []
supergridseeds = [];

% now define the seeds
if isequal(seedmode,0)

  % compute typical seeds (one large and one small at center)
  seeds = [(1+res(1))/2 (1+res(2))/2 resmx/4*sqrt(0.5)  typicalgain 0.5;
           (1+res(1))/2 (1+res(2))/2 resmx/16*sqrt(0.5) typicalgain 0.5;
           (1+res(1))/2 (1+res(2))/2 resmx/64*sqrt(0.5) typicalgain 0.5];
           
  % tack on HRF
  seeds(:,end+(1:numinhrf)) = repmat(flatten(hrf),[size(seeds,1) 1]);

elseif isequal(seedmode,1)

  % compute super-grid seeds
  supergridseeds = computesupergridseeds(res,stimulus,data,modelfun,hrf,maxpolydeg,dimdata,dimtime);

  % make a function that individualizes the seeds and tacks on HRF
  seeds = @(vx) [subscript(squish(supergridseeds,dimdata),{vx ':'}) flatten(hrf)];

elseif isequal(seedmode,2)

  % compute random seed (angle is random between 0 and 2*pi; ecc is random between 0 and resmx/2)
  tempfun = @(ang0) [-sin(ang0) cos(ang0)];
  seeds = @(vx) [[(1+res(1))/2 (1+res(2))/2] + (rand*resmx/2)*tempfun(rand*2*pi) resmx/4*sqrt(0.5) 10 0.5 flatten(hrf)];

elseif isequal(seedmode,3)

  % compute super-grid seeds
  supergridseeds = computesupergridseeds(res,stimulus,data,modelfun,hrf,maxpolydeg,dimdata,dimtime,typicalgain);

  % make a function that individualizes the seeds and tacks on HRF
  seeds = @(vx) cat(2,[subscript(squish(supergridseeds,dimdata),{vx ':'});
                 (1+res(1))/2 (1+res(2))/2 resmx/4*sqrt(0.5)  typicalgain 0.5;
                 (1+res(1))/2 (1+res(2))/2 resmx/16*sqrt(0.5) typicalgain 0.5;
                 (1+res(1))/2 (1+res(2))/2 resmx/64*sqrt(0.5) typicalgain 0.5],repmat(flatten(hrf),[4 1]));

else

  % make a function that individualizes the seeds and tacks on HRF
  seeds = @(vx) [subscript(squish(seedmode,      dimdata),{vx ':'}) flatten(hrf)];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE RESAMPLING STUFF

% define wantresampleruns and resampling
switch xvalmode
case 0
  wantresampleruns = [];
  resampling = 0;
case 1
  wantresampleruns = 1;
  half1 = copymatrix(zeros(1,length(data)),1:round(length(data)/2),1);
  half2 = ~half1;
  resampling = [(1)*half1 + (-1)*half2;
                (-1)*half1 + (1)*half2];
case 2
  wantresampleruns = 0;
  resampling = [];
  for p=1:length(data)
    half1 = copymatrix(zeros(1,size(data{p},2)),1:round(size(data{p},2)/2),1);
    half2 = ~half1;
    resampling = cat(2,resampling,[(1)*half1 + (-1)*half2;
                                   (-1)*half1 + (1)*half2]);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE STIMULUS AND DATA

%%%%% CLUSTER CASE

if usecluster

  % save stimulus and transport to the remote server
  while 1
    filename0 = sprintf('stim%s.mat',randomword(5));  % file name
    localfile0 = [tempdir '/' filename0];             % local path to file
    remotefile0 = [remotedir '/' filename0];          % remote path to file
  
    % redo if file already exists locally or remotely
    if exist(localfile0) || 0==unix(sprintf('ssh %s ls %s',remotelogin,remotefile0))
      continue;
    end
  
    % save file and transport it
    save(localfile0,'stimulus');
    assert(0==unix(sprintf('rsync -av %s %s:"%s/"',localfile0,remotelogin,remotedir)));
  
    % record
    localfilestodelete{end+1} = localfile0;
    remotefilestodelete{end+1} = remotefile0;
  
    % stop
    break;
  end
  clear stimulus;  % don't let it bleed through anywhere!

  % define stimulusINPUT
  stimulusINPUT = @() loadmulti(remotefile0,'stimulus');

  % save data and transport to the remote server
  while 1
    filename0 = sprintf('data%s',randomword(5));   % directory name that will contain 001.bin, etc.
    localfile0 = [tempdir '/' filename0];          % local path to dir
    remotefile0 = [remotedir '/' filename0];       % remote path to dir
  
    % redo if dir already exists locally or remotely
    if exist(localfile0) || 0==unix(sprintf('ssh %s ls %s',remotelogin,remotefile0))
      continue;
    end
  
    % save files and transport them
    assert(mkdir(localfile0));
    for p=1:numruns
      savebinary([localfile0 sprintf('/%03d.bin',p)],'single',squish(data{p},dimdata)');  % notice squish
    end
    assert(0==unix(sprintf('rsync -av %s %s:"%s/"',localfile0,remotelogin,remotedir)));

    % record
    localfilestodelete{end+1} = localfile0;
    remotefilestodelete{end+1} = remotefile0;

    % stop
    break;
  end
  clear data;

  % define dataINPUT
  binfiles = cellfun(@(x) [remotefile0 sprintf('/%03d.bin',x)],num2cell(1:numruns),'UniformOutput',0);
  dataINPUT = @(vxs) cellfun(@(x) double(loadbinary(x,'single',[0 numvxs],-vxs)),binfiles,'UniformOutput',0);

  % prepare the output directory name
  while 1
    filename0 = sprintf('prfresults%s',randomword(5));
    localfile0 = [tempdir '/' filename0];
    remotefile0 = [remotedir2 '/' filename0];
    if exist(localfile0) || 0==unix(sprintf('ssh %s ls %s',remotelogin,remotefile0))
      continue;
    end
    localfilestodelete{end+1} = localfile0;
    localfilestodelete{end+1} = [localfile0 '.mat'];  % after consolidation
    remotefilestodelete{end+1} = remotefile0;
    break;
  end
  outputdirlocal = localfile0;
  outputdirremote = remotefile0;
  outputdirINPUT = outputdirremote;

%%%%% NON-CLUSTER CASE

else

  stimulusINPUT = {stimulus};
  dataINPUT = @(vxs) cellfun(@(x) subscript(squish(x,dimdata),{vxs ':'})',data,'UniformOutput',0);
  outputdirINPUT = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE OPTIONS

% construct the options struct
if iscell(noisereg)
  noisereg0 = {noisereg};
else
  noisereg0 = noisereg;
end
opt = struct( ...
  'outputdir',outputdirINPUT, ...
  'stimulus',stimulusINPUT, ...
  'data',dataINPUT, ...
  'vxs',vxs, ...
  'model',{model}, ...
  'seed',seeds, ...
  'optimoptions',{{'Display' display 'Algorithm' 'levenberg-marquardt' 'MaxIter' maxiter}}, ...  % 'iter'
  'wantresampleruns',wantresampleruns, ...
  'resampling',resampling, ...
  'metric',@calccod, ...
  'maxpolydeg',maxpolydeg, ...
  'wantremovepoly',1, ...
  'extraregressors',noisereg0, ...
  'wantremoveextra',0, ...
  'dontsave',{{'modelfit' 'opt' 'vxsfull' 'modelpred' 'testdata'}});  % 'resnorms' 'numiters' 

%  'outputfcn',@(a,b,c,d) pause2(.1) | outputfcnsanitycheck(a,b,c,1e-6,10) | outputfcnplot(a,b,c,1,d), ...
        %'outputfcn',@(a,b,c,d) pause2(.1) | outputfcnsanitycheck(a,b,c,1e-6,10) | outputfcnplot(a,b,c,1,d));
        %   % debugging:
        %   chpcstimfile = '/stone/ext1/knk/HCPretinotopy/conimagesB.mat';
        %   chpcdatadir2 = outputdir2;  % go back
        %   opt.outputdir='~/temp1';
        %   profile on;
        %   results = fitnonlinearmodel(opt,100,100);
        %   results = fitnonlinearmodel(opt,1,715233);
        %   profsave(profile('info'),'~/inout/profile_results');
        % %   modelfit = feval(modelfun,results.params,feval(stimulusINPUT));
        % %   thedata = feval(dataINPUT,52948);
        % %   pmatrix = projectionmatrix(constructpolynomialmatrix(304,0:3));
        % %   figure; hold on;
        % %   plot(pmatrix*thedata,'k-');
        % %   plot(pmatrix*modelfit,'r-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIT MODEL

%%%%% CLUSTER CASE

if usecluster

  % submit jobs
  jobnames = {};
  jobnames = [jobnames {makedirid(opt.outputdir,1)}];
  jobids = [];
  jobids = [jobids chpcrun(jobnames{end},'fitnonlinearmodel',numperjob,1,ceil(length(vxs)/numperjob),[],{'data' 'stimulus' 'bad' 'd' 'xx' 'yy' 'modelfun' 'model' 'stimulusINPUT' 'dataINPUT'})];

  % record additional files to delete
  for p=1:length(jobnames)
    remotefilestodelete{end+1} = sprintf('~/sgeoutput/job_%s.*',jobnames{p});  % .o and .e files
    remotefilestodelete{end+1} = sprintf('~/mcc/job_%s.mat',jobnames{p});
    localfilestodelete{end+1} = sprintf('~/mcc/job_%s.mat',jobnames{p});
  end

  % wait for jobs to finish
  sgewaitjobs(jobnames,jobids,remotelogin,remoteuser);

  % download the results
  assert(0==unix(sprintf('rsync -a %s:"%s" "%s/"',remotelogin,outputdirremote,tempdir)));

  % consolidate the results
  fitnonlinearmodel_consolidate(outputdirlocal);

  % load the results
  a1 = load([outputdirlocal '.mat']);

%%%%% NON-CLUSTER CASE

else

  a1 = fitnonlinearmodel(opt);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE OUTPUT

% calc
numfits = size(a1.params,1);

% init
clear results;
results.ang =      NaN*zeros(numvxs,numfits);
results.ecc =      NaN*zeros(numvxs,numfits);
results.expt =     NaN*zeros(numvxs,numfits);
results.rfsize =   NaN*zeros(numvxs,numfits);
results.R2 =       NaN*zeros(numvxs,numfits);
results.gain =     NaN*zeros(numvxs,numfits);
results.resnorms = cell(numvxs,1);
results.numiters = cell(numvxs,1);
results.hrf =      NaN*zeros(numvxs,numinhrf,numfits);

% massage model parameters for output and put in 'results' struct
results.ang(vxs,:) =    permute(mod(atan2((1+res(1))/2 - a1.params(:,1,:), ...
                                          a1.params(:,2,:) - (1+res(2))/2),2*pi)/pi*180,[3 1 2]);
results.ecc(vxs,:) =    permute(sqrt(((1+res(1))/2 - a1.params(:,1,:)).^2 + (a1.params(:,2,:) - (1+res(2))/2).^2),[3 1 2]);
results.expt(vxs,:) =   permute(posrect(a1.params(:,5,:)),[3 1 2]);
results.rfsize(vxs,:) = permute(abs(a1.params(:,3,:)) ./ sqrt(posrect(a1.params(:,5,:))),[3 1 2]);
results.R2(vxs,:) =     permute(a1.trainperformance,[2 1]);
results.gain(vxs,:) =   permute(posrect(a1.params(:,4,:)),[3 1 2]);
results.resnorms(vxs) = a1.resnorms;
results.numiters(vxs) = a1.numiters;
results.hrf(vxs,:,:) =  permute(a1.params(:,5+(1:numinhrf),:),[3 2 1]);

% reshape
results.ang = reshape(results.ang,[xyzsize numfits]);
results.ecc = reshape(results.ecc,[xyzsize numfits]);
results.expt = reshape(results.expt,[xyzsize numfits]);
results.rfsize = reshape(results.rfsize,[xyzsize numfits]);
results.R2 = reshape(results.R2,[xyzsize numfits]);
results.gain = reshape(results.gain,[xyzsize numfits]);
results.resnorms = reshape(results.resnorms,[xyzsize 1]);
results.numiters = reshape(results.numiters,[xyzsize 1]);
results.hrf = reshape(results.hrf,[xyzsize numinhrf numfits]);

% add some more stuff
results.meanvol = meanvol;
results.noisereg = noisereg;
results.params = a1.params;

% save 'results' to a temporary file so we don't lose these precious results!
file0 = [tempname '.mat'];
fprintf('saving results to %s (just in case).\n',file0);
save(file0,'results');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLEAN UP

%%%%% CLUSTER CASE

if usecluster

  % delete local files and directories
  for p=1:length(localfilestodelete)
    if exist(localfilestodelete{p},'dir')  % first dir, then file
      rmdir(localfilestodelete{p},'s');
    elseif exist(localfilestodelete{p},'file')
      delete(localfilestodelete{p});
    end
  end

  % delete remote files and directories
  for p=1:length(remotefilestodelete)
    assert(0==unix(sprintf('ssh %s "rm -rf %s"',remotelogin,remotefilestodelete{p})));
  end

%%%%% NON-CLUSTER CASE

else

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REPORT

fprintf('*** analyzePRF: ended at %s (%.1f minutes). ***\n', ...
        datestr(now),etime(clock,stime)/60);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function seeds = computesupergridseeds(res,stimulus,data,modelfun,hrf,maxpolydeg,dimdata,dimtime,typicalgain)

% function seeds = computesupergridseeds(res,stimulus,data,modelfun,hrf,maxpolydeg,dimdata,dimtime,typicalgain)
%
% <res> is [R C] with the resolution of the stimuli
% <stimulus> is a cell vector of time x (pixels+1)
% <data> is a cell vector of X x Y x Z x time
% <modelfun> is a function that accepts parameters (pp) and stimuli (dd) and outputs predicted time-series (time x 1)
% <hrf> is T x 1 with the HRF
% <maxpolydeg> is a vector of degrees
% <dimdata> is number of dimensions that pertain to voxels
% <dimtime> is the dimension that is the time dimension
% <typicalgain>
%
% return a matrix of dimensions X x Y x Z x parameters with the best seed from the super-grid.
% the seed does NOT include the HRF parameters.
%
% internal notes:
% - note that the gain seed is fake (it is not set correctly)

% define
eccs = [0 0.00551 0.014 0.0269 0.0459 0.0731 0.112 0.166 0.242 0.348 0.498 0.707 1];
angs = linspacecircular(0,2*pi,16);
expts = [0.5 0.25 0.125];

% calc
resmx = max(res);

% calculate sigma gridding (pRF size is 1 px, 2 px, 4 px, ..., up to resmx)
maxn = floor(log2(resmx));
ssindices = 2.^(0:maxn);

% construct full list of seeds (seeds x params)
allseeds = zeros(length(eccs)*length(angs)*length(ssindices)*length(expts),5);
cnt = 1;
for p=1:length(eccs)
  for q=1:length(angs)
    if p==1 && q>1  % for the center-of-gaze, only do the first angle
      continue;
    end
    for s=1:length(ssindices)
      for r=1:length(expts)
        allseeds(cnt,:) = [(1+res(1))/2 - sin(angs(q)) * (eccs(p)*resmx) ...
                           (1+res(2))/2 + cos(angs(q)) * (eccs(p)*resmx) ...
                           ssindices(s)*sqrt(expts(r)) 1 expts(r)];
        cnt = cnt + 1;
      end
    end
  end
end
allseeds(cnt:end,:) = [];  % chop because of the omission above

% generate predicted time-series [note that some time-series are all 0]
predts = zeros(sum(cellfun(@(x) size(x,1),stimulus)),size(allseeds,1),'single');  % time x seeds
temp = catcell(1,stimulus);
fprintf('generating super-grid time-series...'); tic
parfor p=1:size(allseeds,1)
  predts(:,p) = modelfun([allseeds(p,:) flatten(hrf)],temp);
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

% construct polynomials and projection matrix
polyregressors = {};
for p=1:length(maxpolydeg)
  polyregressors{p} = constructpolynomialmatrix(size(data{p},dimtime),0:maxpolydeg(p));
end
pmatrix = projectionmatrix(blkdiag(polyregressors{:}));

% project out polynomials and scale to unit length
predts = unitlength(pmatrix*predts,                                1,[],0);  % time x seeds   [NOTE: some are all NaN]
datats = unitlength(pmatrix*squish(catcell(dimtime,data),dimdata)',1,[],0);  % time x voxels

% compute correlation and find maximum for each voxel
chunks = chunking(1:size(datats,2),100);
bestseedix = {};
parfor p=1:length(chunks)
  % voxels x 1 with index of the best seed (max corr)
  [mx,bestseedix{p}] = max(datats(:,chunks{p})' * predts,[],2);  % voxels x seeds -> max corr along dim 2 [NaN is ok]
end
bestseedix = catcell(1,bestseedix);  % voxels x 1

% massage output
temp = allseeds(bestseedix,:);
temp(:,4) = typicalgain;  % set gain to typical gain
seeds = reshape(temp,[sizefull(data{1},dimdata) size(allseeds,2)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUNK

    % OLD: before wantfithrf
    % % define the model (parameters are R C S G N)
    % modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));
    % model = {{[] [1-res(1)+1 1-res(2)+1 0   0   NaN;
    %               2*res(1)-1 2*res(2)-1 Inf Inf Inf] modelfun} ...
    %          {@(ss)ss [1-res(1)+1 1-res(2)+1 0   0   0;
    %                    2*res(1)-1 2*res(2)-1 Inf Inf Inf] @(ss)modelfun}};
