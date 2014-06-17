function results = analyzePRF(stimulus,data,tr,options)

% function results = analyzePRF(stimulus,data,tr,options)
%
% <stimulus> provides the apertures as a cell vector of R x C x time.
%   values should be in [0,1].  the number of time points can differ across runs.
% <data> provides the data as a cell vector of voxels x time.  can also be
%   X x Y x Z x time.  the number of time points should match the number of 
%   time points in <stimulus>.
% <tr> is the TR in seconds (e.g. 1.5)
% <options> (optional) is a struct with the following fields:
%   <vxs> (optional) is a vector of voxel indices to analyze.  (we automatically
%     sort the voxel indices and ensure uniqueness.)  default is 1:V where 
%     V is total number of voxels.  to reduce computational time, you may want 
%     to create a binary brain mask, perform find() on it, and use the result as <vxs>.
%   <wantglmdenoise> (optional) is whether to use GLMdenoise to determine
%     nuisance regressors to add into the PRF model.  note that in order to use
%     this feature, there must be at least two runs (and conditions must repeat
%     across runs).  we automatically determine the GLM design matrix based on
%     the contents of <stimulus>.  special case is to pass in the noise regressors 
%     directly (e.g. from a previous call).  default: 0.
%   <hrf> (optional) is a column vector with the hemodynamic response function (HRF)
%     to use in the model.  the first value of <hrf> should be coincident with the onset
%     of the stimulus, and the HRF should indicate the timecourse of the response to
%     a stimulus that lasts for one TR.  default is to use a canonical HRF (calculated
%     using getcanonicalhrf(tr,tr)').
%   <maxpolydeg> (optional) is a non-negative integer indicating the maximum polynomial
%     degree to use for drift terms.  can be a vector whose length matches the number
%     of runs in <data>.  default is to use round(L/2) where L is the number of minutes
%     in the duration of a given run.
%   <seedmode> (optional) is a vector consisting of one or more of the
%     following values (we automatically sort and ensure uniqueness):
%       0 means use generic large PRF seed
%       1 means use generic small PRF seed
%       2 means use best seed based on super-grid
%     default: [0 1 2].
%   <xvalmode> (optional) is
%     0 means just fit all the data
%     1 means two-fold cross-validation (first half of runs; second half of runs)
%     2 means two-fold cross-validation (first half of each run; second half of each run)
%     default: 0.  (note that we round when halving.)
%   <numperjob> (optional) is
%     [] means to run locally (not on the cluster)
%     N where N is a positive integer indicating the number of voxels to analyze in each 
%       cluster job.  this option requires a customized computational setup!
%     default: [].
%   <maxiter> (optional) is the maximum number of iterations.  default: 500.
%   <display> (optional) is 'iter' | 'final' | 'off'.  default: 'iter'.
%   <typicalgain> (optional) is a typical value for the gain in each time-series.
%     default: 10.
%
% Analyze pRF data and return the results.
%
% The results structure contains the following fields:
% <ang> contains pRF angle estimates.  Values range between 0 and 360 degrees.
%   0 corresponds to the right horizontal meridian, 90 corresponds to the upper vertical
%   meridian, and so on.
% <ecc> contains pRF eccentricity estimates.  Values are in pixel units with a lower
%   bound of 0 pixels.
% <rfsize> contains pRF size estimates.  pRF size is defined as sigma/sqrt(n) where
%   sigma is the standard of the 2D Gaussian and n is the exponent of the power-law
%   function.  Values are in pixel units with a lower bound of 0 pixels.
% <expt> contains pRF exponent estimates.
% <gain> contains pRF gain estimates.  Values are in the same units of the data
%   and are constrained to be non-negative.
% <R2> contains R^2 values that indicate the goodness-of-fit of the model to the data.
%   Values are in percentages and generally range between 0% and 100%.  The R^2 values
%   are computed after projecting out polynomials from both the data and the model fit.
%   (Because of this projection, R^2 values can sometimes drop below 0%.)  Note that
%   if cross-validation is used (see <xvalmode>), the interpretation of <R2> changes
%   accordingly.
% <resnorms> and <numiters> contain optimization details (residual norms and 
%   number of iterations, respectively).
% <meanvol> contains the mean volume, that is, the mean of each voxel's time-series.
% <noisereg> contains a record of the noise regressors used in the model.
% <params> contains a record of the raw parameter estimates that are obtained internally
%   in the code.  These raw parameters are transformed to a more palatable format for
%   the user (as described above).
% <options> contains a record of the options used in the call to analyzePRF.m.
%
% Details on the pRF model:
% - Before analysis, we zero out any voxel that has a non-finite value or has all zeros
%   in at least one of the runs.  This prevents weird issues due to missing or bad data.
% - The pRF model that is fit is similar to that described in Dumoulin and Wandell (2008),
%   except that a static power-law nonlinearity is added to the model.  This new model, 
%   called the Compressive Spatial Summation (CSS) model, is described in Kay, Winawer, 
%   Mezer, & Wandell (2013).
% - The model involves computing the dot-product between the stimulus and a 2D isotropic
%   Gaussian, raising the result to an exponent, scaling the result by a gain factor,
%   and then convolving the result with a hemodynamic response function (HRF).  Polynomial
%   terms are included (on a run-by-run basis) to model the baseline signal level.
% - The 2D isotropic Gaussian is scaled such that the summation of the values in the
%   Gaussian is equal to one.  This eases the interpretation of the gain of the model.
% - The exponent parameter in the model is constrained to be non-negative.
% - The gain factor in the model is constrained to be non-negative; this aids the 
%   interpretation of the model (e.g. helps avoid voxels with negative BOLD responses
%   to the stimuli).
% - The workhorse of the analysis is fitnonlinearmodel.m, which is essentially a wrapper 
%   around routines in the MATLAB Optimization Toolbox.  We use the Levenberg-Marquardt 
%   algorithm for optimization, minimizing squared error between the model and the data.
% - A two-stage optimization strategy is used whereby all parameters excluding the
%   exponent parameter are first optimized (holding the exponent parameter fixed) and 
%   then all parameters are optimized (including the exponent parameter).  This 
%   strategy helps avoid local minima.
%
% Regarding GLMdenoise:
% - If the <wantglmdenoise> option is specified, we derive noise regressors using
%   GLMdenoise prior to model fitting.  This is done by creating a GLM design matrix
%   based on the contents of <stimulus> and then using this design matrix in conjunction
%   with GLMdenoise to analyze the data.  The noise regressors identified by GLMdenoise
%   are then used in the fitting of the pRF models (the regressors enter the model
%   additively, just like the polynomial regressors).
%
% Regarding seeding issues:
% - To minimize the impact of local minima, the default strategy is to perform full 
%   optimizations starting from three different initial seeds.
% - The first seed is a generic large pRF that is centered with respect to the stimulus,
%   has a pRF size equal to 1/4th of the stimulus extent (thus, +/- 2 pRF sizes matches
%   the stimulus extent), and has an exponent of 0.5.
% - The second seed is a generic small pRF that is just like the first seed except has
%   a pRF size that is 10 times smaller.
% - The third seed is a "supergrid" seed that is identified by performing a quick grid
%   search prior to optimization (similar in spirit to methods described in Dumoulin and 
%   Wandell, 2008).  In this procedure, a list of potential seeds is constructed by 
%   exploring a range of eccentricities, angles, and exponents.  For each potential 
%   seed, the model prediction is computed, and the seed that produces the closest 
%   match to the data is identified.  Note that the supergrid seed may be different
%   for different voxels.
%
% history:
% 2014/06/17 - version 1.1
% 2014/06/15 - add inputs <hrf> and <maxpolydeg>.
% 2014/06/10 - version 1.0
% 2014/04/27 - gain seed is now set to 0; add gain to the output
% 2014/04/29 - use typicalgain now (default 10). allow display input.

% internal notes:
% - for cluster mode, need to make sure fitnonlinearmodel is compiled (compilemcc.m)
% - to check whether local minima are a problem, can look at results.resnorms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REPORT

fprintf('*** analyzePRF: started at %s. ***\n',datestr(now));
stime = clock;  % start time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTERNAL CONSTANTS

% define
remotedir = '/scratch/knk/input/';
remotedir2 = '/scratch/knk/output/';
remotelogin = 'knk@login2.chpc.wustl.edu';
remoteuser = 'knk';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP AND PREPARATION

% massage cell inputs
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

% deal with inputs
if ~exist('options','var') || isempty(options)
  options = struct();
end
if ~isfield(options,'vxs') || isempty(options.vxs)
  options.vxs = 1:numvxs;
end
if ~isfield(options,'wantglmdenoise') || isempty(options.wantglmdenoise)
  options.wantglmdenoise = 0;
end
if ~isfield(options,'hrf') || isempty(options.hrf)
  options.hrf = [];
end
if ~isfield(options,'maxpolydeg') || isempty(options.maxpolydeg)
  options.maxpolydeg = [];
end
if ~isfield(options,'seedmode') || isempty(options.seedmode)
  options.seedmode = [0 1 2];
end
if ~isfield(options,'xvalmode') || isempty(options.xvalmode)
  options.xvalmode = 0;
end
if ~isfield(options,'numperjob') || isempty(options.numperjob)
  options.numperjob = [];
end
if ~isfield(options,'maxiter') || isempty(options.maxiter)
  options.maxiter = 500;
end
if ~isfield(options,'display') || isempty(options.display)
  options.display = 'iter';
end
if ~isfield(options,'typicalgain') || isempty(options.typicalgain)
  options.typicalgain = 10;
end

% massage
options.seedmode = union(options.seedmode(:),[]);

% calc
usecluster = ~isempty(options.numperjob);

% prepare stimuli
for p=1:length(stimulus)
  stimulus{p} = squish(stimulus{p},2)';  % frames x pixels
  stimulus{p} = [stimulus{p} p*ones(size(stimulus{p},1),1)];  % add a dummy column to indicate run breaks
  stimulus{p} = single(stimulus{p});  % make single to save memory
end

% deal with data badness (set bad voxels to be always all 0)
bad = cellfun(@(x) any(~isfinite(x),dimtime) | all(x==0,dimtime),data,'UniformOutput',0);  % if non-finite or all 0
bad = any(cat(dimtime,bad{:}),dimtime);  % badness in ANY run
for p=1:numruns
  data{p}(repmat(bad,[ones(1,dimdata) size(data{p},dimtime)])) = 0;
end

% calc mean volume
meanvol = mean(catcell(dimtime,data),dimtime);

% what HRF should we use?
if isempty(options.hrf)
  options.hrf = getcanonicalhrf(tr,tr)';
end
numinhrf = length(options.hrf);

% what polynomials should we use?
if isempty(options.maxpolydeg)
  options.maxpolydeg = cellfun(@(x) round(size(x,dimtime)*tr/60/2),data);
end
if isscalar(options.maxpolydeg)
  options.maxpolydeg = repmat(options.maxpolydeg,[1 numruns]);
end
fprintf('using the following maximum polynomial degrees: %s\n',mat2str(options.maxpolydeg));

% initialize cluster stuff
if usecluster
  localfilestodelete = {};
  remotefilestodelete = {};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE OUT NOISE REGRESSORS

if isequal(options.wantglmdenoise,1)
  noisereg = analyzePRFcomputeGLMdenoiseregressors(stimulus,data,tr);
elseif isequal(options.wantglmdenoise,0)
  noisereg = [];
else
  noisereg = options.wantglmdenoise;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE MODEL

% pre-compute some cache
[d,xx,yy] = makegaussian2d(resmx,2,2,2,2);

% define the model (parameters are R C S G N)
modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),options.hrf,dd(:,prod(res)+1));
model = {{[] [1-res(1)+1 1-res(2)+1 0    0   NaN;
              2*res(1)-1 2*res(2)-1 Inf  Inf Inf] modelfun} ...
         {@(ss)ss [1-res(1)+1 1-res(2)+1 0    0   0;
                   2*res(1)-1 2*res(2)-1 Inf  Inf Inf] @(ss)modelfun}};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE SEEDS

% init
seeds = [];

% generic large seed
if ismember(0,options.seedmode)
  seeds = [seeds;
           (1+res(1))/2 (1+res(2))/2 resmx/4*sqrt(0.5) options.typicalgain 0.5];
end

% generic small seed
if ismember(1,options.seedmode)
  seeds = [seeds;
           (1+res(1))/2 (1+res(2))/2 resmx/4*sqrt(0.5)/10 options.typicalgain 0.5];
end

% super-grid seed
if ismember(2,options.seedmode)
  supergridseeds = analyzePRFcomputesupergridseeds(res,stimulus,data,modelfun, ...
                                                   options.maxpolydeg,dimdata,dimtime, ...
                                                   options.typicalgain);
end

% make a function that individualizes the seeds
if exist('supergridseeds','var')
  seedfun = @(vx) [[seeds];
                   [subscript(squish(supergridseeds,dimdata),{vx ':'})]];
else
  seedfun = @(vx) [seeds];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE RESAMPLING STUFF

% define wantresampleruns and resampling
switch options.xvalmode
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

  % define stimulus
  stimulus = @() loadmulti(remotefile0,'stimulus');

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

  % define data
  binfiles = cellfun(@(x) [remotefile0 sprintf('/%03d.bin',x)],num2cell(1:numruns),'UniformOutput',0);
  data = @(vxs) cellfun(@(x) double(loadbinary(x,'single',[0 numvxs],-vxs)),binfiles,'UniformOutput',0);

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
  outputdir = outputdirremote;

%%%%% NON-CLUSTER CASE

else

  stimulus = {stimulus};
  data = @(vxs) cellfun(@(x) subscript(squish(x,dimdata),{vxs ':'})',data,'UniformOutput',0);
  outputdir = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE OPTIONS

% last-minute prep
if iscell(noisereg)
  noiseregINPUT = {noisereg};
else
  noiseregINPUT = noisereg;
end

% construct the options struct
opt = struct( ...
  'outputdir',outputdir, ...
  'stimulus',stimulus, ...
  'data',data, ...
  'vxs',options.vxs, ...
  'model',{model}, ...
  'seed',seedfun, ...
  'optimoptions',{{'Display' options.display 'Algorithm' 'levenberg-marquardt' 'MaxIter' options.maxiter}}, ...
  'wantresampleruns',wantresampleruns, ...
  'resampling',resampling, ...
  'metric',@calccod, ...
  'maxpolydeg',options.maxpolydeg, ...
  'wantremovepoly',1, ...
  'extraregressors',noiseregINPUT, ...
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
  jobids = [jobids chpcrun(jobnames{end},'fitnonlinearmodel',options.numperjob, ...
                           1,ceil(length(options.vxs)/options.numperjob),[], ...
                           {'data' 'stimulus' 'bad' 'd' 'xx' 'yy' 'modelfun' 'model'})];

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

% massage model parameters for output and put in 'results' struct
results.ang(options.vxs,:) =    permute(mod(atan2((1+res(1))/2 - a1.params(:,1,:), ...
                                                  a1.params(:,2,:) - (1+res(2))/2),2*pi)/pi*180,[3 1 2]);
results.ecc(options.vxs,:) =    permute(sqrt(((1+res(1))/2 - a1.params(:,1,:)).^2 + ...
                                             (a1.params(:,2,:) - (1+res(2))/2).^2),[3 1 2]);
results.expt(options.vxs,:) =   permute(posrect(a1.params(:,5,:)),[3 1 2]);
results.rfsize(options.vxs,:) = permute(abs(a1.params(:,3,:)) ./ sqrt(posrect(a1.params(:,5,:))),[3 1 2]);
results.R2(options.vxs,:) =     permute(a1.trainperformance,[2 1]);
results.gain(options.vxs,:) =   permute(posrect(a1.params(:,4,:)),[3 1 2]);
results.resnorms(options.vxs) = a1.resnorms;
results.numiters(options.vxs) = a1.numiters;

% reshape
results.ang =      reshape(results.ang,      [xyzsize numfits]);
results.ecc =      reshape(results.ecc,      [xyzsize numfits]);
results.expt =     reshape(results.expt,     [xyzsize numfits]);
results.rfsize =   reshape(results.rfsize,   [xyzsize numfits]);
results.R2 =       reshape(results.R2,       [xyzsize numfits]);
results.gain =     reshape(results.gain,     [xyzsize numfits]);
results.resnorms = reshape(results.resnorms, [xyzsize 1]);
results.numiters = reshape(results.numiters, [xyzsize 1]);

% add some more stuff
results.meanvol =  meanvol;
results.noisereg = noisereg;
results.params =   a1.params;
results.options = options;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUNK

% % define the model (parameters are R C S G N [HRF])
% modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),pp(5+(1:numinhrf))',dd(:,prod(res)+1));
% model = {{[] [1-res(1)+1 1-res(2)+1 0    0   NaN repmat(NaN,[1 numinhrf]);
%               2*res(1)-1 2*res(2)-1 Inf  Inf Inf repmat(Inf,[1 numinhrf])] modelfun} ...
%          {@(ss)ss [1-res(1)+1 1-res(2)+1 0    0   0   repmat(NaN,[1 numinhrf]);
%                    2*res(1)-1 2*res(2)-1 Inf  Inf Inf repmat(Inf,[1 numinhrf])] @(ss)modelfun} ...
%          {@(ss)ss [1-res(1)+1 1-res(2)+1 0    0   0   repmat(-Inf,[1 numinhrf]);
%                    2*res(1)-1 2*res(2)-1 Inf  Inf Inf repmat(Inf,[1 numinhrf])] @(ss)modelfun}};
% 
% % if not fitting the HRF, exclude the last model step
% if ~wantfithrf
%   model = model(1:2);
% end
%wantfithrf = 0;    % for now, leave at 0

% results.hrf =      NaN*zeros(numvxs,numinhrf,numfits);
% results.hrf(options.vxs,:,:) =  permute(a1.params(:,5+(1:numinhrf),:),[3 2 1]);
% results.hrf =      reshape(results.hrf,      [xyzsize numinhrf numfits]);
