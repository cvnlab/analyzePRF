function noisereg = analyzePRFcomputeGLMdenoiseregressors(stimulus,data,tr)

% function noisereg = analyzePRFcomputeGLMdenoiseregressors(stimulus,data,tr)
%
% <stimulus>,<data>,<tr> are from analyzePRF
%
% this is an internal function called by analyzePRF.m.  this function constructs a 
% GLM design matrix based on <stimulus> and then uses GLMdenoisedata.m to determine 
% the optimal set of noise regressors.  these regressors are returned in <noisereg>.
%
% in case you want to see the results of GLMdenoisedata.m, the figure output and 
% the results of GLMdenoisedata.m are saved to temporary files (as reported to the
% command window).

% put up a warning
warning('use of GLMdenoiseregressors in analyzePRF is experimental and not recommended');

% internal constants
corrthresh = .9;  % used in determining which apertures are the same

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figure out GLM design matrix

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
  spots = find(corrs > corrthresh);  % any frame with r > corrthresh counts as the same
%%%  fprintf('cnt=%d: numspots=%d\n',cnt,length(spots));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% run GLMdenoise to get the noise regressors

% what directory to save results to?
glmdenoisedir = [tempname];
assert(mkdir(glmdenoisedir));

% call GLMdenoise
fprintf('using GLMdenoise figure directory %s\n',[glmdenoisedir '/GLMdenoisefigures']);

% prep
hrfmodel = 'assume';
hrfknobs = [];

% run it
results = GLMdenoisedata(glmdesign,data,tr,tr,hrfmodel,hrfknobs, ...
                         struct('numboots',0), ...
                         [glmdenoisedir '/GLMdenoisefigures']);

% get the noise regressors
noisereg = cellfun(@(x) x(:,1:results.pcnum),results.pcregressors,'UniformOutput',0);

% save 'results' to a file
file0 = [glmdenoisedir '/GLMdenoise.mat'];
fprintf('saving GLMdenoise results to %s (in case you want them).\n',file0);
results = rmfield(results,{'pcweights' 'models' 'modelse'});  % remove some boring fields
save(file0,'results','noisereg','glmdesign');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% quick visualization

for p=1:length(glmdesign)
  figureprep([100 100 700 700]); hold on;
  imagesc(glmdesign{p}); colormap(gray);
  axis image tight;
  set(gca,'YDir','reverse');
  title(sprintf('run %02d',p));
  figurewrite(sprintf('glmdesign%02d',p),[],[],glmdenoisedir);
end
