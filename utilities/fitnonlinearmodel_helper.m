function results = fitnonlinearmodel_helper(opt,stimulus)

% This is a helper function for fitnonlinearmodel.m.  Not for external use!
%
% Notes:
% - opt.data is always a cell vector and contains only one voxel
% - in the nonlinear case, the seed to use has been hacked into model{1}{1} and may have multiple rows


% Set up the seeds
ismultipleseeds = size(opt.model{1}{1},1) > 1;
ismultiplemodels = length(opt.model) > 1;


% calc

numparams = size(opt.model{end}{2},2);

data = opt.data{1}{1};

% init
results = struct;
results.params = zeros(1,numparams);
results.modelpred = cell(1,1);  % but converted to a matrix at the end
results.modelfit =  cell(1,1);  % but converted to a matrix at the end



% loop over resampling cases

datastd = std(data);
if datastd == 0
    datastd = 1;
end
data = data / datastd;


options = opt.optimoptions;
if ~isempty(opt.outputfcn)
    if nargin(opt.outputfcn) == 4
        options.OutputFcn = @(a,b,c) feval(opt.outputfcn,a,b,c,data);
    else
        options.OutputFcn = opt.outputfcn;
    end
end

% loop over seeds
params = [];
for ss=1:size(opt.model{1}{1},1)
    
    % loop through models
    for mm=1:length(opt.model)
        
        % which parameters are we actually fitting?
        ix = ~isnan(opt.model{mm}{2}(1,:));
        
        % calculate seed, model, and transform
        if mm==1
            seed = opt.model{mm}{1}(ss,:);
            model = opt.model{mm}{3};
            transform = opt.model{mm}{4};
        else
            seed = feval(opt.model{mm}{1},params0);
            model = feval(opt.model{mm}{3},params0);
            transform = feval(opt.model{mm}{4},params0);
        end
        
        % figure out bounds to use
        if isequal(options.Algorithm,'levenberg-marquardt')
            lb = [];
            ub = [];
        else
            lb = opt.model{mm}{2}(1,ix);
            ub = opt.model{mm}{2}(2,ix);
        end
        
        % define the final model function
        fun = @(pp) feval(model,copymatrix(seed,ix,pp));
        
        % perform the fit (NOTICE THE DIVISION BY DATASTD, THE NAN PROTECTION, THE CONVERSION TO DOUBLE)
        if ~any(ix)
            params0 = seed;   % if no parameters are to be optimized, just return the seed
            resnorm = NaN;
            output = [];
            output.iterations = NaN;
        else
            [params0,resnorm,residual,exitflag,output] = ...
                lsqcurvefit(@(x,y) double(nanreplace(feval(fun,x) / datastd,0,2)),seed(ix),[],double(data),lb,ub,options);
            params0 = copymatrix(seed,ix,params0);
        end
        
        % record
        results.numiters(rr,ss,mm) = output.iterations;
        
    end
    
    % record
    results.resnorms(rr,ss) = resnorm;  % final resnorm
    params(ss,:) = params0;  % final parameters
    
end

% which seed produced the best results?
[d,mnix] = min(results.resnorms(rr,:));
finalparams = params(mnix,:);


% record the results
results.params(rr,:) = finalparams;

% report
if ~islinear && ismultipleseeds
    if strcmp(verbosity,'all')
        fprintf('    seed %d was best. final estimated parameters are [',mnix); ...
            fprintf('%.3f ',finalparams); fprintf('].\n');
    end
end

% prepare data and model fits
% [NOTE: in the nonlinear case, this inherits model, transform, and trainstimTRANSFORM from above!!]
traindatatemp = trainS*traindata;
if islinear
    modelfittemp = trainS*(trainstim*finalparams');
else
    modelfittemp = nanreplace(trainS*feval(model,finalparams,trainstimTRANSFORM),0,2);
end

if isempty(testdata)  % handle this case explicitly, just to avoid problems
    results.testdata{rr} = [];
    results.modelpred{rr} = [];
else
    results.testdata{rr} = testS*testdata;
    if islinear
        results.modelpred{rr} = testS*(teststim*finalparams');
    else
        results.modelpred{rr} = nanreplace(testS*feval(model,finalparams,feval(transform,teststim)),0,2);
    end
end

% prepare modelfit
if wantmodelfit
    if islinear
        results.modelfit{rr} = (allstim*finalparams')';
    else
        results.modelfit{rr} = nanreplace(feval(model,finalparams,feval(transform,allstim)),0,2)';
    end
else
    results.modelfit{rr} = [];  % if not wanted by user, don't bother computing
end

% compute metrics
results.trainperformance(rr) = feval(opt.metric,modelfittemp,traindatatemp);
if isempty(results.testdata{rr})  % handle this case explicitly, just to avoid problems
    results.testperformance(rr) = NaN;
else
    results.testperformance(rr) = feval(opt.metric,results.modelpred{rr},results.testdata{rr});
end


% compute aggregated metrics
results.testdata = catcell(1,results.testdata);
results.modelpred = catcell(1,results.modelpred);
results.modelfit = catcell(1,results.modelfit);
if isempty(results.testdata)
    results.aggregatedtestperformance = NaN;
else
    results.aggregatedtestperformance = feval(opt.metric,results.modelpred,results.testdata);
end
