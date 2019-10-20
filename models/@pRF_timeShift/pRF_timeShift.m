classdef pRF_timeShift < handle
    
    properties (Constant)
        dimdata = 1;
        dimtime = 2;
        nParams = 6;
        nStages = 2;
        floatSet = {[1 2 4],[1 2 3 4 6]};
        fixSet = {[3 5 6],[5]};
        description = ...
            ['A pRF mapping approach that assumes a circular symmetric \n' ...
             'Gaussian receptive field and a fixed, compressive non- \n' ...
             'linearity. The model adjusts the time-to-peak of the HRF \n',...
             'and thus requires that the stimulus play in forward and \n',...
             'time-reversed directions. A two-stage non-linear search \n', ...
             'is performed, fist across pRF center and gain, and then \n', ...
             'across the entire parameter set. \n'];
    end
    
    % Private properties
    properties (GetAccess=private)        
        xx
        yy
        hrf
        T
    end
    
    % Fixed after object creation
    properties (SetAccess=private)
        stimulus
        res
        tr
        payload
        nAcqs
        nTRsPerAcq
    end
    
    % These may be modified after object creation
    properties (SetAccess=public)
        hrfParams
        seedScale
        typicalGain
        forceBounds
        verbose
    end
    
    methods

        % Constructor
        function obj = pRF_timeShift(data,stimulus,tr,varargin)
                        
            % instantiate input parser
            p = inputParser; p.KeepUnmatched = false;
            
            % Required
            p.addRequired('data',@iscell);
            p.addRequired('stimulus',@iscell);
            p.addRequired('tr',@isscalar);
            
            p.addParameter('payload',{},@iscell);
            p.addParameter('hrfParams',[6 12 10 20],@isvector);
            p.addParameter('typicalGain',30,@isscalar);
            p.addParameter('seedScale','medium',@ischar);
            p.addParameter('forceBounds',true,@islogical);
            p.addParameter('verbose',true,@islogical);

            % parse
            p.parse(data, stimulus, tr, varargin{:})
            
            % Derive properties from the data variable and then clear
            obj.nAcqs = length(data);
            obj.nTRsPerAcq = cellfun(@(x) size(x,2),data);
            clear data
            
            % Obtain the dimensions of the stimulus frames and store
            res = [size(stimulus{1},1) size(stimulus{1},2)];
            obj.res = res;
            
            % Vectorize the stimuli. Add a dummy column to indicate run
            % breaks. Concatenate the cells and store
            for ii=1:length(stimulus)
                stimulus{ii} = squish(stimulus{ii},2)';
                stimulus{ii} = [stimulus{ii} ii*ones(size(stimulus{ii},1),1)];
            end
            obj.stimulus = catcell(1,stimulus);
            clear stimulus
            
            % Distribute other params to obj properties
            obj.tr = tr;
            obj.payload = p.Results.payload;
            obj.hrfParams = p.Results.hrfParams;
            obj.typicalGain = p.Results.typicalGain;
            obj.seedScale = 'medium';
            obj.forceBounds = p.Results.forceBounds;
            obj.verbose = p.Results.verbose;

            % Create and cache the hrf
            obj.genhrf
            
            % Create and cache the projection matrix
            obj.cacheProjectionMatrix;
            
            % Create and cache the 2D Gaussian in a private property
            [~,obj.xx,obj.yy] = makegaussian2d(max(res),2,2,2,2);
            
        end
        
        % Set methods
        function set.hrfParams(obj, value)
            obj.hrfParams = value;
            obj.genhrf;
        end
        
        % Methods
        rawData = prep(obj,rawData)
        cacheProjectionMatrix(obj)
        x0 = initial(obj)
        [lb, ub] = bounds(obj)
        signal = clean(obj, signal)
        fit = forward(obj, params)
        metric = metric(obj, signal, x)
        seeds = seeds(obj, data, vxs)
        results = results(obj, params, metric)
    end
end