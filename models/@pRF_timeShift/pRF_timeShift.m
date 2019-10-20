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
        T
    end
    
    % Fixed after object creation
    properties (SetAccess=private)
        stimulus
        res
        hrf
        tr
        payload
        nAcqs
        nTRsPerAcq
    end
    
    % These may be modified after object creation
    properties (SetAccess=public)
        seedScale
        typicalGain
        forceBounds
        verbose
    end
    
    methods

        % Constructor
        function obj = pRF_timeShift(data,stimulus,res,hrf,tr,varargin)
                        
            % instantiate input parser
            p = inputParser; p.KeepUnmatched = false;
            
            % Required
            p.addRequired('data',@iscell);
            p.addRequired('stimulus',@iscell);
            p.addRequired('res',@isvector);
            p.addRequired('hrf',@isvector);
            p.addRequired('tr',@isscalar);
            
            p.addParameter('payload',{},@iscell);
            p.addParameter('typicalGain',30,@isscalar);
            p.addParameter('forceBounds',true,@islogical);
            p.addParameter('verbose',true,@islogical);
            
            % parse
            p.parse(data,stimulus,res,hrf,tr, varargin{:})
            
            % Derive properties from the data variable and then clear
            obj.nAcqs = length(data);
            obj.nTRsPerAcq = cellfun(@(x) size(x,2),data);
            clear data
            
            % Distribute passed params to obj properties
            obj.stimulus = catcell(1,stimulus);
            obj.res = res;
            obj.hrf = hrf;
            obj.tr = tr;
            obj.payload = p.Results.payload;
            
            % Set defaults
            obj.typicalGain = p.Results.typicalGain;
            obj.seedScale = 'medium';
            obj.forceBounds = p.Results.forceBounds;
            obj.verbose = p.Results.verbose;

            % Create and cache the projection matrix
            obj.cacheProjectionMatrix;
            
            % Create and cache the 2D Gaussian in a private property
            [~,obj.xx,obj.yy] = makegaussian2d(max(res),2,2,2,2);
            
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