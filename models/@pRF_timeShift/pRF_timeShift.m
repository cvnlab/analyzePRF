classdef pRF_timeShift < handle
    
    properties (Constant)
        nParams = 6;
        nStages = 2;
        floatSet = {[1 2 3 4],[1 2 3 4 5 6]};
        fixSet = {[5 6],[]};
    end
    
    % Private properties
    properties (GetAccess=private)
        xx
        yy
    end
    
    % Fixed after object creation
    properties (SetAccess=private)
        stimulus
        res
        hrf
    end
    
    % These may be modified after object creation
    properties (SetAccess=public)
        fixed
        seedScale
        gain
        verbose
    end
    
    methods
        % Constructor
        function obj = pRF_timeShift(stimulus,res,hrf)
            % Add some parsing in here
            obj.stimulus = stimulus;
            obj.res = res;
            obj.hrf = hrf;
            
            % Set by default
            obj.gain = 30;
            obj.fixed = [];
            obj.seedScale = 'large';
            obj.verbose = true;
            
            % Create and cache the 2D Gaussian in a private property
            [~,obj.xx,obj.yy] = makegaussian2d(max(res),2,2,2,2);
        end
        
        % Methods
        x0 = initial(obj)
        [lb, ub] = bounds(obj)
        signal = prep(obj, signal)
        fit = forward(obj, params)
        metric = metric(obj, signal, x)
        seeds = seeds(obj, data, vxs)
        results = results(obj, params, metric)
    end
end