classdef pRF_timeShift < handle
    
    properties (Constant)
        nParams = 6;
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
            
            % Create and cache the 2D Gaussian in a private property
            [~,obj.xx,obj.yy] = makegaussian2d(max(res),2,2,2,2);
        end
        
        % Methods
        x0 = initial(obj)
        [lb, ub] = bounds(obj)
        signal = forward(obj, params)
        
    end
end