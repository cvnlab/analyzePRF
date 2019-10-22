classdef deriveHRF < handle
    
    properties (Constant)
        
        % The identity of the dimensions of the data variable
        dimdata = 1;
        dimtime = 2;
        
        % THe number of parameters in the model
        nParams = 4;
        
        % The model is executed as a single stage search.
        nStages = 1;
        floatSet = {[1 2 3 4 5 6]};
        fixSet = {[]};
        
        % A description of the model
        description = ...
            ['The deriveHRF model'];
    end
    
    % Private properties
    properties (GetAccess=private)
        % The projection matrix used to regress our nuisance effects
        T
    end
    
    % Seen, but not touched
    properties (SetAccess=private)
        
        % The stimulus vector, concatenated across acquisitions. Thus it
        % will have the dimensions totalTRs x 1
        stimulus
        
        % A vector of the length totalTRs x 1 that has an index value to
        % indicate which acquisition (1, 2, 3 ...) this TR is from.
        acqGroups
        
        % TR of the data in seconds
        tr
        
        % A cell array that contains things that the model might want
        payload
        
        % The number of acquisitions
        nAcqs
        
        % A vector with the number of TRs in each acquisition.
        nTRsPerAcq

    end
    
    % These may be modified after object creation
    properties (SetAccess=public)
        
        % The number of low frequencies to be removed from each acquisition
        polyDeg
        
        % Typical amplitude of the BOLD fMRI response in the data
        typicalGain
        
        % HRF duration to model in seconds
        duration
        
        % Verbosity
        verbose
    end
    
    methods

        % Constructor
        function obj = deriveHRF(data,stimulus,tr,varargin)
                        
            % instantiate input parser
            p = inputParser; p.KeepUnmatched = false;
            
            % Required
            p.addRequired('data',@iscell);
            p.addRequired('stimulus',@iscell);
            p.addRequired('tr',@isscalar);
            
            p.addParameter('payload',{},@iscell);
            p.addParameter('polyDeg',[],@isscalar);
            p.addParameter('typicalGain',300,@isscalar);
            p.addParameter('duration',50,@isscalar);
            p.addParameter('verbose',true,@islogical);

            % parse
            p.parse(data, stimulus, tr, varargin{:})
            
            % Derive properties from the data variable and then clear
            obj.nAcqs = length(data);
            obj.nTRsPerAcq = cellfun(@(x) size(x,2),data);
            clear data
                        
            % Vectorize the stimuli. Add a dummy column to indicate run
            % breaks. Concatenate the cells and store
            for ii=1:length(stimulus)
                stimulus{ii} = squish(stimulus{ii},2)';
                acqGroups{ii} = ii*ones(size(stimulus{ii},1),1);
            end
            obj.stimulus = catcell(1,stimulus);
            obj.acqGroups = catcell(1,acqGroups);
            clear stimulus acqGroups
            
            % Distribute other params to obj properties
            obj.tr = tr;
            obj.payload = p.Results.payload;
            obj.polyDeg = p.Results.polyDeg;
            obj.typicalGain = p.Results.typicalGain;
            obj.duration = p.Results.duration;
            obj.verbose = p.Results.verbose;
            
            % Create and cache the projection matrix
            obj.genprojection;
                        
        end
        
        % Set methods
        function set.polyDeg(obj, value)
            obj.polyDeg = value;
            obj.genprojection;
        end
        
        % Methods
        rawData = prep(obj,rawData)
        genprojection(obj)
        x0 = initial(obj)
        [lb, ub] = bounds(obj)
        signal = clean(obj, signal)
        [fit, hrf] = forward(obj, pp)
        metric = metric(obj, signal, x)
        seeds = seeds(obj, data, vxs)
        results = results(obj, params, metric)
    end
end