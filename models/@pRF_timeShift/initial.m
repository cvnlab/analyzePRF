function x0 = initial(obj)
% Returns initial guess for the model parameters
%
% Syntax:
%   x0 = obj.initial()
%
% Description:
%   Initial values for the prf_timeShift model. Rationale is as follows:
%       x, y :  Center of the stimulus
%       sigma:  1 or 10 pixels, depending upon obj.scale
%       gain :  Set by obj.typicalGain
%       exp  :  Locked to 0.05, following Benson et al, 2018, HCP 7T data
%       shift:  Zero HRF temporal shift      
%
% Inputs:
%   none
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   x0                    - 1xnParams vector.
%


% Obj variables
res = obj.res;
typicalGain = obj.typicalGain;
scale = obj.seedScale;
nParams = obj.nParams;

% Assign the x0 variable
x0 = zeros(1,nParams);

% Adjust the size of the receptive field sigma
switch scale
    case 'large'
        scaler = 1;
    case 'small'
        scaler = 10;
    otherwise
        scaler = 5;
end

% Assemble X0
x0(1) = (1+res(1))/2;               % xPosition
x0(2) = (1+res(2))/2;               % yPosition
x0(3) = max(res)/4*sqrt(0.5)/scaler;   % sigma
x0(4) = typicalGain;                % typical gain (amplitude)
x0(5) = 0.05;                       % compressive exponent
x0(6) = 0;                          % HRF temporal shift (in TRs)

end

