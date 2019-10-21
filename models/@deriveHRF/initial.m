function x0 = initial(obj)
% Returns initial guess for the model parameters
%
% Syntax:
%   x0 = obj.initial()
%
% Description:
%   Initial values for the deriveHRF model.
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
nParams = obj.nParams;
typicalGain = obj.typicalGain;

% Assign the x0 variable
x0 = zeros(1,nParams);

% Assemble X0
x0(1) = 4;               % gamma1
x0(2) = 2;               % gamma2 = gamma1 * x(2) 
x0(3) = 10;              % ratio of gamma 1 / gamma 2 amplitudes
x0(4) = typicalGain;     % typical gain (amplitude)

end

