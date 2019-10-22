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
x0(1) = 6;               % Aratio
x0(2) = 7;               % alpha1 
x0(3) = 1;               % beta1
x0(4) = 16;              % alpha2
x0(5) = 1;               % beta2
x0(6) = typicalGain;     % typical gain (amplitude)

end

