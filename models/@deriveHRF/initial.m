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
typicalGain = typicalGain;

% Assign the x0 variable
x0 = zeros(1,nParams);

% Assemble X0
x0(1) = 4;               % xPosition
x0(2) = 8;               % yPosition
x0(3) = 10;              % sigma
x0(4) = typicalGain;     % typical gain (amplitude)

end

