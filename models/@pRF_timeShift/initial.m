function x0 = initial(obj)
% Initial guess for the model
%
% Inputs:
%   res                   - Scalar. the resolution of the stimulus in
%                           pixels (typically ~100).
%   gain                  - The typical gain (i.e., aplitude) of the
%                           response to the stimulus modulation. In raw
%                           BOLD fMRI data this value maybe 50. If the data
%                           are expessed in percentage change terms, then
%                           a value of 1 is appropriate, and if in terms of
%                           proportion response, then 0.01.
%   scale                 - Char vector. Indicates how big to make the
%                           receptive fields.
%                           

res = obj.res;
gain = obj.gain;
scale = obj.seedScale;
nParams = obj.nParams;

x0 = zeros(1,nParams);

% Handle incomplete arguments
if nargin==0
    return
end

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
x0(4) = gain;                       % typical gain (amplitude)
x0(5) = 0.5;                        % compressive exponent
x0(6) = 0;                          % HRF temporal shift (in TRs)

end

