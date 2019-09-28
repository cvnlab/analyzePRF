function X0 = modelX0(res,resmx,gain,scale)
% Initial guess for the model
%
% Inputs:
%   res                   - Scalar. the resolution of the stimulus in
%                           pixels (typically ~100).
%   resmx
%   gain                  - The typical gain (i.e., aplitude) of the
%                           response to the stimulus modulation. In raw
%                           BOLD fMRI data this value maybe 50. If the data
%                           are expessed in percentage change terms, then
%                           a value of 1 is appropriate, and if in terms of
%                           proportion response, then 0.01.
%   scale                 - Char vector. Indiates how big to make the
%                           receptive fields.
%                           

nParams = 7;
X0 = zeros(1,nParams);

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
X0(1) = (1+res(1))/2;               % xPosition
X0(2) = (1+res(2))/2;               % yPosition
X0(3) = resmx/4*sqrt(0.5)/scaler;   % sigma
X0(4) = gain;                       % typical gain (amplitude)
X0(5) = 0.5;                        % compressive exponent
X0(6) = 0;                          % HRF temporal shift (in TRs)
X0(7) = 0;                          % unused

end

