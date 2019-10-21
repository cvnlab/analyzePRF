function y = fshift(x,s)
% FSHIFT Fractional circular shift
%   Syntax:
%
%       >> y = fshift(x,s)
%
%   FSHIFT circularly shifts the elements of vector x by a (possibly
%   non-integer) number of elements s. FSHIFT works by applying a linear
%   phase in the spectrum domain and is equivalent to CIRCSHIFT for integer
%   values of argument s (to machine precision).

% Author:   François Bouffard
%           fbouffard@gmail.com

needtr = 0; 
if size(x,1) == 1; 
    x = x(:); 
    needtr = 1; 
end;
N = size(x,1); 
r = floor(N/2)+1; 
f = ((1:N)-r)/(N/2); 
p = exp(-1j*s*pi*f).';
if ~mod(N,2)
    % N is even. This becomes important for complex signals.
    % Thanks to Ahmed Fasih for pointing out the bug.
    % For even signals, f(1) = -1 and phase is sampled at -pi. 
    % The correct value for p(1) should be the average of the f = -1 and
    % f = 1 cases. Since f has antisymmetric phase and unitary (and thus
    % symmetric) magnitude, the average is the real part of p.
    p(1) = real(p(1));
end
y = ifft(fft(x).*ifftshift(p)); 
if isreal(x);
    y = real(y); 
end;
if needtr; 
    y = y.'; 
end;