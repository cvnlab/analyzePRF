function BOLDHRF = createCanonicalHRF(t,gamma1,gamma2,gammaScale)

% function HRF = createCanonicalHRF(gamma1,gamma2,gammaScale)
%
% creates canonical HRF

BOLDHRF = gampdf(t, gamma1, 1) - gampdf(t, gamma2, 1)/gammaScale ;
% scale to unit sum to preserve amplitude of y following convolution
BOLDHRF = BOLDHRF/sum(BOLDHRF) ;

gribble = 1;