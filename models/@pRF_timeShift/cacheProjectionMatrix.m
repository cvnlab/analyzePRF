function cacheProjectionMatrix(obj)
% Generates and stores a projection matrix that removes nuisance effects
%
% Syntax:
%   obj.cacheProjectionMatrix
%
% Description:
%   Generates a set of polynomial regressors and a corresponding
%   partialling matrix. Matrix multiplication of the projection matrix by
%   the time-series results in a residual time-series that has had the
%   effects modeled in the projection matrix removed.
%
% Inputs:
%   none
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   none
%


% Obj variables
tr = obj.tr;
nAcqs = obj.nAcqs;
nTRsPerAcq = obj.nTRsPerAcq;

% Decide how many low frequencies to filter from each acquisition
maxpolydeg = round(nTRsPerAcq.*tr/60/2);

% Construct polynomial regressors matrix
polyregressors = {};
for pp=1:nAcqs
    if isnan(maxpolydeg(pp))
        polyregressors{pp} = zeros(nTRsPerAcq(pp),0);
    else
        polyregressors{pp} = constructpolynomialmatrix(nTRsPerAcq(pp),0:maxpolydeg(pp));
    end
end

% Construct total regressors matrix for fitting purposes. This is the point
% at which extra regressors could be added. First generate the entire
% matrix.
tmatrix = {};
for pp=1:nAcqs
    tmatrix{pp} = cat(2,polyregressors{pp});
end
% Then, separate them using blkdiag
temp = blkdiag(tmatrix{:});
cnt = 0;
for pp=1:nAcqs
    tmatrix{pp} = temp(cnt+(1:size(tmatrix{pp},1)),:);
    cnt = cnt + size(tmatrix{pp},1);
end
clear temp;

% Store the projection matrix
obj.T = projectionmatrix(catcell(1,tmatrix));

end
