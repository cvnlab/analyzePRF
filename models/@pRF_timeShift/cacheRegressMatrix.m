function cacheRegressMatrix(obj)

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

% construct total regressors matrix for fitting purposes
% (i.e. both polynomials and extra regressors)
% first, construct the run-wise regressors
tmatrix = {};
for pp=1:nAcqs
    tmatrix{pp} = cat(2,polyregressors{pp});
end
% then, separate them using blkdiag
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
