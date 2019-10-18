
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE OUTPUT

% depending on which analysis we did (quick or full optimization),
% we have to get the outputs in a common format
if wantquick
    paramsA = permute(squish(supergridseeds,dimdata),[3 2 1]);  % fits x parameters x voxels
    rA = squish(rvalues,dimdata)';                              % fits x voxels
else
    paramsA = a1.params;                                        % fits x parameters x voxels
    rA = a1.trainperformance;                                   % fits x voxels
end

% calc
numfits = size(paramsA,1);

% init
clear results;
results.ang =      NaN*zeros(numvxs,numfits);
results.ecc =      NaN*zeros(numvxs,numfits);
results.expt =     NaN*zeros(numvxs,numfits);
results.rfsize =   NaN*zeros(numvxs,numfits);
results.R2 =       NaN*zeros(numvxs,numfits);
results.gain =     NaN*zeros(numvxs,numfits);
results.hrfshift = NaN*zeros(numvxs,numfits);
results.resnorms = cell(numvxs,1);
results.numiters = cell(numvxs,1);

% massage model parameters for output and put in 'results' struct
results.ang(options.vxs,:) =    permute(mod(atan2((1+res(1))/2 - paramsA(:,1,:), ...
    paramsA(:,2,:) - (1+res(2))/2),2*pi)/pi*180,[3 1 2]);
results.ecc(options.vxs,:) =    permute(sqrt(((1+res(1))/2 - paramsA(:,1,:)).^2 + ...
    (paramsA(:,2,:) - (1+res(2))/2).^2),[3 1 2]);
results.expt(options.vxs,:) =   permute(posrect(paramsA(:,5,:)),[3 1 2]);
results.hrfshift(options.vxs,:) =   permute(paramsA(:,6,:),[3 1 2]);
results.rfsize(options.vxs,:) = permute(abs(paramsA(:,3,:)) ./ sqrt(posrect(paramsA(:,5,:))),[3 1 2]);
results.R2(options.vxs,:) =     permute(rA,[2 1]);
results.gain(options.vxs,:) =   permute(posrect(paramsA(:,4,:)),[3 1 2]);
if ~wantquick
    results.resnorms(options.vxs) = a1.resnorms;
    results.numiters(options.vxs) = a1.numiters;
end

% reshape
results.ang =      reshape(results.ang,      [xyzsize numfits]);
results.ecc =      reshape(results.ecc,      [xyzsize numfits]);
results.expt =     reshape(results.expt,     [xyzsize numfits]);
results.rfsize =   reshape(results.rfsize,   [xyzsize numfits]);
results.R2 =       reshape(results.R2,       [xyzsize numfits]);
results.gain =     reshape(results.gain,     [xyzsize numfits]);
results.hrfshift = reshape(results.hrfshift, [xyzsize numfits]);
results.resnorms = reshape(results.resnorms, [xyzsize 1]);
results.numiters = reshape(results.numiters, [xyzsize 1]);

% add some more stuff
results.meanvol =  meanvol;
results.noisereg = noisereg;
results.params =   paramsA;
results.options = options;