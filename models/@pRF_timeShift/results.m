function results = results(obj, params, metric)

% Obj variables
res = obj.res;

% Stimulus center
rCenter = (1+res(1))/2;
cCenter = (1+res(2))/2;

% Map params and metric to a results structure
r = params(:,1);
c = params(:,2);
results.ang = ...
    mod( atan2( rCenter - r, c - cCenter ), 2*pi ) / pi*180;
results.ecc = ...
    sqrt( (rCenter - r).^2 + (c - cCenter).^2);

results.rfsize =   params(:,3);
results.gain =     params(:,4);
results.expt =     params(:,5);
results.hrfshift = params(:,6);
results.R2 =       metric;

end