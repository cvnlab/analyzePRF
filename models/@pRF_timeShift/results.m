function results = results(obj, params, metric)
% Packages the model outputs into a results structure
%
% Syntax:
%   results = obj.results(params)
%
% Description:
%   The output of the model is a matrix of parameters fits at each voxel /
%   vertex. This routine performs post-processing of the parameter values
%   and arranges them in a human-readable structure for output.
%
% Inputs:
%   params                - A [v nParams] matrix of parameter values, where
%                           v is the number of voxels/vertices in data.
%   metric                - A [1 v] vector of the metric of the fit at each
%                           voxel or vertex.
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   results               - Structure, with fields for each of the
%                           parameters, the metric, and some metat data.
%


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
results.rfsize =   abs(params(:,3) ./ sqrt(posrect(params(:,5))));
results.gain =     posrect(params(:,4));
results.expt =     posrect(params(:,5));
results.hrfshift = params(:,6);
results.R2 =       metric;

% Add the params themselves
results.params =   params;

% Identify the color scale to be used for plotting the different components
[lb, ub] = obj.bounds;
results.meta.mapField = {'ang','ecc','rfsize','gain','expt','hrfshift','R2'};
results.meta.mapName = {'angle','eccentricity','rfsize','gain','exponent','hrfshift','R2'};
results.meta.mapScale = {'pol','ecc','logJet','linearJet','linearJet','blueRed','grayRed'};
results.meta.mapLabel = {'Polar angle [deg]','Eccentricity [deg]',...
    'Sigma [deg]','response gain [T2* units]',...
    'compressive exponent [au]','shift hrf peak time [secs]','R^2'};
results.meta.mapBounds = {[-180 180],[1 90],[0.1 ub(3)],[lb(4) ub(4)],[lb(5) ub(5)],[lb(6) ub(6)],[0 1]};

end