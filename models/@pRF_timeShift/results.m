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
pixelsPerDegree = obj.pixelsPerDegree;
screenMagnification = obj.screenMagnification;

% Stimulus center
rCenter = (1+res(1))/2;
cCenter = (1+res(2))/2;

% Rows and columns
r = params(:,1);
c = params(:,2);

% Map params and metric to a results structure
results.cartX = c - cCenter;
results.cartY = rCenter - r;
results.angle = ...
    mod( atan2( rCenter - r, c - cCenter ), 2*pi ) / pi*180;
results.eccen = ...
    sqrt( (rCenter - r).^2 + (c - cCenter).^2);
results.sigma =   abs(params(:,3) ./ sqrt(params(:,5)));
results.gain =     params(:,4);
results.expt =     params(:,5);
results.hrfshift = params(:,6);
results.R2 =       metric;

% Convert the result data to units of visual angle in degrees. In this
% step, account as well for the effect of corrective lenses. The stimulus
% screen is subject to magnification / minification if the subject is
% wearing corrective lenses. We account for this effect here. While in
% principle this could be rolled into the pixelsPerDegree variable, we
% prefer to keep these separate to aid clear book-keeping.
fieldsToAdjust = {'cartX','cartY','eccen','sigma'};
for ii = 1:length(fieldsToAdjust)
    results.(fieldsToAdjust{ii}) = ...
        ( results.(fieldsToAdjust{ii}) ./pixelsPerDegree) .* screenMagnification;
end

% Add the params themselves
results.params =   params;

% Identify the color scale to be used for plotting the different components
[lb, ub] = obj.bounds;
results.meta.mapField = {'cartX','cartY','angle','eccen','sigma','gain','expt','hrfshift','R2'};
results.meta.mapScale = {'blueRed','blueRed','angle','eccen','logJet','linearJet','linearJet','blueRed','grayRed'};
results.meta.mapLabel = {'Horizontal [deg]','Vertical [deg]','Polar angle [deg]','Eccentricity [deg]',...
    'Sigma [deg]','response gain [T2* units]',...
    'compressive exponent [au]','shift hrf peak time [secs]','R^2'};
results.meta.mapBounds = {...
    [lb(1) ub(1)]./pixelsPerDegree,...
    [lb(2) ub(2)]./pixelsPerDegree,...
    [-180 180],...
    [1 90],...
    [0.5 20],...
    [lb(4) ub(4)],...
    [lb(5) ub(5)],...
    [lb(6) ub(6)],...
    [0 1]};

end