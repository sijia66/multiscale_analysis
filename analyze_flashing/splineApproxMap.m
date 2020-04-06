function [splApprox, Map_hat] = splineApproxMap(X, Y, Map, options)

if ~exist('options', 'var') || isempty(options)
    
    options = struct([]);
end

parameters = {'order', 'numKnts', 'kntsX', 'kntsY'}; %rec length in s
defaults   = {4,         6,         [],      []}; %#ok<NASGU>

fn = fieldnames(options);
for i=1:length(parameters)
    if ~ismember(parameters{i}, fn)
        eval(['options.', parameters{i}, '= defaults{i};'])
    end
end


%if knots not specified, use uniform spacing:
if isempty(options.kntsX)
    unifX = linspace(X(1), X(end), options.numKnts);
    options.kntsX = augknt(unifX, options.order);
end

if isempty(options.kntsY)
    unifY = linspace(Y(1), Y(end), options.numKnts);
    options.kntsY = augknt(unifY, options.order);
end


%fit the spline
splApprox = spap2({options.kntsX, options.kntsY}, {options.order, options.order}, {Y, X}, Map);


%also generate map estimated by the spline fine
Map_hat = fnval(splApprox, {Y, X});