function [Palign_b2a, map_a, map_b, map_b_rereg] = splineMapAlign_trans(spl_a, spl_b, x, y)

%%%%note: check x vs y definitions in spline maps!!!


%evaluate splines at specified {x,y} coordinates
map_a = fnval(spl_a, {y, x});

map_b = fnval(spl_b, {y, x});


%get the transformation that optimally co-reigsters the 2 maps
metric = registration.metric.MeanSquares;
optimizer = registration.optimizer.RegularStepGradientDescent;
optimizer.MaximumIterations = 400;
optimizer.MinimumStepLength = 1e-4;
optimizer.MaximumStepLength = 0.05;

[tform] = imregtform(map_b, map_a, 'translation', optimizer, metric);


%get the x/y transformation values
Palign_b2a = tform.T(3,1:2);

%tform is in unit of pixel shifts. Transform to x/y spacing by assuming
%uniform grid
Palign_b2a = Palign_b2a.*[mode(diff(x)) mode(diff(y))];

%also compute transformed map_b
map_b_rereg = imwarp(map_b, tform, 'OutputView', imref2d(size(map_a)));
