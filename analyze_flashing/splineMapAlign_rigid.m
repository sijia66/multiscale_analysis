function [Palign_b2a, thetaAlign_b2a, map_a, map_b, map_b_rereg, tform] = splineMapAlign_rigid(spl_a, spl_b, x, y)



%evaluate splines at specified {x,y} coordinates
map_a = fnval(spl_a, {y, x});

map_b = fnval(spl_b, {y, x});


%get the transformation that optimally co-reigsters the 2 maps
metric = registration.metric.MeanSquares;
optimizer = registration.optimizer.RegularStepGradientDescent;
optimizer.MaximumIterations = 400;
optimizer.MinimumStepLength = 1e-4;
optimizer.MaximumStepLength = 0.04;

[tform] = imregtform(map_b, map_a, 'rigid', optimizer, metric);


%get the x/y transformation values
Palign_b2a = tform.T(3,1:2);
thetaAlign_b2a = atan2d(tform.T(1,2),tform.T(1,1));

%tform is in unit of pixel shifts. Transform to x/y spacing by assuming
%uniform grid
Palign_b2a = Palign_b2a.*[mode(diff(x)) mode(diff(y))];

%also compute transformed map_b
map_b_rereg = imwarp(map_b, tform, 'OutputView', imref2d(size(map_a)));
