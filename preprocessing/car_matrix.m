function refMatrix = car_matrix(P_chamber)
% this function returns a list recordings referenced by neiboring
%Inputs

% electrodes

eN = size(P_chamber,1);
refMatrix = eye(eN) - 1/eN * ones(eN);

end