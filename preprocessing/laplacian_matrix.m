function [refMatrix, refLocation,electrode_index] = laplacian_matrix(P_chamber,distThres)
% this function returns a list recordings referenced by neiboring
%Inputs
%   dataForReferencing has two fields
%       ePos: N electrodes * 2 (x,y coordinates)
% electrodes


eN = size(P_chamber,1);

%this is going to be a search algorithm
refMatrix = [];
refLocation = []; 
referenced = [];

unused_electrodes = [];
for ei = 1:eN
    
        
    interElecD = vecnorm(repmat(P_chamber(ei,:),length(P_chamber),1)...
        - P_chamber, 2, 2);
    
    interElecD_ind = (interElecD <= distThres)...
        & interElecD > 0; %return a logical index
    
    %if no adjacent electrodes, mark it and throw it later
    if all(interElecD_ind == 0 )
        unused_electrodes = [unused_electrodes;ei];
    end
    %subtract from surroundings
    refMatrix(ei,:) = - double(interElecD_ind) / sum(interElecD_ind);
    refMatrix(ei,ei) = 1;    
end

%delete unused electrodes and P_chamber info 
electrode_index = 1:eN;
electrode_index(unused_electrodes) = [];
refMatrix(unused_electrodes,:) = [];
refLocation = P_chamber;
refLocation(unused_electrodes,:) = [];

end