function [refMatrix, refLocation] = bipolar_matrix(P_chamber,distThres)
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


for ei = 1:eN
     
    for et = (ei+1):eN %only looking at unseen electrodes
        if ismember(et,referenced) || ismember(ei,referenced)
            break
        end
        
        interElecD = norm(P_chamber(ei,:) - P_chamber(et,:));
        
        if(interElecD <= distThres)
            % build up the subtraction matrix
            new_row  = zeros(1,length(P_chamber));
            new_row(ei) = 1;
            new_row(et) = -1;
            
            refMatrix = [refMatrix;new_row];
            
            %build location matrix
            new_row  = zeros(1,length(P_chamber));
            new_row(ei) = 1 / 2;
            new_row(et) = 1 /2 ;
            refLocation = [refLocation;new_row];
            
            %skip the rest if finds one
            %save to referenced matrix
            referenced = [referenced;et;ei];
           
        end
    end
end

end