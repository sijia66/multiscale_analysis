function resultData = multi_rereferencing(dataForReferencing,distThres)
% this function returns a list recordings referenced by neiboring
%Inputs
%   dataForReferencing has two fields
%       data : number of trials * N electrodes * time points
%       ePos: N electrodes * 2 (x,y coordinates)
% electrodes

data = dataForReferencing.data;
P_chamber = dataForReferencing.ePos;

eN = size(P_chamber,1);

%this is going to be a search algorithm



eNew = 1;

for ei = 1:eN
     
    for et = (ei+1):eN %only looking at unseen electrodes
        interElecD = norm(P_chamber(ei,:) - P_chamber(et,:));
        
        if(interElecD <= distThres)
            % subtract and record the location distance
            dataTemp = data(:,ei,:) - data(:,et,:);
            posTemp = (P_chamber(ei,:) + P_chamber(et,:)) / 2;
            
            refData(:,eNew,:) = dataTemp;
            
            reflocations(eNew,:) = posTemp;
            
            eNew = eNew + 1;
        end
    end
end

resultData.data = refData;
resultData.ePos = reflocations;

end