function P_chamber_coh = cal_coh_dist_inter(P_chamber_temp1,P_chamber_temp2,cohSeq_firstFrame,distThres)
%cal_coh_dist calculates the average coherence given a distance threshold
%   The coherence between two drives within the distance
%  is averaged
%   Input
%       P_chamber_temp1 is a nE by 2 matrix of positions for first drive
%       P_chamber_temp2 is a nE by 2 matrix of positions for second drive
%       cohSeq_firstFrame is a nE by nE coherence matrix
%       distThres is the distance threshold
%   Output
%       P_chamber_coh is a nE by 3 matrix of the positions and the average
%       coherence
%   
%   Under construction:
%       need to do error checking

%P_chamber_temp = P_chamber_SC32

%distThres = 1.5; % search thresold, in mm

nE1 = size(P_chamber_temp1,1);
nE2 = size(P_chamber_temp2,1);


for ei = 1:nE1
    coh_temp = [];
    for et = 1:nE2
        interElecD = norm(P_chamber_temp1(ei,:) - P_chamber_temp2(et,:));
        
        if(interElecD <= distThres)
            coh_temp = [coh_temp; cohSeq_firstFrame(ei,et)];
        end
    end
    
    if(isnan(nanmean(coh_temp)))
        P_chamber_coh(ei,:) = [P_chamber_temp1(ei,:) 0];
    else
        P_chamber_coh(ei,:) = [P_chamber_temp1(ei,:) nanmean(coh_temp)];
    end
end

end

