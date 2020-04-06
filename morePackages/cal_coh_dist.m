function P_chamber_coh = cal_coh_dist(P_chamber_temp,cohSeq_firstFrame,distThres)
%cal_coh_dist calculates the average coherence given a distance threshold
%   The coherence between each electrode and electrodes within the distance
%  is averaged
%   Input
%       P_chamber_temp is a nE by 2 matrix of positions
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

nE = size(P_chamber_temp,1);

for ei = 1:nE
    coh_temp = [];
    for et = 1:nE
        interElecD = norm(P_chamber_temp(ei,:) - P_chamber_temp(et,:));
       
        if(interElecD <= distThres)
            if(ei <= et) % only looking at a upper triangular matrix
                coh_temp = [coh_temp; cohSeq_firstFrame(ei,et)];
            else
                coh_temp = [coh_temp; cohSeq_firstFrame(et,ei)];
            end
        end
    end
    
    if(isnan(nanmean(coh_temp)))
        P_chamber_coh(ei,:) = [P_chamber_temp(ei,:) 0];
    else
        P_chamber_coh(ei,:) = [P_chamber_temp(ei,:) nanmean(coh_temp)];
    end
end

end

