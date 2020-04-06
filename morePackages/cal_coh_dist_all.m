function P_chamber_coh = cal_coh_dist_all(P_chamber_temp,cohSeq_firstFrame)
%cal_coh_dist calculates the average coherence given a distance threshold
%   The coherence between each electrode and electrodes within the distance
%  is averaged
%   Input
%       P_chamber_temp is a nE by 2 matrix of positions
%       cohSeq_firstFrame is a nE by nE coherence matrix
%   Output
%       P_chamber_coh is a nE*(nE+1) / 2 by 2 matrix of the positions and  average
%       coherence
%   
%   Under construction:
%       need to do error checking

%P_chamber_temp = P_chamber_SC32

%distThres = 1.5; % search thresold, in mm

nE = size(P_chamber_temp,1);

P_chamber_coh = [];

for ei = 1:nE
    for et = (ei+1):nE
        
        interElecD = norm(P_chamber_temp(ei,:) - P_chamber_temp(et,:));
        coh_temp = [interElecD cohSeq_firstFrame(ei,et)];
        P_chamber_coh = [ P_chamber_coh; coh_temp];
        
    end
end

%     if(isnan(nanmean(coh_temp)))
%         P_chamber_coh(ei,:) = [P_chamber_temp(ei,:) 0];
%     else
%         P_chamber_coh(ei,:) = [P_chamber_temp(ei,:) nanmean(coh_temp)];
%     end

