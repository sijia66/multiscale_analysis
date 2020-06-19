function targ_predict_accuracy = lda_sample_time(trLfpData,timeVec,trial_TargetAssigned, win_sizes,start_times)
%this function extracts features averaging a sliding window 
%   Detailed explanation goes here


targ_predict_accuracy = nan(length(win_sizes),length(start_times));


for wt = 1:length(win_sizes)
    for st_i =  1:length(start_times)
        
        timeIdx = (timeVec >= start_times(st_i)) & (timeVec <=  start_times(st_i) + win_sizes(wt));
        if(start_times(st_i) + win_sizes(wt) > max(timeVec))
            continue
        end
        %neural_feature: num of trials by number of feature or recordings.
        neural_feature = squeeze(mean(trLfpData(:,:,timeIdx),3));
        target_hat = runLeaveOneOutClassification(neural_feature,trial_TargetAssigned);
        targ_predict_accuracy(wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
       
        
    end
    %disp(['Finished:',num2str(wt/length(win_sizes))])
end


end

