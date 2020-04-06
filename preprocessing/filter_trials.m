function [trialInfo_new,Trials] = filter_trials(trialInfo,Trials)
% oh this function just filters the trials based on predefined criterion
% currently supports
%   random sequence:  set trialInfo.sequence_random = 1
%   fitlter long acq time: trialInfo.filterAcq = 1 and filterAcq_time
% example usage:
%   trialInfo.sequence_random = 1;
%   trialInfo.filterAcq = 1;
%   trialInfo.filterAcq_time = 5; %s
%   [trialInfo,trials] = filter_trials(trialInfo,trials);

    if trialInfo.sequence_random
        v = [Trials.SequenceRandom];
        Trials = Trials(v ==1);
    end
    if trialInfo.filterAcq
        stopToAcq = [Trials.ReachStop] - [Trials.TargAq];
        v = abs(stopToAcq) <= trialInfo.filterAcq_time;
        Trials = Trials(v);
    end
    trialInfo.num_trials = length(Trials);
    
    trialInfo_new = trialInfo;



end

