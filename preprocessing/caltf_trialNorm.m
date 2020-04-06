function [s_z_trials,s_z_info] = caltf_trialNorm(data,dataTrials,trialInfo,freqParams)
% this function calculates time frequency decomposition for an electrode of
% a given trial
% Inputs:
%   data: number of trials by time points
%   dataTrials: number of trials
%   trialInfo.timeReachStart: alignment to reachStart
%   freqParams: parameters to do frequency decomposition
%Packages used
%   tfspec from amy's post doc lab
%   cal_index_freq
%Output
%   s_z normalize spectrum

s_z_trials = [];
for ti = 1:size(data,1)
    
    dataTemp = squeeze(data(ti,:));
    %first normalize
    dataTemp = dataTemp / norm(dataTemp);
    
    timeAqGo = dataTrials(ti).StartAq - dataTrials(ti).ReachStart + trialInfo.timeReachStart;
    timeGo = dataTrials(ti).Go - dataTrials(ti).ReachStart + trialInfo.timeReachStart; %relative to reach start
    timeTargsOn = dataTrials(ti).TargsOn - dataTrials(ti).ReachStart + trialInfo.timeReachStart;
    dataTemp = dataTemp(timeAqGo:timeTargsOn);
    
    %do frequency calculation
    [s, freq] = tfspec(dataTemp, freqParams.tapers, freqParams.Fs, freqParams.dn, [], [], [], [], freqParams.contFlag);
    
    s_n_mean = squeeze(mean(s));
    s_n_std = squeeze(std(s));
    
    % now look at the trials
    % we are gonna get
    dataTemp = squeeze(data(ti,:));
    %(timeReachStart - timeInterest  + 1) : (timeReachStart + timeInterest) ));
    %first normalize
    dataTemp = dataTemp / norm(dataTemp);
    
    [s, freq] = tfspec(dataTemp, freqParams.tapers, freqParams.Fs, freqParams.dn, [], [], [], [], freqParams.contFlag);
    
    %get z-score
    sN = size(s,1);
    s_z = (s - repmat(s_n_mean,sN,1)) ./ repmat(s_n_std,sN,1);
    
    s_z_trials(ti,:,:) = s_z;
    
end

s_length = size(s_z,1);
s_timePoints = freqParams.tapers(1)* 0.5 +(1:s_length)*freqParams.tapers(1) * 0.1;

s_z_info.xlabel = s_timePoints;
s_z_info.ylabel = freq;
end