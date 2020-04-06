function [s_z_trials,s_z_info] = caltf_trial(data,freqParams)
% this function calculates time frequency decomposition for an electrode of
% a given trial work along the second dimension
% Inputs:
%   data: number of trials by time points
%   dataTrials: number of trials
%   freqParams: parameters to do frequency decomposition
%Packages used
%   tfspec from amy's post doc lab
%   cal_index_freq
%Output
%   s_z normalize spectrum

s_z_trials = [];
for ti = 1:size(data,1)
        
    % now look at the trials
    % we are gonna get
    dataTemp = squeeze(data(ti,:));
    %(timeReachStart - timeInterest  + 1) : (timeReachStart + timeInterest) ));
    
    [s, freq] = tfspec(dataTemp, freqParams.tapers, freqParams.Fs, freqParams.dn, [], [], [], [], freqParams.contFlag);
   
    
    s_z_trials(ti,:,:) = s;
    
end

s_length = size(s,1);
s_timePoints = freqParams.tapers(1)* 0.5 +(1:s_length)*freqParams.tapers(1) * 0.1;

s_z_info.xlabel = s_timePoints;
s_z_info.ylabel = freq;
end