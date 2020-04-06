function [s_z_trials,s_z_info] = caltf_trialNorm_long(data,s_n_mean,s_n_std,freqParams)
% this function calculates time frequency decomposition for an electrode of
% a given trial normalized by precalculated mean and std for entire
% recording
% Inputs:
%   data: number of trials by time points
%   s_n_mean: 1 by frequencies
%   s_n_std: 1 by frequencies

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
    %first normalize

    
    [s, freq] = tfspec(dataTemp, freqParams.tapers, freqParams.Fs, freqParams.dn, freqParams.freqRange, [], [], [], freqParams.contFlag);
    
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