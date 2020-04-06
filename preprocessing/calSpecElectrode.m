
function [s_cum_byTrials_raw, s_cum_byTrials,s_cum,s_cum_raw,freq] = calSpecElectrode(dat,s_mean_ETs,s_std_ETs,ei,freqParams)
%This function calculates the spectrum for a given electrode
%Inputs:
%   dat:            trial * ei * time
%   ei:             the electrode to analyze
%   s_mean_ETs:     Average Normalization spectrum for all electrodes 
%                   ei * frequency
%
%   s_std_ETs:      Standard deviation of the specra
%                   ei * frequency
%   freqParams:     parameters for calculting the spectrum
%
%Outputs:
%   s_cum_byTrials_raw: raw spectrum sorted by trials, trial * time * frequency
%   s_cum_byTrials:     z-scored raw data, trial * time * frequency
%   s_cum:              stacked z-scored raw data, time * frequency

s_mean = s_mean_ETs(ei,:)';
s_std = s_std_ETs(ei,:)';

s_cum_raw = [];
s_cum = [];


if exist('freqParams', 'var')
    fn = fieldnames(freqParams);
else
    fn = {};
    freqParams = struct([]);
end
parameters = {'tapers','Fs','contFlag','Fk'}; %rec length in s
defaults   = {[1 2],1000,1,500}; %#ok<NASGU>

%load default
for i=1:length(parameters)
    if ~ismember(parameters{i}, fn)
        eval(['freqParams(1).', parameters{i}, '= defaults{i};'])
    end
end


for ti = 1:size(dat,1)
    datTemp = squeeze((dat(ti,ei,:)))';
    [s, freq] = tfspec(datTemp, freqParams.tapers, freqParams.Fs, [], freqParams.Fk, [], [], [], freqParams.contFlag);
    
    %do normalization using all trials from the same day
    s = s';
    s_cum_byTrials_raw(ti,:,:) = s;
    s_cum_raw = [s_cum_raw s];
    
    %normalize with respect input data
        s_n = (s - repmat(s_mean,1,size(s,2))) ./ repmat(s_std,1,size(s,2));
        s_cum_byTrials(ti,:,:) = s_n;
        s_cum = [s_cum s_n];
end

end