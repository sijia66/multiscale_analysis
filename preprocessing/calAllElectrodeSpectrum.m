function [s_mean_ETs,s_std_ETs,freq] = calAllElectrodeSpectrum(data,freqParams)
%   Detailed explanation goes here

% trLen = 1;% in s, this determines our resolution
% Fs = 1000; 
% contflag = 1;
% taperFreq = 1;

if exist('freqParams', 'var')
    fn = fieldnames(freqParams);
else
    fn = {};
    freqParams = struct([]);
end
parameters = {'tapers','Fs','contFlag'}; %rec length in s
defaults   = {[1 4],1000,1}; %#ok<NASGU>

for i=1:length(parameters)
    if ~ismember(parameters{i}, fn)
        eval(['freqParams(1).', parameters{i}, '= defaults{i};'])
    end
end


trialN = size(data,1);
eN = size(data,2);
pN = size(data,3); % number of electrodes

disp('Start crunching the numbers')
disp('Go grab a coffee')

for ei = 1:eN
   
    %pool together all results and get the spectrum
    datTemp = squeeze(data(:,ei,:));
    dat = reshape(datTemp',trialN * pN,1)';
    [s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, [], [], [], [], [], freqParams.contFlag);
    
    s_mean = mean(s)';
    s_std = std(s)';
    
    s_mean_ETs(ei,:) = s_mean';
    s_std_ETs(ei,:) = s_std';
    
    if (mod(ei,10) == 0)
        disp(strcat('Finished:',num2str(ei)))
    end
end

disp('Done calculating')
end

