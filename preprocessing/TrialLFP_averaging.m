% find the more varied electrodes

ti = 20;


data_trial = squeeze(dataTarget.target1(:,ti,:));

%averaging over electrode to estimate the common signal across the
%across the entire brain area
data_trial_norm = zeros(size(data_trial));
%normalize 
for ci = 1:size(data_trial,1)
    data_trial_norm(ci,:) = data_trial(ci,:) / norm(data_trial(ci,:));
end

%average across all electrodes
data_trial_mean = squeeze(mean(data_trial_norm));


figure

hold on
plot(data_trial_norm',":")
plot((data_trial_mean), 'LineWidth',2)
hold off
title(strcat('Electrode',num2str(ti)))
xlabel('time (ms)')

%%
% use standard deviation to find suspicious electrodes
ti = 1;
loweredE_SC32 = 1; 
%trialInfo.badECOG = [47,59,163];

timeRange = 1801:2200;
timeRange_display = timeRange - 2000;

%averaging over electrode to estimate the common signal across the
%across the entire brain area

dataTemp = dataTarget_SC32.target1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if loweredE_SC32
    Ndrives = size(Trials(ti).Depth,2);
    depthProfile = (Trials(1).Depth{1,2})'; % in micron
    
    ei = 1:length(depthProfile);
    ei(depthProfile ==0) = [];
    data_trial = squeeze(dataTemp(ti,ei,:));
    
    disp('Only selected lowered electrodes')
else
    data_trial = squeeze(dataTemp(ti,:,:));
end

data_trial_norm = zeros(size(data_trial));
%normalize 
for ci = 1:size(data_trial,1)
    data_trial_norm(ci,:) = data_trial(ci,:) / norm(data_trial(ci,:));
end


data_electrode_mean = mean(data_trial_norm);
dataTemp = data_trial_norm(:,timeRange);

data_electrode_std = sum((dataTemp - ...
         repmat(data_electrode_mean(timeRange),size(dataTemp,1),1)) .^2,2) ...
            / size(dataTemp,1);
        
figure
subplot(2,1,1)
hold on
plot(timeRange_display, data_trial_norm(:,timeRange)',":")
plot(timeRange_display,data_electrode_mean(timeRange), 'LineWidth',2)
hold off
title('SC32 - LFP')
% xlabel('Time from  Movement Onset (ms)')
ylim([-0.06 0.06])
%ylabel('Normalized Amp.')

%ECOG
dataTemp = dataTarget.target1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if loweredE_SC32
    Ndrives = size(Trials(ti).Depth,2);
    depthProfile = (Trials(1).Depth{1,2})'; % in micron
    
    ei = 1:length(depthProfile);
    ei(depthProfile ==0) = [];
    data_trial = squeeze(dataTemp(ti,ei,:));
    
    disp('Only selected lowered electrodes')
else
    data_trial = squeeze(dataTemp(ti,:,:));
end

data_trial_norm = zeros(size(data_trial));
%normalize 
for ci = 1:size(data_trial,1)
    data_trial_norm(ci,:) = data_trial(ci,:) / norm(data_trial(ci,:));
end


data_electrode_mean = mean(data_trial_norm);
dataTemp = data_trial_norm(:,timeRange);

data_electrode_std = sum((dataTemp - ...
         repmat(data_electrode_mean(timeRange),size(dataTemp,1),1)) .^2,2) ...
            / size(dataTemp,1);

subplot(2,1,2)
hold on
plot(timeRange_display, data_trial_norm(:,timeRange)',":")
plot(timeRange_display,data_electrode_mean(timeRange), 'LineWidth',2)
hold off
title('ECOG')
xlabel('Time from  Movement Onset (ms)')
%ylabel('Normalized Amp.')
ylim([-0.06 0.06])

%%

[pca_coeff,pca_score,pca_latent] = pca(data_trial_norm');
pca_latentCum = cumsum(pca_latent) / sum(pca_latent);

% figure
% plot(data_electrode_std)
% title('Standard Deviation for Each Electrode')
% xlabel('Electrode Number')
% 
% 
% figure
% subplot(1,2,1)
% plot(pca_latentCum)
% title('Cumulative Variance')
% 
% subplot(1,2,2)
% plot(pca_score(:,1))
% title('First PC')


%%
freqParams.tapers = [1 1];
freqParams.Fs = 1000;
freqParam.dn = 0.1;
freqParams.contFlag = 1;

freqRange = [0 50];


dat = data_electrode_mean;
[s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, freqParam.dn, [], [], [], [], freqParams.contFlag);
freqRangeIndex = cal_index_freq(freq ,freqRange(1),freqRange(2));

figure
plot(freq(freqRangeIndex ),10*log(s(:,freqRangeIndex)') )
xlabel('Frequency')
ylabel('Power Amplitude (dB)')
title('Frequency Spectrum over the course of a reach task')

figure
imagesc(10*log(s(:,freqRangeIndex)'))
colorbar
axis xy


%%
% look at average response ti

plot()







