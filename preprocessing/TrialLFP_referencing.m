
%%
% separate out the grid
% using all frequency plots to identify electrodes of the same grid


data_temp = data_twoDrives{1};
ti = 1:length(Trials); 
outlierElectrodes = [47,59,163];
fprintf('Set the following electrodes to NAN:'); outlierElectrodes

data_temp(:,outlierElectrodes,:) = nan;

gridCoor_1 = 1:61;
gridCoor_2 = 62:90;
gridCoor_3 = 91:151;
gridCoor_4 = 152:211;

ei = 1; % example electrode to look at 
gi = 2;


dataForReferencing(1).data = (data_temp(ti,gridCoor_1,:));
dataForReferencing(1).ePos = P_chamber(gridCoor_1,:); %check if this is the right chamber
dataForReferencing(1).mean = nanmean(dataForReferencing(1).data,2);
dataForReferencing(1).data_ref_mean = dataForReferencing(1).data - dataForReferencing(1).mean;

dataForReferencing(2).data = (data_temp(ti,gridCoor_2,:));
dataForReferencing(2).ePos = P_chamber(gridCoor_2,:); %check if this is the right chamber
dataForReferencing(2).mean = nanmean(dataForReferencing(2).data,2);
dataForReferencing(2).data_ref_mean = dataForReferencing(2).data - dataForReferencing(2).mean;

dataForReferencing(3).data = (data_temp(ti,gridCoor_3,:));
dataForReferencing(3).ePos = P_chamber(gridCoor_3,:); %check if this is the right chamber
dataForReferencing(3).mean = nanmean(dataForReferencing(3).data,2);
dataForReferencing(3).data_ref_mean = dataForReferencing(3).data - dataForReferencing(3).mean;

dataForReferencing(4).data = (data_temp(ti,gridCoor_4,:));
dataForReferencing(4).ePos = P_chamber(gridCoor_4,:); %check if this is the right chamber
dataForReferencing(4).mean = nanmean(dataForReferencing(4).data,2);
dataForReferencing(4).data_ref_mean = dataForReferencing(4).data - dataForReferencing(4).mean;

data_ref_mean(:,gridCoor_1,:) = dataForReferencing(1).data_ref_mean;
data_ref_mean(:,gridCoor_2,:) = dataForReferencing(2).data_ref_mean;
data_ref_mean(:,gridCoor_3,:) = dataForReferencing(3).data_ref_mean;
data_ref_mean(:,gridCoor_4,:) = dataForReferencing(4).data_ref_mean;

set(0,'defaultAxesFontSize',20)
%plotting the average plot
subplot(4,1,1)
imagesc(squeeze(dataForReferencing(1).mean))
colorbar
title('Average LFP for Grid 1 over All Trials')

subplot(4,1,2)
imagesc(squeeze(dataForReferencing(2).mean))
colorbar
title('Grid 2')

subplot(4,1,3)
imagesc(squeeze(dataForReferencing(3).mean))
colorbar
title('Grid 3')

subplot(4,1,4)
imagesc(squeeze(dataForReferencing(4).mean))
xlabel('time (ms)')
ylabel('Trial Number')
colorbar
title('Grid 4')

%look at some example plot

figure
subplot(3,1,1)
plot(squeeze(dataForReferencing(1).data(gi,ei,:)))
title(sprintf('First Trial LFP for Electrode %d of Grid %d',ei,gi))

subplot(3,1,2)
plot(squeeze(dataForReferencing(1).mean(gi,:)))
title(sprintf('Average Grid LFP for all Electrodes in Grid %d',gi))

subplot(3,1,3)
plot(squeeze(dataForReferencing(1).data_ref_mean(gi,ei,:)))
title(sprintf('Grid Mean Subtracted LFP for Electrode %d in Grid %d',ei,gi))
xlabel('Time (ms)')


%% look at the redundancy in the data using PCA


for qi = 1:4
    dataQ = dataForReferencing(qi).data;
    
    clear 'dataQ_PCAresults'
    for ti = 1:size(dataQ,1)
        
        dataQ_trial = squeeze(dataQ(ti,:,:));
        [pca_coeff,pca_score,pca_latent] = pca(dataQ_trial');
        
        
        dataQ_PCAresults.coeff(ti,:,:) = pca_coeff;
        dataQ_PCAresults.score(ti,:,:) = pca_score;
        dataQ_PCAresults.latent(ti,:,:) = pca_latent;
    end
    
    dataQ_PCAresults.scoreMean = squeeze(mean(dataQ_PCAresults.score));
    dataQ_PCAresults.scoreStd = squeeze(std(dataQ_PCAresults.score));
    
    dataQ_PCAresults.latentMean = squeeze(mean(dataQ_PCAresults.latent));
    dataQ_PCAresults.latentStd = squeeze(std(dataQ_PCAresults.latent));
    
    dataQ_PCAresults.latentCum = cumsum(dataQ_PCAresults.latentMean) / sum(dataQ_PCAresults.latentMean);
    
    data_PCAresults(qi) = dataQ_PCAresults;
end

%
dataQ_PCAresults = data_PCAresults(1);
figure
plot(dataQ_PCAresults.latentCum)
hold on
dataQ_PCAresults = data_PCAresults(2);plot(dataQ_PCAresults.latentCum);
dataQ_PCAresults = data_PCAresults(3);plot(dataQ_PCAresults.latentCum);
dataQ_PCAresults = data_PCAresults(4);plot(dataQ_PCAresults.latentCum);
hold off
xlabel('Number of Factors')
ylabel('Cumulative Percent of Variance')
legend('1st Quardrant','2nd Quardrant','3rd Quardrant','4th Quardrant')

%compare across PCs
figure
pci = 1; dataQ_PCAresults = data_PCAresults(1);


plot(dataQ_PCAresults.scoreMean(:,pci))
hold on
pci = 2; plot(dataQ_PCAresults.scoreMean(:,pci))
pci = 3; plot(dataQ_PCAresults.scoreMean(:,pci))
xlabel('time (ms)')
ylabel('Amplitude')
legend('1st PC','2nd PC','3rd PC')
hold off 
%need to do

%look at some principle components

figure
pci = 1;dataQ_PCAresults = data_PCAresults(1);
plot(dataQ_PCAresults.scoreMean(:,pci),'LineWidth',4)
hold on

plot(squeeze(dataQ_PCAresults.score(:,:,pci))',":")

%plot(dataQ_PCAresults.scoreMean(:,pci) + dataQ_PCAresults.scoreStd(:,pci))
%plot(dataQ_PCAresults.scoreMean(:,pci) - dataQ_PCAresults.scoreStd(:,pci),':')
xlabel('time (ms)')
ylabel('Amplitude')
hold off 
%need to do


figure
imagesc(squeeze(dataQ_PCAresults.score(:,:,pci)))

%compare between regions
figure
pci = 1;
dataQ_PCAresults = data_PCAresults(1);plot(dataQ_PCAresults.scoreMean(:,pci))
hold on
dataQ_PCAresults = data_PCAresults(2);plot(dataQ_PCAresults.scoreMean(:,pci))
dataQ_PCAresults = data_PCAresults(3);plot(dataQ_PCAresults.scoreMean(:,pci))
dataQ_PCAresults = data_PCAresults(4);plot(dataQ_PCAresults.scoreMean(:,pci))
xlabel('time (ms)')
ylabel('Amplitude')
legend('1st PC 1st Quardrant','1st PC 2nd Quardrant','1st PC 3rd Quardrant', '1st PC 4th Quardrant')
hold off 
%need to do

%%
% examine the properties of the identified PCs
pci = 1;
qi = 1;

freqParams.tapers = [1 1];
freqParams.Fs = 1000;
freqParams.contFlag = 1;

freqRange = [13 25];

dataQ_PCAresults = data_PCAresults(qi);

dat = dataQ_PCAresults.scoreMean(:,pci)';
[s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, [], [], [], [], [], freqParams.contFlag);
freqRangeIndex = cal_index_freq(freq ,freqRange(1),freqRange(2));

figure
plot(freq(freqRangeIndex ),10*log(s(:,freqRangeIndex)') )

figure
imagesc(10*log(s(:,freqRangeIndex)'))
xlabel('Frequency')
ylabel('Power Amplitude (dB)')
title('Frequency Spectrum over the course of a reach task')
colorbar
axis xy

%%
refData = multi_rereferencing(dataForReferencing(1),0.75);

scatter(dataForReferencing(1).ePos(:,1),dataForReferencing(1).ePos(:,2))

hold on
scatter(refData(1).ePos(:,1),refData(1).ePos(:,2))
hold off

legend('Original Electrode Positions','Neighbor Rereferenced Positions')

%brute force


%%
%look at the data
plot(squeeze(mean(refData.data(:,3,:))))





