%% 0  Notes on required input variables and functions 
%variables

%experiment: to get the drive coordinates
%cohDat: calculated from spontCoh
%specDat: from spontCoh
%freqRange: to specify the frequency to analyze


%% 1.2 Generate a data structure for other scripts
lfpType = 'lfp';
Fs_lfp = 1000;
dataTemp = data_ref_mean;

dataCoh(1).tsDat_tr = permute(dataTemp,[2 3 1]);
dataCoh(1).driveName = driveNameAnalyze;
dataCoh(1).datType = lfpType;
dataCoh(1).Fs = Fs_lfp;
dataCoh(1).driveNum = 1; %for ECOG

%% 1.2.2 then load SC32
dataTemp = data;
dataCoh(2).tsDat_tr = permute(dataTemp,[2 3 1]);
dataCoh(2).driveName = driveNameAnalyze; 
dataCoh(2).datType = lfpType;
dataCoh(2).Fs = Fs_lfp;
dataCoh(2).driveNum = 2; %for SC32

%% 
fileName = '180328_specCoh_ECOG_gridMeanSsubtracted_SC32_100trials_Sept4';

tic   
%we only need to look   at lower frequencies  

cohParams.tapers = [0.1 10] ;  
datParams.fileType = 'lfp';
[cohDat, specDat, data] = spontCoh('', [], '', cohParams, datParams, dataCoh);

toc  
save(fileName,'cohDat','specDat','data')
    


%% 2 plot as an image 

% load the experiment from the recording
P_chamber_SC32 = ...
    drmap_coordTransform_microdrive2chamber(experiment, 'LM1_SC32_1');
P_chamber_SC32 = P_chamber_SC32';  

P_chamber_ECOG_3 = ...
    drmap_coordTransform_microdrive2chamber(experiment, 'LM1_ECOG_3');
P_chamber_ECOG_3 = P_chamber_ECOG_3';


cohSeq_ECOG = cohDat(1).coh; %1 for ECOG32
cohSeq_SC32 = cohDat(2).coh; %2 for SC32
cohSeq_ECOG_SC32 = cohDat(3).coh; %3 interDrive spectral coherency 

% %need to check frequency range
% freqRangeFull = [13 30];
% 
% freqRangeIndexFull  = cal_index_freq(specDat(1).freq ,freqRangeFull(1),freqRangeFull(2));
% 
% freqRangeBetaFeature = [20 24];
% freqRangeIndexBetaFeature = cal_index_freq(specDat(1).freq ,freqRangeBetaFeature(1),freqRangeBetaFeature(2));
% freqRange = freqRangeFull;
% freqRangeIndex = setdiff(freqRangeIndexFull,freqRangeIndexBetaFeature);

freqRange = [0 4];
freqRangeIndex  = cal_index_freq(specDat(1).freq ,freqRange(1),freqRange(2));


cohSeq_fr_ECOG = mean(cohSeq_ECOG(:,:,freqRangeIndex),3);
cohSeq_fr_SC32 = mean(cohSeq_SC32(:,:,freqRangeIndex),3);
cohSeq_fr_ECOG_SC32 = mean(cohSeq_ECOG_SC32(:,:,freqRangeIndex),3);

%%
%look at the distributions

perIndex = 99;
fprintf('Percentile used to construct the map: %d\n',perIndex)

cohSeq_fr_temp = cohSeq_fr_ECOG;

[nRow,nCol] = size(cohSeq_fr_temp);
cohSeq_fr_temp = reshape(cohSeq_fr_temp,nRow*nCol,1);

coh_th(1) = prctile(cohSeq_fr_temp, perIndex );
fprintf('ECOG: %.4f\n',coh_th(1));


figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,3,1)
hist(cohSeq_fr_temp)

hold on
line([coh_th(1), coh_th(1)], ylim, 'LineWidth', 2, 'Color', 'r');
hold off
xlabel('ECoG Coherency')
ylabel('Counts')

cohSeq_fr_temp = cohSeq_fr_SC32;
[nRow,nCol] = size(cohSeq_fr_temp);
cohSeq_fr_temp = reshape(cohSeq_fr_temp,nRow*nCol,1);

coh_th(2) = prctile(cohSeq_fr_temp, perIndex );
fprintf('SC32: %.4f\n',coh_th(2));

subplot(1,3,2)
hist(cohSeq_fr_temp)
hold on
line([coh_th(2), coh_th(2)], ylim, 'LineWidth', 2, 'Color', 'r');
hold off
xlabel('SC32 Coherency')
ylabel('Counts')

cohSeq_fr_temp = cohSeq_fr_ECOG_SC32;
[nRow,nCol] = size(cohSeq_fr_temp);
cohSeq_fr_temp = reshape(cohSeq_fr_temp,nRow*nCol,1);

coh_th(3) = prctile(cohSeq_fr_temp, perIndex );
fprintf('ECOG vs. SC32: %.4f\n',coh_th(3));

subplot(1,3,3)
hist(cohSeq_fr_temp)
hold on
line([coh_th(3), coh_th(3)], ylim, 'LineWidth', 2, 'Color', 'r');
hold off
xlabel('ECOG vs. SC32 Coherency')
ylabel('Counts')




%%
%adjust the thresholds
%skip this if using the percentiles
coh_th(1) = 0.9 %ECOG


%%
% generate the maps and the network

figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,3,1)
imagesc(cohSeq_fr_ECOG)
title(sprintf('ECOG Averaged Between %0.2f and %0.2f Hz',freqRange(1),freqRange(2)))
ylabel('ECOG')
xlabel('ECOG')
colorbar

subplot(2,3,2)
imagesc(cohSeq_fr_SC32)
title(sprintf('SC32 Averaged Between %0.2f and %0.2f Hz',freqRange(1),freqRange(2)))
ylabel('SC32')
xlabel('SC32')
colorbar

subplot(2,3,3)
imagesc(cohSeq_fr_ECOG_SC32)
title(sprintf(' Inter-drive Coherency Averaged Between %0.2f and %0.2f Hz',freqRange(1),freqRange(2)))
ylabel('SC32')
xlabel('ECOG')
colorbar


%%plotting the corresponding maps

graphParam.threshold = coh_th(1);
graphParam.thresholdPercentile = perIndex;
graphParam.driveName = 'ECOG';
graphParam.drive_idx =  1; %1 for ECOG; 2 for SC32
cMatrix = cohSeq_fr_ECOG;

graphParam.textShift = 0.1;
P_chamber = P_chamber_ECOG_3;

graphParam.freqRange = freqRange;
graphParam.labelElectrode = true;
graphParam.textShift = 0.1;
graphParam.acList = [experiment.hardware.microdrive(graphParam.drive_idx).electrodes.acquisitionid];

subplot(2,3,4)
mapCoherenceGraph(cMatrix,P_chamber, graphParam)

%%%%%%%%%%%%%%% settings need to be changed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
graphParam.threshold = coh_th(2);
graphParam.thresholdPercentile = perIndex;
graphParam.driveName = 'SC32';
graphParam.drive_idx =  2; %1 for ECOG; 2 for SC32
cMatrix = cohSeq_fr_SC32;
P_chamber = P_chamber_SC32;

graphParam.freqRange = freqRange;
graphParam.labelElectrode = true;
graphParam.textShift = 0.1;
graphParam.acList = [experiment.hardware.microdrive(graphParam.drive_idx).electrodes.acquisitionid];

subplot(2,3,5)
mapCoherenceGraph(cMatrix,P_chamber, graphParam)


%%%%%%%%%%%%%%% settings need to be changed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cMatrix = cohSeq_fr_ECOG_SC32;
graphParam.threshold = coh_th(3);
graphParam.thresholdPercentile = perIndex;

graphParam.driveName1 = 'SC32';
graphParam.P_chamber1 = P_chamber_SC32;

graphParam.driveName2 = 'ECOG';
graphParam.P_chamber2 = P_chamber_ECOG_3;

graphParam.freqRange = freqRange;
graphParam.labelElectrode = true;
graphParam.drive_idx =  2; %1 for ECOG
graphParam.textShift = 0.1;

subplot(2,3,6)

A = cMatrix > graphParam.threshold;
inter_gplot(A,graphParam)



%%
% compare between days and apply the same threshold


freqRange = [13 30];
fileNames = {'/180328/180328_specCoh_ECOG_gridMeanSsubtracted_SC32_101trials_June14.mat',...
             '/180329/results/180329_specCoh_ECOG_gridMeanSsubtracted_SC32_121trials_June16.mat',...
             '/180408/180408_specCoh_ECOG_gridMeanSsubtracted_SC32_100trials_June19'};

perIndex = 99;


for fi = 1:length(fileNames)
    
    %load the files 
    dir = strcat(MONKEYDIR,fileNames{fi});
    load(dir)

    cohSeq_ECOG = cohDat(1).coh; %1 for ECOG32
    cohSeq_SC32 = cohDat(2).coh; %2 for SC32
    cohSeq_ECOG_SC32 = cohDat(3).coh; %3 interDrive spectral coherency 


    freqRangeIndex  = cal_index_freq(specDat(1).freq ,freqRange(1),freqRange(2));
    cohSeq_fr_ECOG = mean(cohSeq_ECOG(:,:,freqRangeIndex),3);
    cohSeq_fr_SC32 = mean(cohSeq_SC32(:,:,freqRangeIndex),3);
    cohSeq_fr_ECOG_SC32 = mean(cohSeq_ECOG_SC32(:,:,freqRangeIndex),3);

    cohSeq_fr_ECOG_pool(:,:,fi) = cohSeq_fr_ECOG;

end


% find the thresholds for entire data set
fprintf('Percentile used to construct the map: %d\n',perIndex)
cohSeq_fr_temp =  cohSeq_fr_ECOG_pool;
[nRow,nCol,nFi] = size(cohSeq_fr_temp);
cohSeq_fr_temp = reshape(cohSeq_fr_temp,nRow*nCol * nFi,1);
coh_th(1) = prctile(cohSeq_fr_temp, perIndex );
fprintf('ECOG: %.4f\n',coh_th(1));
%
hist(cohSeq_fr_temp)
hold on
line([coh_th(1), coh_th(1)], ylim, 'LineWidth', 2, 'Color', 'r');
hold off
xlabel('ECoG Coherency')
ylabel('Counts')



for fi = 1:length(fileNames)

    % calculate connectivty degree
    cMatrix = applyNetworkThreshold(squeeze(cohSeq_fr_ECOG_pool(:,:,fi)),coh_th(1));
    
    %calculate the network degree
    connectDegree(:,fi) = calNetworkDegree(cMatrix,'sym');
    clustercoefficients(:,fi) = calNetworkCluster(cMatrix);
    temppathlengths(:,fi) = calNetworkPathLength(cMatrix);
end




%%
%visulization
subplot(311)
plot(connectDegree)
ylabel('Number of Connections')
title('Connectivity Degree')
legend('180328','180329','180408')

subplot(312)
plot(clustercoefficients)
ylabel('Fraction of Connections')
title('Clustering Coefficient')

subplot(313)
plot(temppathlengths)
xlabel('Electrodes')
ylabel('Number of Nodes')
title('Path length')


%% plot the results in matrix form
% 


figure
subplot(331)
[map,mapInfo]= scatterToMatrix(P_chamber_ECOG_3,connectDegree(:,1));
imagesc(map)
title('180328')
scatterToMatrix_plot_label

subplot(332)
[map,mapInfo]= scatterToMatrix(P_chamber_ECOG_3,connectDegree(:,2));
imagesc(map)
title('180329')
scatterToMatrix_plot_label

subplot(333)
[map,mapInfo]= scatterToMatrix(P_chamber_ECOG_3,connectDegree(:,3));
imagesc(map)
title('180408')
scatterToMatrix_plot_label


subplot(334)
[map,mapInfo]= scatterToMatrix(P_chamber_ECOG_3,clustercoefficients(:,1));
imagesc(map)
title('180328')
scatterToMatrix_plot_label

subplot(335)
[map,mapInfo]= scatterToMatrix(P_chamber_ECOG_3,clustercoefficients(:,2));
imagesc(map)
title('180329')
scatterToMatrix_plot_label

subplot(336)
[map,mapInfo]= scatterToMatrix(P_chamber_ECOG_3,clustercoefficients(:,3));
imagesc(map)
title('180408')
scatterToMatrix_plot_label

subplot(337)
[map,mapInfo]= scatterToMatrix(P_chamber_ECOG_3,temppathlengths(:,1));
imagesc(map)
title('180328')
scatterToMatrix_plot_label

subplot(338)
[map,mapInfo]= scatterToMatrix(P_chamber_ECOG_3,temppathlengths(:,2));
imagesc(map)
title('180329')
scatterToMatrix_plot_label

subplot(339)
[map,mapInfo]= scatterToMatrix(P_chamber_ECOG_3,temppathlengths(:,3));
imagesc(map)
title('180408')
scatterToMatrix_plot_label
%

%% extract common electrodes 
perCentile = 95;

perCols_1 = calTopPercentile(connectDegree(:,1), perCentile);
perCols_2 = calTopPercentile(connectDegree(:,2), perCentile);
perCols_3 = calTopPercentile(connectDegree(:,3), perCentile);

perCols_intersect = intersect(perCols_1,perCols_2);
perCols_intersect = intersect(perCols_intersect,perCols_1);
