%% 1.2 we copied and pasted code from the main code mapping coherence graph
lfpType = 'lfp';
Fs_lfp = 1000;
dataTemp = data_ref_mean;

dataCoh(1).tsDat_tr = permute(dataTemp,[2 3 1]);
dataCoh(1).driveName = 'LM1_ECOG_3';
dataCoh(1).datType = lfpType;
dataCoh(1).Fs = Fs_lfp;
dataCoh(1).driveNum = 1; %for ECOG

%% 1.2.2 then load SC32
dataTemp = data;
dataCoh(2).tsDat_tr = permute(dataTemp,[2 3 1]);
dataCoh(2).driveName = 'LM1_SC32_1'; 
dataCoh(2).datType = lfpType;
dataCoh(2).Fs = Fs_lfp;
dataCoh(2).driveNum = 2; %for SC32


%% then we load

fileName = '180328_specCoh_ECOG_gridMeanSsubtracted_SC32_100trials_Sept9';

tic   
%we only need to look   at lower frequencies  

cohParams.tapers = [0.1 10] ;  
datParams.fileType = 'lfp';
[cohDat, specDat, data] = spontCoh('', [], '', cohParams, datParams, dataCoh);

toc  
save(fileName,'cohDat','specDat','data')


%% september 9, we will only look at 200 ms before movement start and movement
%start

% changed taper frequency to [0.1 10] becuase our time period is 200 ms

fileName = '180328_specCoh_ECOG_gridMeanSsubtracted_SC32_200msAfterGo_100trials_Oct9';

tic
for ti = 1:7 %  7 targets
    
    
    eval(sprintf('dataTemp = dataTarget.target%d;',ti))
    dataCoh(1).tsDat_tr = permute(dataTemp,[2 3 1]);
    dataCoh(1).Fs = Fs_lfp;
    
    eval(sprintf('dataTemp = dataTarget_SC32.target%d;',ti))
    dataCoh(2).tsDat_tr = permute(dataTemp,[2 3 1]);
    dataCoh(2).Fs = Fs_lfp;
    [cohDat, specDat, ~] = spontCoh('', [], '', cohParams, datParams, dataCoh);
    
    cohDat_target{ti} = cohDat;
    
end

disp('Finished')
save(fileName,'cohDat_target','specDat')
toc

%% after that we are going to see some interesting patterns 
    freqRange = [0 4];

for ti = 1:7
    
    cohDat_temp = cohDat_target{ti};
    
    cohSeq_ECOG= cohDat_temp(1).coh; %1 for ECOG32
    cohSeq_SC32 = cohDat_temp(2).coh; %2 for SC32
    cohSeq_ECOG_SC32 = cohDat_temp(3).coh; %3 interDrive spectral coherency
    

    freqRangeIndex  = cal_index_freq(specDat(1).freq ,freqRange(1),freqRange(2));
    
    cohSeq_fr_ECOG_target(ti,:,:) = mean(cohSeq_ECOG(:,:,freqRangeIndex),3);
    cohSeq_fr_SC32_target(ti,:,:) = mean(cohSeq_SC32(:,:,freqRangeIndex),3);
    cohSeq_fr_ECOG_SC32_target(ti,:,:) = mean(cohSeq_ECOG_SC32(:,:,freqRangeIndex),3);

end

%% plot out the matrix

%plotted 
figure
for ti = 1:7
    subplot(2,4,ti)
    imagesc(squeeze(cohSeq_fr_ECOG_target(ti,:,:)))
    title(sprintf('Target %d',ti))
end 

%% calculate out the image display 


perIndex = 99;
fprintf('Percentile used to construct the map: %d\n',perIndex)

figure    
for ti = 1:7
    
    cohSeq_fr_temp = squeeze(cohSeq_fr_ECOG_target(ti,:,:));
    
    [nRow,nCol] = size(cohSeq_fr_temp);
    cohSeq_fr_temp = reshape(cohSeq_fr_temp,nRow*nCol,1);
    
    coh_th_temp = prctile(cohSeq_fr_temp, perIndex );
    
    
    
    graphParam.threshold = coh_th_temp;
    graphParam.thresholdPercentile = perIndex;
    graphParam.driveName = 'ECOG';
    graphParam.drive_idx =  1; %1 for ECOG; 2 for SC32
    cMatrix =  squeeze(cohSeq_fr_ECOG_target(ti,:,:));
    
    graphParam.textShift = 0.1;
    P_chamber = P_chamber_ECOG_3;
    
    graphParam.freqRange = freqRange;
    graphParam.labelElectrode = true;
    graphParam.textShift = 0.1;
    graphParam.acList = [experiment.hardware.microdrive(graphParam.drive_idx).electrodes.acquisitionid];
    
    subplot(2,4,ti)
    mapCoherenceGraph(cMatrix,P_chamber, graphParam)
    
end

% this is not too useful, we need to the occurance probability matrix. 


%% calculate tf trial by trial basis

fileName = '180328_specCoh_ECOG_gridMeanSsubtracted_SC32_100trials_Sept10';

cohParams.tapers = [0.1 10] ;  
cohParams.sampling = 1000; %Hz
cohParams.dn = 0.05 ;% s

datParams.trMin = 1; % can we even do single trial coherence

coh_matrix_trials = [];
tic
for ti = 1:7%  7 targets
 
    eval(sprintf('data_target_temp = dataTarget.target%d;',ti));
    eval(sprintf('data_target_sc32_temp  = dataTarget_SC32.target%d;',ti));
    
    data_target_temp = cat(2, data_target_temp,data_target_sc32_temp);
    
    n_trials = size(data_target_temp,1 );
    
    for trial_i = 1:n_trials
        
        %load ECOG data
        dataTemp = squeeze(data_target_temp(trial_i,:,:));
        [coh_matrix,f] = cal_coh(dataTemp,cohParams);
        
        coh_matrix_trials(ti,trial_i,:,:,:) = coh_matrix;
        
    end
    
    fprintf('Finished %d\n',ti)
end
   
disp('Finished calculate trial by trial coherency')
save(fileName,'coh_matrix_trials','specDat')
toc

% the format of coh_matrix
% target * trial * electrode * electrode * freq
% we need to compress target and trial dimensions into one 
% we then need to average the frequencies

%%
% we write our own function
% let's test it out 

[m,~]= cal_coh(zeros(101,1),cohParams);

%% we don't need to sort the targets at all, we can redo this

fileName = '180328_specCoh_ECOG_gridMeanSsubtracted_SC32_100trials_imaginaryCOh_Sept17';

cohParams.tapers = [0.1 10] ;
cohParams.sampling = 1000; %Hz
cohParams.dn = 0.05 ;% s

datParams.trMin = 1; % can we even do single trial coherence
coh_matrix_trials = [];
data_target_temp = cat(2, data_ref_mean, data_twoDrives{1,2});
n_trials = size(data_target_temp,1 );
     
tic
for trial_i = 1:n_trials
    
    %load ECOG data
    dataTemp = squeeze(data_target_temp(trial_i,:,:));
    [coh_matrix,f] = cal_coh(dataTemp,cohParams);
    
    coh_matrix_trials(trial_i,:,:,:) = coh_matrix;
    
    fprintf('Finished %d of %d\n',trial_i,n_trials)
    
end            

disp('Finished calculate trial by trial coherency')
toc

% coh_matrix_trials: trial * electrode * electrode * freq
% the target info is stored in trial_TargetAssigned

save(fileName, 'coh_matrix_trials', 'f', '-v7.3')
%% compress the frequency information
% need to calculate the frequency information for different bands

freqRange_all = [0 4;...
    4 13;...
    13 30;...
    30 80;...
    80 300];
freq_n = size(freqRange_all,1);

for fi = 1:freq_n
    %calculate the frequency range
    freqRange = freqRange_all(fi,:); %does it make sense to do this?
    freqRangeIndex = cal_index_freq(f ,freqRange(1),freqRange(2));
    
    %do the averaging
    coh_matrix_trials_freq = ...
        squeeze(mean(coh_matrix_trials(:,:,:,freqRangeIndex),4));
    
    % rearrange the matrix into the form ## * ## * number of trials
    
    coh_matrix_trials_freq = permute(coh_matrix_trials_freq,[2 3 1]);
    
    %imagesc(sq(coh_matrix_trials_freq(:,:,1)))
    
    % shovel it to annova analysis
    
    feature_p_values_all = fs_annova(coh_matrix_trials_freq,trial_TargetAssigned);
    
    feature_p_values_all_frequency(fi,:,:) = feature_p_values_all;
    
    fprintf('Finished compressing the coh matrix to frequency range of %d to %d\n',freqRange(1),freqRange(2))
end


%% process the matrix for later use
fi = 5;
feature_p_values_all  = sq( feature_p_values_all_frequency(fi,:,:));
freqRange = sq(freqRange_all(fi,:)); %does it make sense to do this?


p_cutoff = 0.01;

% convert upper triangular to ones
mat_matrix_zero = zeros(size(feature_p_values_all));
mat_matrix_ones = ones(size(feature_p_values_all));

%we transpose feature_p_values  to upper triangular and append the lower
%diagnal with ones (ones are useless in p_values)
mat_matrix_temp = tril(feature_p_values_all,-1)' + tril(mat_matrix_ones);

%clean up the matrix
outlierElectrodes = [47,59,163];
lowered_electrodes = Trials(1).Depth{1,2};
lowered_electrodes = find(lowered_electrodes == 0)+ 211; % shift by 211

%attach to the outlier matrix
%this list is a list 
outlierElectrodes = [outlierElectrodes lowered_electrodes];

%set these outlier 
mat_matrix_temp(outlierElectrodes,:) = 1;
mat_matrix_temp(:,outlierElectrodes) = 1;

feature_p_values_all_proc = mat_matrix_temp;

figure
title('Log Transformed P-value Connections')
imagesc(mat_matrix_temp.*(mat_matrix_temp < p_cutoff))
colorbar

length_mat = size(mat_matrix_temp,1);
mat_matrix_temp = mat_matrix_temp + mat_matrix_temp' - ones(length_mat,length_mat);

figure 
title('Number of Significant Connections per Electrode')
scatter(1:size(mat_matrix_temp,1),sum(mat_matrix_temp < p_cutoff))
xlabel('Electrodes (ECOG + SC32)')
ylabel('Number of Significant Connections including InterDrive')


%% first thing first, statistical distribution

mat_matrix_temp = feature_p_values_all;

[nRow,nCol] = size(mat_matrix_temp);
feature_p_values_all_column = reshape(mat_matrix_temp,nRow*nCol,1);

%set all zero elements to Na
feature_p_values_all_column(feature_p_values_all_column == 1) = NaN;

%normalize to thet total number of counts
figure
histogram(feature_p_values_all_column, 'Normalization','pdf')
xlabel('P-values')
ylabel('pdf')

%% we now need to plot this 

%% we need to look at occurance probability to estimate the robustness of the connections
index_range = 1:211;

cMatrix = -log10(feature_p_values_all_proc(index_range,index_range));
graphParam.threshold = 2;

graphParam.thresholdPercentile = 0.01;
graphParam.driveName = 'ECOG';
graphParam.drive_idx =  1; %1 for ECOG; 2 for SC32

graphParam.textShift = 0.1;
P_chamber = P_chamber_ECOG_3;

graphParam.freqRange = freqRange;
graphParam.labelElectrode = true;
graphParam.textShift = 0.1;
graphParam.acList = [experiment.hardware.microdrive(graphParam.drive_idx).electrodes.acquisitionid];

figure
subplot(1,3,1)
mapCoherenceGraph(cMatrix,P_chamber, graphParam)
title('ECOG')
% what this shows that is that by looking at instantaneous synchronicity
% between electrodes, we can statistically differentiate between targets


%now SC32
index_range = 212:243;
cMatrix = -log10(feature_p_values_all_proc(index_range,index_range));
P_chamber = P_chamber_SC32;
graphParam.drive_idx =  2; %1 for ECOG; 2 for SC32
graphParam.acList = [experiment.hardware.microdrive(graphParam.drive_idx).electrodes.acquisitionid];

subplot(1,3,2)
mapCoherenceGraph(cMatrix,P_chamber, graphParam)
title('SC32')

%now interdrive
index_range_1 = 1:211;
index_range_2 = 212:243;

%need to transpose the matrix
cMatrix = -log10(feature_p_values_all_proc(index_range_1,index_range_2)');

graphParam.driveName1 = 'SC32';
graphParam.P_chamber1 = P_chamber_SC32;

graphParam.driveName2 = 'ECOG';
graphParam.P_chamber2 = P_chamber_ECOG_3;

graphParam.drive_idx =  2; %1 for ECOG
graphParam.textShift = 0.1;

subplot(1,3,3)
A = cMatrix > graphParam.threshold;
inter_gplot(A,graphParam)

%% afterwards, 






      







    