% load some example recordings
global MONKEYDIR
MONKEYDIR = 'D:\SiJia\OneDrive - UW\projects\Brain EEG\data';
%MONKEYDIR = 'C:\Users\Si Jia\OneDrive - UW\projects\Brain EEG\data'
selectedfile = 'D:\SiJia\OneDrive - UW\projects\Brain EEG\data\posTargets_180328.mat';
addpath(genpath([MONKEYDIR '/m/']))

% look at the sessions
drive_base = 'LM1_ECOG';
driveSessions = makeDriveDatabase(drive_base,{'180328','180328'});
Sessions = cat(1, driveSessions{:});
driveNames = Sessions(:,3);
driveNameAnalyze = {'LM1_ECOG_3'};
useSess          = ismember(driveNames, driveNameAnalyze);
DayAnalyze = Sessions(useSess,1);

% define processing information
trialInfo.sequence_random = 1;

%yes, we are filtering out the long acq time
trialInfo.filterAcq = 1;
trialInfo.filterAcq_time = 5000; %ms

trialInfo.tBefore = 0; %ms
trialInfo.tAfter  = 6e2;
trialInfo.timeVec = linspace(trialInfo.tBefore,trialInfo.tAfter, trialInfo.tAfter - trialInfo.tBefore);
trialInfo.timeReachStart = abs(trialInfo.tBefore);
trialInfo.trig = 'ReachStart';
trialInfo.lfpType = 'lfp';
trialInfo.Fs_lfp = 1000;
trialInfo.badECoGs = [47 59 163];
trialInfo.ECoG_offset = 211;
trialInfo.goodECoGs = setdiff(1:211,trialInfo.badECoGs);


%we are going to look at all electrode,  211 ECOG and 32 SC32 electrodes
trialInfo.proc_electrodes = 1:243;

%use constant target defination
%selectedfile = 'C:\Users\Si Jia\OneDrive - UW\projects\Brain EEG\code\posTargets_180328.mat';
load(selectedfile)

%we load the regular trials
SessAnalyze = Sessions(useSess,:);
nSess = size(SessAnalyze,1);

for iD=1:nSess
    %load the trials and get the depth information
    trFN = [MONKEYDIR '/' DayAnalyze{iD} '/mat/Trials.mat'];
    load(trFN,'Trials')
    trialInfo.sessName =  DayAnalyze{iD};
    trialInfo.depthProfile = (Trials(1).Depth{1,2})'; % in micron
    trialInfo.goodSC32 = (find( trialInfo.depthProfile > 0))'; 
    trialInfo.goodE = [trialInfo.goodECoGs trialInfo.goodSC32+trialInfo.ECoG_offset];
    trialInfo.badE = setdiff(trialInfo.proc_electrodes,trialInfo.goodE);
    
    %filter out trials
    [trialInfo,Trials] = filter_trials(trialInfo,Trials);
    
    %first load an experiment trial
    expFN = [MONKEYDIR '/' DayAnalyze{iD} '/' Trials(1).Rec '/rec' Trials(1).Rec '.experiment.mat'];
    load(expFN,'experiment')
    
    %load data
    [trLfpData,trialInfo] = load_data(Trials, experiment,trialInfo,MONKEYDIR);
    
    trLfpData(:,trialInfo.badECoGs,:) = nan;
    
        %do target sorting
    [Pos, targids] = calcSeqReachTarg(Trials);
    Pos_Seq_1 = squeeze(Pos(:,2,:));
    trial_TargetAssigned = assignTaskNumber(Pos_Seq_1_unique,Pos_Seq_1);
end


%% load visual evoked potentials

%setting basic info
%clear all

addpath(genpath([MONKEYDIR '/m/']))
drive_base = 'LM1_ECOG';

%trial-alignment parameters          
bn = [trialInfo.tBefore trialInfo.tAfter];
trig = 'StartOn';
lfpType = 'lfp';

%get visEP 'sessions'
taskSessions = VisEP_Database;

%get sessions for 'drive' of interest
driveSessions = makeDriveDatabase(drive_base,{'180328','180328'});
%merge to get sessions w/ labview, split by drive
Sessions = mergeTaskDriveSessions(taskSessions, driveSessions);

%convert to a cell (#sessions x N) rather than nested cells b/c easier  to
%work wiith
Sessions = cat(1, Sessions{:});

%get driveName and theta
driveNames = Sessions(:,3);
driveTheta = [Sessions{:,7}];

%find unique driveNames and the theta for each one
[unique_driveNames, inds] = unique(driveNames);
unique_driveTheta = driveTheta(inds);
nDrives = length(unique_driveNames);

%sub-sampling of sessions if desired
%driveNameAnalyze = unique_driveNames; %all data
driveNameAnalyze = {'LM1_ECOG_3'};
useSess          = ismember(driveNames, driveNameAnalyze);

SessAnalyze       = Sessions(useSess,:);
DayAnalyze        = Sessions(useSess,1);
DriveNamesAnalyze = driveNames(useSess);

trialInfo.loweredE_SC32 = 1; 
trialInfo.exampleSC32 = 6;

%load all relevant trials
nSess = size(SessAnalyze,1);
eTrs = [];

for iD=1:nSess
    %we load the regular trials
    trFN = [MONKEYDIR '/' DayAnalyze{iD} '/mat/Trials.mat'];
    load(trFN,'Trials')
    
    %find out the loweredElectrodes
    if trialInfo.loweredE_SC32
        trialInfo.Ndrives = size(Trials(1).Depth,2);
        trialInfo.depthProfile = (Trials(1).Depth{1,2})'; % in micron
        var_temp = 1:length(trialInfo.depthProfile);
        var_temp(trialInfo.depthProfile > 0) = [];
        trialInfo.badSC32 = var_temp;
    else
        trialInfo.badSC32 = [];
    end
    trialInfo.goodSC32 = setdiff(1:32,trialInfo.badSC32);
    
    
    %load error trials
    eTrFN = [MONKEYDIR '/' DayAnalyze{iD} '/mat/ErrorTrials.mat'];
    load(eTrFN, 'ErrorTrials');
    sess_trial_size(nSess) = size(ErrorTrials,2);
    
    %visEPs have taskcode = 1 (touch task used)
    visStimTrials = ErrorTrials( getTaskCode(ErrorTrials)==1 );
    
    %patch work to fix the channel ID problem
    if ~isfield(visStimTrials,'ChannelID')
        for i = 1:length(visStimTrials)
            visStimTrials(i).ChannelID(1) = Trials(1).ChannelID(1);
            visStimTrials(i).ChannelID(2) = Trials(1).ChannelID(2);
        end
        disp('Filled the channelID structure')
    end
   
    dayTrials = visStimTrials;

    
    %load all trial data + position info of electrodes
    [trLfpData_ECOG32, eIDs, t, Pchamber_byDay{iD}] = loadTrLFP_session(dayTrials, SessAnalyze(iD,:), 'LM1_ECOG_3', trig, bn, lfpType);
    %load SC32 trials
    [trLfpData_SC32, eIDs, t, Pchamber_byDay{iD}] = loadTrLFP_session(dayTrials, SessAnalyze(iD,:), 'LM1_SC32_1', trig, bn, lfpType);
    
    trialInfo.N_ECOG = size(trLfpData_ECOG32,2);
    trialInfo.goodE = [1:trialInfo.N_ECOG trialInfo.goodSC32];
    trialInfo_all(iD) = trialInfo;
    
    %stack arrays together
    trLfpData_VESP = cat(2,trLfpData_ECOG32,trLfpData_SC32);
    
    %wrap up and clean up the variables
    clear 'trLfpData_ECOG32'
    clear 'trLfpData_SC32'
    clear 'trLfpData'
end

%% look at a trial across all electrodes

ti = 1;
data_temp = squeeze(trLfpData(ti, trialInfo.goodE,:));

[coeff,score,latent]  = pca(data_temp'); %the rows are observations, and the columns are variables

plot(cumsum(latent)/sum(latent))
xlabel('Number of Components')
ylabel('Variance')
colorbar

%% breakdown into ECoG and SC32 recordings

ti = 1;
data_temp_ECoG = squeeze(trLfpData(ti, trialInfo.goodECoGs,:));
data_temp_SC32 = squeeze(trLfpData(ti, ...
                    trialInfo.ECoG_offset+ trialInfo.goodSC32,:));

[coeff_ECoG,score_ECoG,latent_ECoG]  = pca(data_temp_ECoG'); %the rows are observations, and the columns are variables
[coeff_SC32,score_SC32,latent_SC32]  = pca(data_temp_SC32'); %the rows are observations, and the columns are variables


figure
plot(cumsum(latent_ECoG)/sum(latent_ECoG))
hold on
plot(cumsum(latent_SC32)/sum(latent_SC32))


xlabel('Number of Components')
ylabel('Fraction of variance')
axis([0 20 0 1])
legend({'ECoG','SC32'})

saveas(gca, ['PCA_ECoG_SC32_variance.png'])

%% next, the different frequency components
%https://www.mathworks.com/help/signal/ref/designfilt.html
freq_bands = [0 12;... %low frequency, alpha, delta
              12 30;...%beta 
              30 80;...%gammma
              80 150]; %high gamma band
filter_order = 20;
Fs = trialInfo.Fs_lfp;

%build a filter bank,  and promise no magic numbers
filter_lpFilt = designfilt('lowpassiir','FilterOrder',filter_order, ...
         'PassbandFrequency',freq_bands(1,2),'PassbandRipple',0.2, ...
         'SampleRate',Fs);
%beta band
filter_beta = designfilt('bandpassiir','FilterOrder',filter_order, ...
         'HalfPowerFrequency1',freq_bands(2,1),...
         'HalfPowerFrequency2',freq_bands(2,2), ...
         'SampleRate',Fs);

%gamma band
filter_gamma = designfilt('bandpassiir','FilterOrder',filter_order, ...
         'HalfPowerFrequency1',freq_bands(3,1),...
         'HalfPowerFrequency2',freq_bands(3,2), ...
         'SampleRate',Fs);

%high gamma band
filter_hi_gamma = designfilt('bandpassiir','FilterOrder',filter_order, ...
         'HalfPowerFrequency1',freq_bands(4,1),...
         'HalfPowerFrequency2',freq_bands(4,2), ...
         'SampleRate',Fs);

%save this to a bank 
filter_bank(1) = filter_lpFilt;
filter_bank(2) = filter_beta;
filter_bank(3) = filter_gamma;
filter_bank(4) = filter_hi_gamma;

%visualize the filter responses
for i = 1:length(filter_bank)
    fvtool(filter_bank(i))
end

%apply this to an example  recordings
ti = 1;
ei = 1;
data_temp_ECoG = squeeze(trLfpData(ti, ei,:));

figure
tVec = trialInfo.timeVec; 
subplot(1+length(filter_bank),1,1)
plot(tVec, data_temp_ECoG)
title('Original signal')

for i = 1:length(filter_bank)
    subplot(1+length(filter_bank),1,i+1)
    dataOut = filter(filter_bank(i),data_temp_ECoG );
    plot(tVec, dataOut)
    title([num2str(freq_bands(i,1)),...
           '-', ...
           num2str(freq_bands(i,2)),...
           ' Hz'
                   ])
end
xlabel('time relative to movement start (ms)')
saveas(gca, ['example electrode filtering.png'])
save('filter_bank.mat','filter_bank','freq_bands')

%% for the same recording, we examine its dimensionality at different frequency bands
%use the same frequency bands as above
%filter_bank 

ti = 1;
data_temp_ECoG = squeeze(trLfpData(ti, trialInfo.goodECoGs,:));
data_temp_SC32 = squeeze(trLfpData(ti, ...
                    trialInfo.ECoG_offset+ trialInfo.goodSC32,:));




figure
subplot(2,1,1)
[coeff_ECoG,score_ECoG,latent_ECoG]  = pca(data_temp_ECoG'); 
%the rows are observations, and the columns are variables



plot(cumsum(latent_ECoG)/sum(latent_ECoG),'LineWidth',2)
hold on


for i = 1:length(filter_bank)

    dataOut = filter(filter_bank(i),data_temp_ECoG' ); 
    dataOut = dataOut'; %same as the input dimension
    %2 for along the rows
    [coeff_ECoG,score_ECoG,latent_ECoG]  = pca(dataOut'); 
    plot(cumsum(latent_ECoG)/sum(latent_ECoG),'LineWidth',2)
    

end


ylabel('Fraction of variance', 'fontSize' , 14)
axis([0 20 0 1])
legend({'ECoG original 0 - 200 Hz)',...
        ['ECoG ',num2str(freq_bands(1,1)),'-',num2str(freq_bands(1,2)),' Hz'],...
        ['ECoG ',num2str(freq_bands(2,1)),'-',num2str(freq_bands(2,2)),' Hz'],...
        ['ECoG ',num2str(freq_bands(3,1)),'-',num2str(freq_bands(3,2)),' Hz'],...
        ['ECoG ',num2str(freq_bands(4,1)),'-',num2str(freq_bands(4,2)),' Hz'],...
        }, 'Location', 'best')
    

[coeff_SC32,score_SC32,latent_SC32]  = pca(data_temp_SC32'); 
subplot(2,1,2)
plot(cumsum(latent_SC32)/sum(latent_SC32),'LineWidth',2)

hold on
for i = 1:length(filter_bank)
    dataOut = filter(filter_bank(i),data_temp_SC32' ); 
    dataOut = dataOut'; %same as the input dimension
    %2 for along the rows
    [coeff_SC32,score_SC32,latent_SC32]  = pca(dataOut'); 
    plot(cumsum(latent_SC32)/sum(latent_SC32), 'LineWidth',2)

end
xlabel('Number of Components', 'fontSize' , 14)
ylabel('Fraction of variance', 'fontSize' , 14)
axis([0 20 0 1])
legend({'SC32 original 0 - 200 Hz)',...
        ['SC32 ',num2str(freq_bands(1,1)),'-',num2str(freq_bands(1,2)),' Hz'],...
        ['SC32 ',num2str(freq_bands(2,1)),'-',num2str(freq_bands(2,2)),' Hz'],...
        ['SC32 ',num2str(freq_bands(3,1)),'-',num2str(freq_bands(3,2)),' Hz'],...
        ['SC32 ',num2str(freq_bands(4,1)),'-',num2str(freq_bands(4,2)),' Hz'],...
        }, 'Location', 'best')

saveas(gca, 'pca_time_at_different_frequencies.png')

%%
%do this across different trials. 
dims_ECoG_freq = [];
dims_SC32_freq = [];

for ti = 1:length(Trials)
    data_temp_ECoG = squeeze(trLfpData(ti, trialInfo.goodECoGs,:));
    data_temp_SC32 = squeeze(trLfpData(ti, ...
        trialInfo.ECoG_offset+ trialInfo.goodSC32,:));
    
    [coeff_ECoG,score_ECoG,latent_ECoG]  = pca(data_temp_ECoG');
    [coeff_SC32,score_SC32,latent_SC32]  = pca(data_temp_SC32');
    dims_ECoG_freq(ti,1)= cal_eff_dim(latent_ECoG);
    dims_SC32_freq(ti,1)= cal_eff_dim(latent_SC32);
    
    for i = 1:length(filter_bank)
        dataOut = filter(filter_bank(i),data_temp_ECoG' );
        dataOut = dataOut'; %same as the input dimension
        [coeff_ECoG,score_ECoG,latent_ECoG]  = pca(dataOut');
        
        
        dataOut = filter(filter_bank(i),data_temp_SC32' );
        dataOut = dataOut'; %same as the input dimension
        %2 for along the rows
        [coeff_SC32,score_SC32,latent_SC32]  = pca(dataOut');
        
        dims_ECoG_freq(ti,i+1)= cal_eff_dim(latent_ECoG);
        dims_SC32_freq(ti,i+1)= cal_eff_dim(latent_SC32);
    end
end
disp('finished')

%%
%plot the hell of it
%creating grouping variables
C = dims_ECoG_freq(:,[1;1]*(1:size(dims_ECoG_freq,2)));
C(:,1:2:end) = dims_SC32_freq;
positions_temp  = (1:size(C,2)/2)';
positions_temp = [positions_temp positions_temp+0.25];
temp = positions_temp';
positions = temp(:);

figure
boxplot(C, ...
        'Positions', positions)
ylabel('Effective dimensionality')

set(gca,'xtick', mean(positions_temp'))
set(gca,'xticklabel',...
    {'Broadband',...
        [num2str(freq_bands(1,1)),'-',num2str(freq_bands(1,2)),' Hz'],...
        [num2str(freq_bands(2,1)),'-',num2str(freq_bands(2,2)),' Hz'],...
        [num2str(freq_bands(3,1)),'-',num2str(freq_bands(3,2)),' Hz'],...
        [num2str(freq_bands(4,1)),'-',num2str(freq_bands(4,2)),' Hz'],...
        })

color = ['b', 'y', 'b', 'y','b', 'y','b', 'y','b', 'y'];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');

hleg1 = legend(c(1:2), 'LFP', 'ECoG' );

set(gca,'FontSize',18, 'XTickLabelRotation',-30)

saveas(gca, 'pca_dim_freq_box.png')


%% do dim using different rereferencing techniques
% we compare the dimensionality with common average referencing
% our hypothesis is that common average referencing actually increases
% dimensionionality. whether that is useful is up to debate. 

ti = 1;
data_temp_ECoG = squeeze(trLfpData(ti, trialInfo.goodECoGs,:));
data_temp_SC32 = squeeze(trLfpData(ti, ...
                    trialInfo.ECoG_offset+ trialInfo.goodSC32,:));

% %do filtering
% dim = 3;
% [data_temp_ECoG_car, data_mean_ECoG]= ...
%     cal_car(data_temp_ECoG,dim);
% 
% [data_temp_SC32_car, data_mean_SC32]= ...
%     cal_car(data_temp_ECoG,dim);


%% break down into diff frequencies
dim = 1; % choose the 1st dim to apply car

dims_ECoG_freq = [];
dims_SC32_freq = [];
dims_ECoG_car_freq = [];
dims_SC32_car_freq = [];

for ti = 1:length(Trials)
    data_temp_ECoG = squeeze(trLfpData(ti, trialInfo.goodECoGs,:));
    data_temp_SC32 = squeeze(trLfpData(ti, ...
        trialInfo.ECoG_offset+ trialInfo.goodSC32,:));
    
    %add the car filters
    [data_temp_ECoG_car, data_mean_ECoG]= ...
                    cal_car(data_temp_ECoG,dim);
    [data_temp_SC32_car, data_mean_SC32]= ...
                    cal_car(data_temp_SC32,dim);
    
    
    [coeff_ECoG,score_ECoG,latent_ECoG]  = pca(data_temp_ECoG');
    [coeff_SC32,score_SC32,latent_SC32]  = pca(data_temp_SC32');
    [coeff_ECoG,score_ECoG,latent_ECoG_car]  = pca(data_temp_ECoG_car');
    [coeff_SC32,score_SC32,latent_SC32_car]  = pca(data_temp_SC32_car');
    
    dims_ECoG_freq(ti,1)= cal_eff_dim(latent_ECoG);
    dims_SC32_freq(ti,1)= cal_eff_dim(latent_SC32);
    dims_ECoG_car_freq(ti,1)= cal_eff_dim(latent_ECoG_car);
    dims_SC32_car_freq(ti,1)= cal_eff_dim(latent_SC32_car);
    
    
    for i = 1:length(filter_bank)
        dataOut = filter(filter_bank(i),data_temp_ECoG' );
        dataOut = dataOut'; %same as the input dimension
        [coeff_ECoG,score_ECoG,latent_ECoG]  = pca(dataOut');
        
        
        dataOut = filter(filter_bank(i),data_temp_SC32' );
        dataOut = dataOut'; %same as the input dimension
        %2 for along the rows
        [coeff_SC32,score_SC32,latent_SC32]  = pca(dataOut');
        
        %repeat this for the filters
        dataOut = filter(filter_bank(i),data_temp_ECoG_car' );
        dataOut = dataOut'; %same as the input dimension
        [coeff_ECoG,score_ECoG,latent_ECoG_car]  = pca(dataOut');
        
        
        dataOut = filter(filter_bank(i),data_temp_SC32_car' );
        dataOut = dataOut'; %same as the input dimension
        %2 for along the rows
        [coeff_SC32,score_SC32,latent_SC32_car]  = pca(dataOut');
        
        dims_ECoG_freq(ti,i+1)= cal_eff_dim(latent_ECoG);
        dims_SC32_freq(ti,i+1)= cal_eff_dim(latent_SC32);
        dims_ECoG_car_freq(ti,i+1)= cal_eff_dim(latent_ECoG_car);
        dims_SC32_car_freq(ti,i+1)= cal_eff_dim(latent_SC32_car);
    end
end
disp('finished')

%% plot the comparision for ECoG

C = dims_ECoG_freq(:,[1;1]*(1:size(dims_ECoG_freq,2)));
C(:,2:2:end) = dims_ECoG_car_freq;
positions_temp  = (1:size(C,2)/2)';
positions_temp = [positions_temp positions_temp+0.25];
temp = positions_temp';
positions = temp(:);

figure
boxplot(C, ...
        'Positions', positions)
ylabel('Effective dimensionality')

set(gca,'xtick', mean(positions_temp'))
set(gca,'xticklabel',...
    {'Original (0 - 200 Hz)',...
        [num2str(freq_bands(1,1)),'-',num2str(freq_bands(1,2)),' Hz'],...
        [num2str(freq_bands(2,1)),'-',num2str(freq_bands(2,2)),' Hz'],...
        [num2str(freq_bands(3,1)),'-',num2str(freq_bands(3,2)),' Hz'],...
        [num2str(freq_bands(4,1)),'-',num2str(freq_bands(4,2)),' Hz'],...
        })

color = ['b', 'y', 'b', 'y','b', 'y','b', 'y','b', 'y'];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');
hleg1 = legend(c(1:2), 'ECoG', 'ECoG-CAR' );
set(gca,'FontSize',10)
saveas(gca, 'pca_dim_ECoG_car_box.png')

%% comparing between SC32 and SC32-CAR

C = dims_SC32_freq(:,[1;1]*(1:size(dims_SC32_freq,2)));
C(:,2:2:end) = dims_SC32_car_freq;
positions_temp  = (1:size(C,2)/2)';
positions_temp = [positions_temp positions_temp+0.25];
temp = positions_temp';
positions = temp(:);

figure
boxplot(C, ...
        'Positions', positions)
ylabel('Effective dimensionality')

set(gca,'xtick', mean(positions_temp'))
set(gca,'xticklabel',...
    {'Original (0 - 200 Hz)',...
        [num2str(freq_bands(1,1)),'-',num2str(freq_bands(1,2)),' Hz'],...
        [num2str(freq_bands(2,1)),'-',num2str(freq_bands(2,2)),' Hz'],...
        [num2str(freq_bands(3,1)),'-',num2str(freq_bands(3,2)),' Hz'],...
        [num2str(freq_bands(4,1)),'-',num2str(freq_bands(4,2)),' Hz'],...
        })

color = ['b', 'y', 'b', 'y','b', 'y','b', 'y','b', 'y'];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');
hleg1 = legend(c(1:2), 'SC32', 'SC32-CAR' );
set(gca,'FontSize',10)
saveas(gca, 'pca_dim_SC32_car_box.png')


%% do dim reduction with multipolar and bipolar rereferencing
%get the chamber information
P_chamber_ECOG = drmap_coordTransform_microdrive2chamber(experiment,1);
P_chamber_ECOG = P_chamber_ECOG';

P_chamber_SC32 = drmap_coordTransform_microdrive2chamber(experiment,2);
P_chamber_SC32 = P_chamber_SC32';

%only look at good electrodes
P_chamber_ECoG = P_chamber_ECOG(trialInfo.goodECoGs,:);
P_chamber_SC32 = P_chamber_SC32(trialInfo.goodSC32,:);

%get common average referencing filters
refMatrix_car = ...
            car_matrix(P_chamber_ECoG);

%get laplacian matrix 
[refMatrix_lap, P_chamber_lap] = laplacian_matrix(P_chamber_ECoG,distThres);

%gen ref matrix
distThres = 0.75;
[refMatrix_multi, refLocation_multi] = ...
            multi_rereferencing_matrix(P_chamber_ECoG,distThres);

[refMatrix_bi, refLocation_bi] = ...
            bipolar_matrix(P_chamber_ECoG,distThres);

tca_data = trLfpData(:,trialInfo.goodECoGs,:);
X = tensor(tca_data);

%get car referencing
Y = ttm(X,refMatrix_car,2); %just only multiply slice section. 
trLfpData_ECoG_car = double(Y);

%get laplacian referencing
Y = ttm(X,refMatrix_lap,2); %just only multiply slice section. 
trLfpData_ECoG_lap = double(Y);

%get bipolar referencing
Y = ttm(X,refMatrix_bi,2); %just only multiply slice section. 
trLfpData_ECoG_bipolar = double(Y);

%get multipolar referencing
Y = ttm(X,refMatrix_multi,2); 
trLfpData_ECoG_multipolar = double(Y);

%% after the polar, we can actually do polar

for ti = 1:length(Trials)
    %grab the data
    data_temp_ECoG = squeeze(trLfpData(ti, trialInfo.goodECoGs,:));
%     data_temp_ECoG_car= squeeze(trLfpData_ECoG_car(ti, :,:));
%     data_temp_ECoG_bipolar = squeeze(trLfpData_ECoG_bipolar(ti, :,:));
%     data_temp_ECoG_multipolar = squeeze(trLfpData_ECoG_multipolar(ti, :,:));
        
    data_temp_ECoG_car= refMatrix_car * data_temp_ECoG;
    data_temp_ECoG_lap = refMatrix_lap * data_temp_ECoG;
    data_temp_ECoG_bipolar = refMatrix_bi * data_temp_ECoG;
    data_temp_ECoG_multipolar = refMatrix_multi * data_temp_ECoG;
        
  
    %cal dim
    [coeff_ECoG,score_ECoG,latent_ECoG]  = pca(data_temp_ECoG');
    [coeff_ECoG,score_ECoG,latent_ECoG_car]  = pca(data_temp_ECoG_car');
    [coeff_ECoG,score_ECoG,latent_ECoG_bipolar]  = pca(data_temp_ECoG_bipolar');
        [coeff_ECoG,score_ECoG,latent_ECoG_lap]  = pca(data_temp_ECoG_lap');
    [coeff_ECoG,score_ECoG,latent_ECoG_multipolar]  = pca(data_temp_ECoG_multipolar');
    
    dims_ECoG(ti,1)= cal_eff_dim(latent_ECoG);
    dims_ECoG_car(ti,1)= cal_eff_dim(latent_ECoG_car);
    dims_ECoG_bipolar(ti,1)= cal_eff_dim(latent_ECoG_bipolar);
    dims_ECoG_lap(ti,1)= cal_eff_dim(latent_ECoG_lap);
    dims_ECoG_multipolar(ti,1)= cal_eff_dim(latent_ECoG_multipolar);
end

figure
%then can do dim reduction 
boxplot([dims_ECoG dims_ECoG_car dims_ECoG_bipolar dims_ECoG_lap dims_ECoG_multipolar],...
        'Labels', {'Raw ECoG','car' ,'bipolar','Laplacian', 'multipolar'})
ylabel('Effectie dim')
xlabel('Common referencing techniques')
saveas(gca, 'comparing re-referencing technqiues.png')
%% do dim reduction with differnet resampling techniques
%% dimensionality of the visual evoked potentials
% get the signal from the previous section.
% do PCA on both ECoG and SC32
NUM_TRIALS = 121;

for ti = 1:length(Trials) %use the same number of trials as the reaching trials
    %grab the data
    data_temp_ECoG = squeeze(trLfpData_VESP(ti, trialInfo.goodECoGs,:));
    data_temp_SC32 = squeeze(trLfpData_VESP(ti, ...
                    trialInfo.ECoG_offset+ trialInfo.goodSC32,:));
   
    
    %cal dim
    [coeff_ECoG,score_ECoG,latent_ECoG]  = pca(data_temp_ECoG');
    [coeff_SC32,score_SC32,latent_SC32]  = pca(data_temp_SC32');
    
    dims_VESP_ECoG(ti,1)= cal_eff_dim(latent_ECoG);
    dims_VESP_SC32(ti,1)= cal_eff_dim(latent_SC32);

end




figure
plot(sq(mean(trLfpData_VESP(1:NUM_TRIALS,1,:))))
hold on 
plot(sq(mean(trLfpData_VESP(1:NUM_TRIALS,211+6,:))))
legend('ECoG E1', 'SC32 E6')

figure
%then can do dim reduction 
boxplot([dims_VESP_ECoG dims_VESP_SC32],...
        'Labels', {'Raw ECoG','Raw SC32'})
ylabel('Effectie dim')
xlabel('Types of Recordings')
saveas(gca, 'effective dims_ECoG_SC32.png')
%% dimensionality across visual evoked potentials

timeVec  = trialInfo.timeVec;
time_win = 100; %ms
time_start = -trialInfo.timeReachStart:...
    time_win:(max(timeVec) - time_win);
dims_ECoG = [];
dims_SC32 = [];
NUM_TRIALS = 121; % this has to match of the reaching trials

for ti = 1:NUM_TRIALS
    
    for i = 1:length(time_start)
        time_roi = find((timeVec >= time_start(i)) ...
            & (timeVec < (time_start(i) + time_win)));
        
        data_temp_ECoG = squeeze(trLfpData_VESP(ti, trialInfo.goodECoGs,time_roi));
        data_temp_SC32 = squeeze(trLfpData_VESP(ti, ...
            trialInfo.ECoG_offset+ trialInfo.goodSC32,time_roi));
        
        [coeff_ECoG,score_ECoG,latent_ECoG]  = pca(data_temp_ECoG'); %the rows are observations, and the columns are variables
        [coeff_SC32,score_SC32,latent_SC32]  = pca(data_temp_SC32'); %the rows are observations, and the columns are variables
        
        dims_ECoG(ti,i)= cal_eff_dim(latent_ECoG);
        dims_SC32(ti,i)= cal_eff_dim(latent_SC32);
        
    end
    
end

figure

p1 = boxplot( dims_ECoG, time_start / 1000);
hold on
p2 = boxplot( dims_SC32,time_start/1000);
plot( median(dims_ECoG), 'lineWidth', 3)
plot(median(dims_SC32), 'lineWidth', 3)
legend({'ECoG','SC32'})
xlabel('starting times relative to flash on (s)')
ylabel('effective dimensionality')
hold off
saveas(gca, ...
      ['effective_dimensionality_VESP_ECoG_SC32_',num2str(time_win),'.png'] )
%% dim across trials at different time points
%% look at variability across trials over different time points
timeVec  = trialInfo.timeVec;
time_win = 10; %ms
time_start = -600:...
    time_win:(max(timeVec) - time_win);
dims_ECoG_trials = [];
dims_SC32_trials = [];
NUM_TRIALS = 121;

for i = 1:length(time_start)
    time_roi = find((timeVec >= time_start(i)) ...
        & (timeVec < (time_start(i) + time_win)));
    
    %3 for averaging across time
    data_temp_ECoG = squeeze(mean(trLfpData_VESP(1:NUM_TRIALS, trialInfo.goodECoGs,time_roi),3));
    data_temp_SC32 = squeeze(mean(trLfpData_VESP(1:NUM_TRIALS, ...
        trialInfo.ECoG_offset+ trialInfo.goodSC32,time_roi),3));
    
    [coeff_ECoG,score_ECoG,latent_ECoG]  = pca(data_temp_ECoG); %the rows are observations, and the columns are variables
    [coeff_SC32,score_SC32,latent_SC32]  = pca(data_temp_SC32); %the rows are observations, and the columns are variables
    
    dims_ECoG_trials(i)= cal_eff_dim(latent_ECoG);
    dims_SC32_trials(i)= cal_eff_dim(latent_SC32);
    
end
figure
plot(time_start / 1000, dims_ECoG_trials)
hold on
plot(time_start / 1000, dims_SC32_trials)
xlabel('starting times (s)')
ylabel('Effective dimensionality')
legend('ECoG' , 'SC32')
saveas(gca, ['pca_dim_VESP_over_trials_timeWindow_',num2str(time_win),'.png'])
%% dimensionality relative to different behaviour times
timeVec  = trialInfo.timeVec;
time_win = 100; %ms
time_start = -trialInfo.timeReachStart:...
    time_win:(max(timeVec) - time_win);
dims_ECoG = [];
dims_SC32 = [];

for ti = 1:length(Trials)
    
    for i = 1:length(time_start)
        time_roi = find((timeVec >= time_start(i)) ...
            & (timeVec < (time_start(i) + time_win)));
        
        data_temp_ECoG = squeeze(trLfpData(ti, trialInfo.goodECoGs,time_roi));
        data_temp_SC32 = squeeze(trLfpData(ti, ...
            trialInfo.ECoG_offset+ trialInfo.goodSC32,time_roi));
        
        [coeff_ECoG,score_ECoG,latent_ECoG]  = pca(data_temp_ECoG'); %the rows are observations, and the columns are variables
        [coeff_SC32,score_SC32,latent_SC32]  = pca(data_temp_SC32'); %the rows are observations, and the columns are variables
        
        dims_ECoG(ti,i)= cal_eff_dim(latent_ECoG);
        dims_SC32(ti,i)= cal_eff_dim(latent_SC32);
        
    end
    
end

figure

p1 = boxplot( dims_ECoG, time_start / 1000);
hold on
p2 = boxplot( dims_SC32,time_start/1000);

plot( median(dims_ECoG), 'lineWidth', 3)
plot(median(dims_SC32), 'lineWidth', 3)
legend({'ECoG','SC32'})
xlabel('starting times relative to movement start (s)')
ylabel('effective dimensionality')

  hold off
saveas(gca, ...
      ['effective_dimensionality_ECoG_SC32_',num2str(time_win),'.png'] )
%% look at variability across trials over different time points
timeVec  = trialInfo.timeVec;
time_win = 10; %ms
time_start = -600:...
    time_win:(max(timeVec) - time_win);
dims_ECoG_trials = [];
dims_SC32_trials = [];

for i = 1:length(time_start)
    time_roi = find((timeVec >= time_start(i)) ...
        & (timeVec < (time_start(i) + time_win)));
    
    %3 for averaging across time
    data_temp_ECoG = squeeze(mean(trLfpData(:, trialInfo.goodECoGs,time_roi),3));
    data_temp_SC32 = squeeze(mean(trLfpData(:, ...
        trialInfo.ECoG_offset+ trialInfo.goodSC32,time_roi),3));
    
    [coeff_ECoG,score_ECoG,latent_ECoG]  = pca(data_temp_ECoG); %the rows are observations, and the columns are variables
    [coeff_SC32,score_SC32,latent_SC32]  = pca(data_temp_SC32); %the rows are observations, and the columns are variables
    
    dims_ECoG_trials(i)= cal_eff_dim(latent_ECoG);
    dims_SC32_trials(i)= cal_eff_dim(latent_SC32);
    
end
figure
plot(time_start / 1000, dims_ECoG_trials)
hold on
plot(time_start / 1000, dims_SC32_trials)
xlabel('starting times (s)')
ylabel('Effective dimensionality')
legend('ECoG' , 'SC32')
saveas(gca, ['pca_dim_over_trials_timeWindow_',num2str(time_win),'.png'])

%% cal_eff_dim functions
function dim_eff = cal_eff_dim(latent)
lamda_j = latent / sum(latent);
dim_eff = 1 / sum(lamda_j.^2);
end



