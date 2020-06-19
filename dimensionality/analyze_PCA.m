% load some example recordings
global MONKEYDIR
%MONKEYDIR = 'D:\SiJia\OneDrive - UW\projects\Brain EEG\data';
MONKEYDIR = 'C:\Users\Si Jia\OneDrive - UW\projects\Brain EEG\data'
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

trialInfo.tBefore = -4e2; %ms
trialInfo.tAfter  = 2e2;
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
selectedfile = 'C:\Users\Si Jia\OneDrive - UW\projects\Brain EEG\code\posTargets_180328.mat';
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

%% get an example spectrum

example_electrodes = [1,11,218,211+13];

ei_d = 218
%set up parameters
params.Fs = 1000;
tfproduct = 20;
params.tapers = [tfproduct 2*tfproduct-1];
params.trialave = 0;
params.fpass = [0 200];
movingwin = [0.2  0.05];
num_ticks = 4;

dat = squeeze(trLfpData(:,ei_d,:));
data = dat'; % times * trials

tf_products = [5 10 20 40];

for i = 1:length(tf_products)
    
    tfproduct = tf_products(i)
    params.tapers = [tfproduct 2*tfproduct-1];
    
    %calculate the  spectrum
    [S,t,f] = mtspecgramc( data, movingwin, params );
    t = t + trialInfo.tBefore/1000; %shift relative to movement start
    S_Z = zscore(S,[],1);
    
    % run PCA analysis on the S_Z
    target_id = 1:7;
    %S_Z_target = S_Z(:,:,trial_TargetAssigned == target_id);
    
    [num_time, num_freq, num_trials] = size(S_Z);
    
    S_Z_flat = reshape(S_Z, num_time*num_freq, num_trials);
    
    [coeff,score,latent,tsquared,explained,mu] = pca(S_Z_flat'); %pca runs on n by p, num of trials by p features
    
    hold on
    plot(cumsum(explained))
    xlabel('Number of components')
    ylabel('Cumulative variance')
    
    
end
legend({'5','10','20','40'})
saveas(gca, 'pca_cum_sum.png')

%%

tfproduct = 40
params.tapers = [tfproduct 2*tfproduct-1];

%calculate the  spectrum
[S,t,f] = mtspecgramc( data, movingwin, params );
t = t + trialInfo.tBefore/1000; %shift relative to movement start
S_Z = zscore(S,[],1);


[num_time, num_freq, num_trials] = size(S_Z);

S_Z_flat = reshape(S_Z, num_time*num_freq, num_trials);

[coeff,score,latent,tsquared,explained,mu] = pca(S_Z_flat'); %pca runs on n by p, num of trials by p features


figure('units','normalized','outerposition',[0 0 1 1])
for num_pca = 1:6
    
    
    S_Z_PCA_comp = squeeze(mean(coeff(:,num_pca),2));
    S_Z_PCA_comp = reshape(S_Z_PCA_comp,num_time, num_freq);
    
    ax = subplot(2,3,num_pca)
    imagesc(S_Z_PCA_comp')
    axis xy
    colorbar
    
    if ei_d > 211
        ei_d_p = ei_d - 211;
        title(['PC ',num2str(num_pca) ,' for SC32 ', num2str(ei_d_p)])
    else
        title(['PC ',num2str(num_pca) ,' for ECoG ', num2str(ei_d)])
    end
    
    colorbar
    axis xy
    add_t_f_axes(ax,t,f,trialInfo,num_ticks)
    
end

saveas(gca, ['pca_time_frequency_tfproduct',num2str(tfproduct),'.png'])



%% update the figures with behaviour info

