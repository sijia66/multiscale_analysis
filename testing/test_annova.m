%%
% load some example recordings
global MONKEYDIR
MONKEYDIR = 'D:\SiJia\OneDrive - UW\projects\Brain EEG\data';
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

trialInfo.tBefore = -2e3;
trialInfo.tAfter  = 1e3;
trialInfo.timeVec = linspace(trialInfo.tBefore,trialInfo.tAfter, trialInfo.tAfter - trialInfo.tBefore);
trialInfo.timeReachStart = abs(trialInfo.tBefore);
trialInfo.trig = 'ReachStart';
trialInfo.lfpType = 'lfp';
trialInfo.Fs_lfp = 1000;
trialInfo.badECoGs = [47 59 163];
trialInfo.goodECoGs = setdiff(1:211,trialInfo.badECoGs);


%we are going to look at all electrode,  211 ECOG and 32 SC32 electrodes
trialInfo.proc_electrodes = 1:243;

%use constant target defination
selectedfile = 'D:\SiJia\OneDrive - UW\projects\Brain EEG\code\posTargets_180328.mat';
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

%% build the feature set
% test which frequency band could carry more information
disp('test which frequency band could carry more information')
example_electrode = 218;

dat = squeeze(trLfpData(:,example_electrode,:));
data = dat'; % times * trials

%set up parameters
params.Fs = 1000;
tfproduct = 5;
params.tapers = [tfproduct 2*tfproduct-1];
params.trialave = 0;
params.fpass = [0 200];
movingwin = [0.2  0.05];

num_ticks = 8;


[S,t,f] = mtspecgramc( data, movingwin, params );
feature_p_values = fs_annova(S,trial_TargetAssigned);

figure
imagesc(-log10(feature_p_values'))
colorbar
axis xy
%add time points
t_length = length(t);
f_length = length(f);

%add movement start
xline(round(t_length * trialInfo.timeReachStart / (trialInfo.timeReachStart + trialInfo.tAfter)), 'LineWidth',2);

s_timePoints = t;
freq = f;
s_timePoints = s_timePoints + trialInfo.tBefore / 1000;

tick_step = round(t_length / num_ticks);
xticks(1:tick_step:t_length)
xticklabels(s_timePoints(1:tick_step:end))

tick_step = round(f_length / num_ticks);
yticks(1:tick_step:f_length)
yticklabels(round(freq(1:tick_step:f_length)))

saveas(gca, ['feature_map_p_values_',num2str(example_electrode), '.png'])


%% test if z-score normalization would change the classification accuracy

disp('test if z score normalization would increase freq contents')
example_electrode = 218;

dat = squeeze(trLfpData(:,example_electrode,:));
data = dat'; % times * trials

%set up parameters
params.Fs = 1000;
tfproduct = 5;
params.tapers = [tfproduct 2*tfproduct-1];
params.trialave = 0;
params.fpass = [0 200];
movingwin = [0.2  0.05];

num_ticks = 8;


[S,t,f] = mtspecgramc( data, movingwin, params );

S_Z = zscore(S,[],1);

feature_p_values = fs_annova(S_Z,trial_TargetAssigned);

figure
imagesc(-log10(feature_p_values'))
colorbar
axis xy
title(['feature map p_values zscore ',num2str(example_electrode)])
%add time points
t_length = length(t);
f_length = length(f);

%add movement start
xline(round(t_length * trialInfo.timeReachStart / (trialInfo.timeReachStart + trialInfo.tAfter)), 'LineWidth',2);

s_timePoints = t;
freq = f;
s_timePoints = s_timePoints + trialInfo.tBefore / 1000;

tick_step = round(t_length / num_ticks);
xticks(1:tick_step:t_length)
xticklabels(s_timePoints(1:tick_step:end))

tick_step = round(f_length / num_ticks);
yticks(1:tick_step:f_length)
yticklabels(round(freq(1:tick_step:f_length)))

saveas(gca, ['feature_map_p_values_zscore_',num2str(example_electrode), '.png'])

figure
imagesc(squeeze(mean(S_Z,3))') %average across trials
colorbar
axis xy
title(['spec map zscore mean e',num2str(example_electrode)])
%add time points
t_length = length(t);
f_length = length(f);

%add movement start
xline(round(t_length * trialInfo.timeReachStart / (trialInfo.timeReachStart + trialInfo.tAfter)), 'LineWidth',2);

s_timePoints = t;
freq = f;
s_timePoints = s_timePoints + trialInfo.tBefore / 1000;

tick_step = round(t_length / num_ticks);
xticks(1:tick_step:t_length)
xticklabels(s_timePoints(1:tick_step:end))

tick_step = round(f_length / num_ticks);
yticks(1:tick_step:f_length)
yticklabels(round(freq(1:tick_step:f_length)))

saveas(gca, ['spec_map_zscore_mean_e',num2str(example_electrode), '.png'])

