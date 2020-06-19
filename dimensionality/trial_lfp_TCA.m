% load some example recordings
% this does not need to change
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

trialInfo.tBefore = -200; %ms
trialInfo.tAfter  = 200;
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

%% prep TCA analysis

%get rid of the corrupted electrodes 
%either noisy electrodes or lowered electrodes

tca_data = trLfpData;
tca_data(:,trialInfo.badE,:) = [];

%we feed this into TCA analysis and try to identify the rank
%we annotate thi data
%num of trials by electrodes by time
options.MinR = 1;
options.MaxR = 12;%you have to supply MinR and MaxR at the same time
rankest(tca_data, options);
%%
R = 3;
[Uhat, output] = cpd(tca_data,R);

figure
subplot(131)
plot(Uhat{1})
title('Trial')
xlabel('Trial number')

subplot(132)
plot(Uhat{2})
title('Electrode')
xlabel('Electrode number')
subplot(133)
plot(trialInfo.timeVec,Uhat{3})
title('Time')
xlabel('ms')

saveas(gca,'TCA_decomp.png')

%% this is not going anywhere
% can we try only looking at one drive
tca_data = trLfpData(:,trialInfo.goodECoGs,:);


%we feed this into TCA analysis and try to identify the rank
%we annotate thi data
%num of trials by electrodes by time
options.MinR = 1;
options.MaxR = 12;%you have to supply MinR and MaxR at the same time
rankest(tca_data, options);

saveas(gca, 'ECoG_error')

%%
R = 12;
[Uhat, output] = cpd(tca_data,R,'Display',1);

figure
subplot(131)
plot(Uhat{1})
title('Trial')
xlabel('Trial number')

subplot(132)
plot(Uhat{2})
title('Electrode')
xlabel('Electrode number')
subplot(133)
plot(trialInfo.timeVec,Uhat{3})
title('Time')
xlabel('ms')

saveas(gca,['TCA_ECoG_decomp_rank_',num2str(R),'.png'])


%% try this only on SC32

% can we try only looking at one drive
tca_data = trLfpData(:,trialInfo.goodSC32,:);


%we feed this into TCA analysis and try to identify the rank
%we annotate thi data
%num of trials by electrodes by time
options.MinR = 1;
options.MaxR = 12;%you have to supply MinR and MaxR at the same time
rankest(tca_data, options);

saveas(gca, 'SC32_error.png')
%this analysis cannot find the L corner that indicates the dominant modes. 
% partial reasons are: 1. the time series are not smoothed enough, maybe
% 2. in Williams' paper, they used non-negative TCA. we don't have that
% ability here. we 
%while, it does not work like magic that we expected, but it does show the
%noisiness of the signa., similar to what we find with PCA. 

%%

R = 12;
[Uhat, output] = cpd(tca_data,R,'Display',1);

figure
subplot(131)
plot(Uhat{1})
title('Trial')
xlabel('Trial number')

subplot(132)
plot(Uhat{2})
title('Electrode')
xlabel('Electrode number')
subplot(133)
plot(trialInfo.timeVec,Uhat{3})
title('Time')
xlabel('ms')

saveas(gca,['TCA_SC32_decomp_rank_',num2str(R),'.png'])

%% single electrode analysis

ei_d = 218
%set up parameters
params.Fs = 1000;
tfproduct = 5;
params.tapers = [tfproduct 2*tfproduct-1];
params.trialave = 0;
params.fpass = [0 200];
movingwin = [0.2  0.05];
num_ticks = 4;

dat = squeeze(trLfpData(:,ei_d,:));
data = dat'; % times * trials
%calculate the  spectrum
[S,t,f] = mtspecgramc( data, movingwin, params );
t = t + trialInfo.tBefore/1000; %shift relative to movement start
S_Z = zscore(S,[],1);

%we will first test out the 
options.MinR = 1;
options.MaxR = 5;%you have to supply MinR and MaxR at the same time
rankest(tca_data, options);

saveas(gca, 'ranktest_tfspec_electrode',num2str(ei_d),'.png')

R = 5;
[Uhat, output] = cpd(S_Z,R,'Display',1);

figure
subplot(131)
plot(t,Uhat{1})
title('Time')
xlabel('ms')

subplot(132)
plot(f,Uhat{2})
title('Frequency')
xlabel('Hz')

subplot(133)
plot(Uhat{3})
title('Trial')
xlabel('Trial number')

saveas(gca,['TCA_ECoG_tfspe_decomp_electrode_', num2str(ei_d),'_rank_',num2str(R),'.png'])

%we 


