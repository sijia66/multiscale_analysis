%%
% the goal of this
%develope the code pipeline to process lfp data
global MONKEYDIR
MONKEYDIR = 'D:\SiJia\OneDrive - UW\projects\Brain EEG\data';
addpath(genpath([MONKEYDIR '/m/']))

% look at the sessions
drive_base = 'LM1_ECOG';
driveSessions = makeDriveDatabase(drive_base,{'180329','180329'});
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

trialInfo.tBefore = -1.5e3;
trialInfo.tAfter  = 1e3;
trialInfo.timeVec = linspace(trialInfo.tBefore,trialInfo.tAfter, trialInfo.tAfter - trialInfo.tBefore);
trialInfo.timeReachStart = abs(trialInfo.tBefore);
trialInfo.trig = 'ReachStart';
trialInfo.lfpType = 'lfp';
trialInfo.Fs_lfp = 1000;

%find the usage electrodes
trialInfo.badECoGs = [47 59 163];
trialInfo.goodECoGs = setdiff(1:211,trialInfo.badECoGs);
%trialInfo.goodSC32 = 211 + find(trialInfo.depthProfile > 0)

%define neural feature
featureROI  = [-50 200]; %in ms

%we are going to look at all electrode,  211 ECOG and 32 SC32 electrodes
trialInfo.proc_electrodes = 1:243;

%use constant target defination
selectedfile = 'D:\SiJia\OneDrive - UW\projects\Brain EEG\code\posTargets_180328.mat';
load(selectedfile)

%we load the regular trials
SessAnalyze = Sessions(useSess,:);
nSess = size(SessAnalyze,1);

%% then we compare the signal obtain to that our CAR

for iD=1:nSess
    %load the trials and get the depth information
    trFN = [MONKEYDIR '/' DayAnalyze{iD} '/mat/Trials.mat'];
    load(trFN,'Trials')
    trialInfo.sessName =  DayAnalyze{iD};
    trialInfo.depthProfile = (Trials(1).Depth{1,2})'; % in micron
    
    %filter out trials
    [trialInfo,Trials] = filter_trials(trialInfo,Trials);
    trialInfo.goodSC32 = 211 + find(trialInfo.depthProfile > 0)
    
    %first load an experiment trial
    expFN = [MONKEYDIR '/' DayAnalyze{iD} '/' Trials(1).Rec '/rec' Trials(1).Rec '.experiment.mat'];
    load(expFN,'experiment')
    
    %do target sorting
    [Pos, targids] = calcSeqReachTarg(Trials);
    Pos_Seq_1 = squeeze(Pos(:,2,:));
    trial_TargetAssigned = assignTaskNumber(Pos_Seq_1_unique,Pos_Seq_1);
    
    %load data
    [trLfpData,trialInfo] = load_data(Trials, experiment,trialInfo,MONKEYDIR);
    [trLfpData,trialInfo] = load_data(Trials, experiment,trialInfo,MONKEYDIR);
    [trLfpData_ECOG_CAR, data_mean_ECOG] = cal_car(trLfpData(:,trialInfo.goodECoGs,:),2);
    [trLfpData_SC32_CAR,data_mean_SC32] = cal_car(trLfpData(:,trialInfo.goodSC32,:),2);

end

%%
%look at sample electrodes
ei = 1;
ti = 1;

plot(trialInfo.timeVec,sq(trLfpData(ti,ei,:)))
hold on
plot(trialInfo.timeVec,sq(data_mean_ECOG(ti,1,:)))
plot(trialInfo.timeVec,sq(trLfpData_ECOG_CAR(ti,ei,:)))
xlabel('Time from movement start (ms)')
legend({'Raw','Electrode Mean','after common average filtering'})
title(['Electrode ' num2str(ti)])

figure

ei = 24;
ti = 1;

plot(trialInfo.timeVec,sq(trLfpData(ti,ei,:)))
hold on
plot(trialInfo.timeVec,sq(data_mean_ECOG(ti,1,:)))
plot(trialInfo.timeVec,sq(trLfpData_ECOG_CAR(ti,ei,:)))
xlabel('Time from movement start (ms)')
legend({'Raw','Electrode Mean','after common average filtering'})
title(['Electrode ' num2str(ti)])

